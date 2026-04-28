# 4. Integrators, NEE and MIS

Reference file: `src/integrators/path/path.cpp` (uni-directional path tracer
with NEE+MIS — this is the integrator you'll test against).

## 4.1 What the integrator hands to your BSDF

Per intersection, in a typical path-tracing iteration:

1. **Russian-roulette / max-depth** check.
2. **NEE — direct illumination from a sampled emitter:**
   ```cpp
   DirectSamplingRecord dRec(its);
   Spectrum value = scene->sampleEmitterDirect(dRec, sampler->next2D());
   BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);
   Spectrum f       = bsdf->eval(bRec);             // ← your eval()
   Float    bsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
                       ? bsdf->pdf(bRec) : 0;        // ← your pdf()
   Float w = miWeight(dRec.pdf, bsdfPdf);
   L += throughput * value * f * w;
   ```
3. **BSDF sampling — extends the path:**
   ```cpp
   BSDFSamplingRecord bRec(its, sampler, ERadiance);  // sampler! → uses extra rand
   Float bsdfPdf;
   Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, sampler->next2D());
   throughput *= bsdfWeight;
   eta *= bRec.eta;
   ```
4. Trace ray in direction `bRec.wo`; if it hits an emitter, MIS-weight the
   contribution against the NEE pdf:
   ```cpp
   Float lumPdf = !(bRec.sampledType & EDelta)
                  ? scene->pdfEmitterDirect(dRec) : 0;
   L += throughput * value * miWeight(bsdfPdf, lumPdf);
   ```

So your BSDF interacts with the integrator via three public faces:

| Face | Used for |
| --- | --- |
| `sample(bRec, pdf, u2)` | extending the camera path. Must return `f*cos/pdf`, set `bRec.wo`, `bRec.eta`, `bRec.sampledType`. |
| `eval(bRec, ESolidAngle)` | scoring the NEE direction. Must return `f*cos`. |
| `pdf(bRec, ESolidAngle)` | MIS weighting in both directions. Must return the same density that `sample()` would have produced. |

## 4.2 The MIS power heuristic in 2 lines

```cpp
inline Float miWeight(Float pdfA, Float pdfB) const {
    pdfA *= pdfA; pdfB *= pdfB;
    return pdfA / (pdfA + pdfB);
}
```

Standard β=2 power heuristic from Veach. Applied separately to each emitter
sample, so it's robust to wildly different PDF magnitudes.

For your LayeredBSDF this is what makes good `pdf()` *necessary*: a wrong
PDF is silently consumed as a wrong MIS weight and produces images that
look "fine but wrong" (correlated noise / fireflies / muted highlights).

## 4.3 What "compatible with NEE/MIS" means *concretely* for LayeredBSDF

* `eval()` must produce a value for **any** `(wi, wo)` pair the integrator
  hands you, including ones your random walk would not naturally produce.
  Your random-walk estimator therefore evaluates the BSDF as **a Monte
  Carlo estimator that treats `wo` as fixed**: you do the walk, and at
  each scatter you accumulate the contribution that exits in direction
  `wo` (a "next-event-style" closure inside the walk). This is exactly
  the position-free / 1D random walk approach in
  *Guo, Hašan, Yan 2018; Belcour 2018; Xia, Hašan, Tausonu, Bickel 2020*.

* `pdf()` must return the **density of `sample()`** at `wo`. Because your
  walk has many internal branches, an exact closed-form pdf is intractable
  — instead you use the same identity as `sample()` did to *define* the
  pdf:
  ```
  pdf(wo)  =  prob(top-Fresnel reflection)            * δ_R(wo)            // delta lobe (if exposed)
            + prob(transmit-and-walk)                  * pdf_walk(wo).
  ```
  `pdf_walk(wo)` is itself a single-sample MC estimate (you re-run the walk
  with one sample to estimate the marginal density of the exit direction).
  In practice you can use a **Newton/Beckmann-style approximation** of
  `pdf_walk` (an effective Lambert + cosine-power lobe whose parameters are
  fit to the layer optical depths) or, more accurately, the bidirectional
  MIS trick used by Guo et al. 2018 §4.4 ("position-free MC").

* Always set `bRec.eta` correctly on transmission. The path tracer uses it
  for refractive Russian-roulette throughput correction
  (`path.cpp:283`), so getting it wrong breaks unbiasedness in scenes that
  have IOR mismatches.

## 4.4 Other integrators (skip on first pass)

* `volpath.cpp` — same as `path.cpp` but supports participating media
  inside scene volumes. **Not needed** for your LayeredBSDF, because the
  layers live *inside* the BSDF — they're not full scene media.
* `bdpt`, `mlt`, `pssmlt` — bidirectional / Markov-chain methods. They
  require your BSDF to produce a valid PDF for the *adjoint* direction
  too; the `mode` field (`ERadiance`/`EImportance`) controls this.
  Keep this in mind but don't worry on first iteration.
* `direct.cpp` — same NEE/MIS logic but only depth=2. Fast smoke test:
  if your LayeredBSDF looks right under `direct`, you've got eval and
  NEE/MIS right.

## 4.5 Validation toolbox

* **`mtsutil testcases`** runs `src/tests/`. The `test_chisquare` testcase
  is the most useful for BSDFs — it draws thousands of `sample()`s,
  histograms them, and chi-square tests against `pdf()`. Add your plugin
  to `data/tests/test_bsdf.xml` and run:
  ```bash
  mtsutil testcases -t test_chisquare data/tests/test_bsdf.xml
  ```
* **White furnace test.** Render your material in a uniform white
  environment; the result must equal the directional-hemispherical
  reflectance × radiance (≈ 1 for an albedo-1 layer). Any energy gain or
  loss is a sign of an error in your eval/pdf.
* **Compare against `coating.cpp`** in the limit of a single dielectric
  layer with a Lambertian bottom — your renderer should match it (modulo
  the specular-reflection energy that `coating.cpp` admits but a true
  random walk also keeps).
* **Compare against `hk.cpp`** in the optically thin (single-scatter)
  limit, with the same `sigmaT`, `albedo`, `phase`, `thickness`.
