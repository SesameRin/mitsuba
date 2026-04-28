# 3. The Emitter system

Reference files: `include/mitsuba/render/emitter.h`,
`src/emitters/spot.cpp`, `src/emitters/directional.cpp`,
`src/emitters/point.cpp`, `src/emitters/area.cpp`.

## 3.1 Why the API has so many sampling methods

A BSDF has one job (scatter), but an emitter is queried by many algorithms:

| Algorithm | Needs to … |
| --- | --- |
| Path tracing (NEE) | Given a point `x` in the scene, sample a direction toward the emitter. → `sampleDirect` |
| Light tracing / BDPT | Sample a starting ray on the emitter. → `sampleRay` |
| Path tracing (BSDF hit) | Evaluate `Le(x, ω)` at a point that the BSDF random-walk happened to hit. → `evalDirection`, `evalPosition`, `pdfDirect` for MIS. |
| Photon mapping | Sample a ray and emit a photon. → `sampleRay` |

So the abstract base class `AbstractEmitter` (`emitter.h`) splits the
emission profile into a position component and a direction component, and
exposes both `sample`/`eval`/`pdf` for each, plus a combined `sampleRay`,
plus a `sampleDirect`/`pdfDirect` for "given a reference point, give me a
sample on the emitter that contributes to this point."

For your task you only need the **emitter side** (`Emitter` is the subclass
that adds an `isCompound`/`isEnvironmentEmitter`/etc. on top of
`AbstractEmitter`; sensors share the same base).

## 3.2 The methods you'll override

For a **delta-position, delta-direction** light (e.g. classical directional
light) the methods reduce to:

```cpp
Spectrum samplePosition(PositionSamplingRecord &pRec, const Point2 &sample, ...);
Spectrum sampleDirection(DirectionSamplingRecord &dRec,
                         PositionSamplingRecord &pRec,
                         const Point2 &sample, ...);
Spectrum sampleRay(Ray &ray,
                   const Point2 &spatialSample,
                   const Point2 &directionalSample,
                   Float time);
Spectrum sampleDirect(DirectSamplingRecord &dRec, const Point2 &sample);
Float    pdfDirect(const DirectSamplingRecord &dRec);

Spectrum evalPosition (const PositionSamplingRecord  &pRec);
Spectrum evalDirection(const DirectionSamplingRecord &dRec,
                       const PositionSamplingRecord  &pRec);
Float    pdfPosition  (const PositionSamplingRecord  &pRec);
Float    pdfDirection (const DirectionSamplingRecord &dRec,
                       const PositionSamplingRecord  &pRec);
```

The constructor must set `m_type` to a combination of
`EDeltaPosition`, `EDeltaDirection`, `EOnSurface` so the integrator knows
how to MIS against the emitter:

| Light | `m_type` |
| --- | --- |
| `point` | `EDeltaPosition` |
| `directional` | `EDeltaDirection` |
| `spot` | `EDeltaPosition` |
| `collimated` | `EDeltaPosition \| EDeltaDirection` |
| `area` | `EOnSurface` |
| `envmap` / `constant` | (neither — env emitter) |

## 3.3 Walk-through: `spot.cpp`

`spot.cpp` is a delta-position emitter (the source is at a single point),
with a smooth angular falloff over a cone:

```cpp
SpotEmitter(Properties &props) : Emitter(props) {
    m_intensity   = props.getSpectrum("intensity", Spectrum(1.0f));
    m_cutoffAngle = degToRad(props.getFloat("cutoffAngle", 20));
    m_beamWidth   = degToRad(props.getFloat("beamWidth", cutoff*0.75));
    m_type        = EDeltaPosition;        // <-- key flag
}
```

The actual emission profile is `falloffCurve(d)`, returning the spectral
intensity along local direction `d`. Then:

* `samplePosition` → returns the (single) emitter location, with `pdf=1` and
  `measure=EDiscrete`. The returned weight is `intensity * 4π` (because the
  internal convention factors out a `1/(4π)` later).
* `sampleDirection` → cosine-cone-sample inside `m_cosCutoffAngle`.
* `sampleDirect` → easy because the position is fixed: `dRec.p = origin`,
  build the direction toward the reference point, return
  `intensity * falloff / dist²`.
* `pdfDirect` → `1.0` if `measure == EDiscrete`, else `0`.

Notice that `sampleDirect` **multiplies by 1/dist²** itself (path tracer
expects this). Compare with `area.cpp` where the light is a surface — there
the PDF is converted from area-measure to solid-angle measure with the same
`dist²/cos` factor.

## 3.4 Walk-through: `directional.cpp`

For a **directional** light (the obvious choice for showcasing layered
materials at a chosen grazing angle):

```cpp
m_type = EDeltaDirection;
```

* The "position" of a directional light is at infinity. Mitsuba represents
  this by sampling a position on a sphere bounding the scene
  (`m_sceneAABB`), so `samplePosition` actually returns a point — that's why
  `Scene::configure` calls into your emitter with the scene bounds.
* `sampleDirect` is delta in direction: just return the light's fixed
  direction, `dist = ∞`, `measure = EDiscrete`, `pdf = 1`.

For the secondary task you have several easy options:

1. Take `directional.cpp` and add a falloff curve in *position*, not
   direction (so the beam is finite-width along the perpendicular plane).
   This is a "laser stripe" that makes Fresnel/grazing effects in your
   layered material trivial to read off.
2. Or extend `spot.cpp` with a tighter angular profile (cosine-power) and a
   parametric IES-style falloff.
3. Or a "collimated rectangular beam" — point + delta direction + a
   spatial mask. This is useful for measuring transmittance through your
   stack at a chosen incidence angle.

## 3.5 How the path tracer uses emitters (so you know the contract)

From `path.cpp`:

```cpp
DirectSamplingRecord dRec(its);
Spectrum L = scene->sampleEmitterDirect(dRec, sampler->next2D());
// L is already L_e * cos / pdf_direct, multiplied by visibility
// (occlusion check is inside sampleEmitterDirect)
if (!L.isZero()) {
    BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);
    Spectrum f = bsdf->eval(bRec);
    Float bsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
                  ? bsdf->pdf(bRec) : 0.0f;
    Float w = miWeight(dRec.pdf, bsdfPdf);   // power heuristic
    Li += throughput * L * f * w;
}
```

Two implications for your custom emitter:

* `sampleDirect` returns `L_e * (the right Jacobian)` already divided by
  `dRec.pdf` is **not** what happens here — actually `sampleEmitterDirect`
  divides by the joint PDF inside `Scene`, so you should follow the
  convention in `spot.cpp` where you return `intensity * falloff / dist²`
  and set `dRec.pdf` to the discrete probability of choosing that emitter
  position.
* If your light is delta-direction or delta-position, `bsdfPdf` is forced
  to 0 in the path tracer above (because `dRec.measure != ESolidAngle`),
  which is correct: BSDF sampling can never hit a delta light, so MIS
  collapses to pure NEE — exactly what you want.

This is also why your **LayeredBSDF must support NEE**: it needs `eval()`
and `pdf()` on arbitrary direction pairs, which means the random-walk
estimator can't just produce samples — it must also be queryable.
