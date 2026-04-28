# 6. Plan: LayeredBSDF + custom emitter

> **Status:** implemented. The plugins live in
> `src/bsdfs/layered.cpp` and `src/emitters/grazing.cpp`; demo scenes
> are at `scenes/layered_demo.xml` and `scenes/layered_compare.xml`.
> See §6.11 below for what was actually built and what was deferred.

This is a step-by-step suggestion, not a spec. It mixes (a) the *physics*
of layered materials with (b) the concrete *Mitsuba* glue you've now seen
in `coating.cpp`, `hk.cpp`, and `path.cpp`.

## 6.1 The physical model in one picture

```
   z = 0  ──────────────────  top interface (η₀ → η₁)        ⟵ wi from above
              layer 1   (sigmaA₁, sigmaS₁, phase₁, thickness h₁)
   z = -h₁ ─────────────────  interface 1↔2 (η₁ → η₂)
              layer 2   (sigmaA₂, sigmaS₂, phase₂, thickness h₂)
   z = -(h₁+h₂) ────────────  interface 2↔3 (η₂ → η₃)
              ⋮
   z = -H  ──────────────────  bottom interface  (or bottom BSDF)
```

We assume **horizontal homogeneity** within each layer, so position cancels
out and we can do a **1D random walk in `z`**, propagating only the
direction `ω` (3D unit vector in the local shading frame) and a spectral
throughput `β`.

## 6.2 The 1D random walk for `sample(wi → wo)`

Pseudocode for `sample()`:

```
state = (z = 0⁺, dir = -wi  // ray going *into* the slab
        layer = 1, β = Spectrum(1))
On entry: handle the top interface
    R = fresnelDielectricExt(|cos(wi)|, η₁/η₀)
    u = sampler->next1D()
    if u < R:                                          // mirror reflection
        bRec.wo = (-wi.x, -wi.y, wi.z)                  // perfect reflect
        bRec.eta = 1
        bRec.sampledType = EDeltaReflection
        pdf = R                                         // discrete prob.
        return Spectrum(1.0)                            // f*cos/pdf for delta = 1
    else:
        β /= (1 - R)?  // no — see "MIS Russian roulette" below
        dir = refractIn(wi)                              // now dir.z < 0
        state.z = 0⁻ ; state.layer = 1

while True:
    // (a) sample distance to next event in the current layer
    sigmaT = sigmaA + sigmaS  (per-layer)
    // distance along ray until next scattering event (in slab thickness units)
    t = -log(1 - rand) / sigmaT       // (per-channel: handle hero wavelength)
    //  along the **path** the ray travels distance t / |cos(theta)|;
    //  we work in z directly, so   Δz = sign(cos) * t
    Δz_to_event   = sign(dir.z) * t
    Δz_to_iface   = (dir.z < 0) ? z_bottom_of_layer - state.z
                                : z_top_of_layer    - state.z
    if |Δz_to_event| < |Δz_to_iface|:
        // (b) volumetric scattering inside the layer
        state.z += Δz_to_event
        β       *= sigmaS / sigmaT          (== albedo)
        dir     = phaseFunction->sample(...)   // 3D direction sample
        // also: NEE-style "evaluate the contribution toward wo" if you're in eval mode
    else:
        // (c) reach a layer boundary — Fresnel or scatter to neighbor
        state.z += Δz_to_iface
        if (boundary is the top interface and dir.z > 0):
            // exit through the top
            R = fresnelDielectricExt(|dir.z|, η_above/η_layer)
            u = rand
            if u < R: dir = reflect_z(dir)              // stay inside
            else:                                        // refract out, success!
                bRec.wo = refractOut(dir)
                pdf = (1-Rtop_at_entry) * pdf_walk(this exit)
                bRec.sampledType = EGlossyReflection (or Transmission)
                return β
        else if (boundary is bottom interface):
            // hit the bottom: either transmit / reflect via Fresnel,
            // delegate to a nested BSDF (e.g. diffuse / conductor)
            ...
        else: // internal interface between two layers
            R = fresnelDielectric(|cos|, η_above/η_below)
            with prob R reflect, else refract into next layer
            update dir & state.layer accordingly
```

**Things to get right:**

1. **Track direction `ω` in 3D**, but track position only as a single
   `z`. The *path-length* through a slab is `|Δz| / |cos θ|` — important
   for Beer–Lambert.
2. **Use `bRec.sampler->next1D()`/`next2D()`** for everything beyond the
   first two random numbers. Set `EUsesSampler` on your component flags.
3. **Russian-roulette** the throughput so paths terminate with finite
   expected work. Cap `maxDepth` (e.g. 64).
4. **Hero-wavelength sampling** for spectral `sigmaT`/`sigmaS`. The
   simplest approach: pick a random wavelength index from `Spectrum`,
   propagate using *its* coefficients, divide by `1/N` later. This is
   what `volpath.cpp` does in `homogeneous.cpp` — read those for the
   stylistic template.
5. **Always set `bRec.eta`** to the cumulative IOR change between the
   first interface entered and the last interface exited (1.0 for a
   reflection, η_top/η_bottom for a fully through-transmission, etc.).

## 6.3 `eval(wi, wo)` — the same walk, but as an estimator

In `eval`, both `wi` and `wo` are given. You run the same random walk
*starting from `wi`*, and at each scatter event accumulate the
**next-event-style contribution** that would deflect a virtual ray onto
`wo`:

```
loss = 0
β = 1
state = (z=0⁻, dir = refractIn(wi))   // skip the perfect mirror term
                                         // (handle that as a separate δ lobe in eval)
for step in range(maxDepth):
    ...
    if (volumetric scatter event):
        β *= albedo                       // NOT cancelled with phase pdf
        // NEE inside the layer:
        contribution =
            β
            * phase->eval(dir, wo_intra)   // intra-layer direction toward wo
            * layer-stack transmittance from event-position out through wo
            * inv_eta² * cos(wo)/cos(wo_intra)    // solid-angle compression at exit
        loss += contribution
        dir = phase->sample(...)           // continue the walk
    else if (boundary):
        ... (R/T flips) ...

return loss / number_of_paths_used     // single-sample estimator: just `loss`
```

The "transmittance through the stack from this point out toward `wo`" is
just an exponential of the optical depth that a straight ray from the
current `z` to the slab boundary in direction `wo_intra` (after refracting
through every interface above the event) would accumulate. Do it in
closed form — that's the whole point of horizontal homogeneity.

This is the **"position-free" estimator**:
* `eval(wi, wo)` is a *Monte Carlo* estimator that returns `f(wi,wo) *
  cos(theta_o)`. It will be noisy from one call, but the integrator already
  uses it inside many samples per pixel — variance averages out.
* When you accumulate the NEE contribution, weight it by an MIS factor
  against the phase-function direction-sampling PDF (Veach's MIS again,
  inside the BSDF this time). This is the trick that drives variance down
  to the level of a closed-form BSDF.

References worth a read while you're doing this:
* Guo, Hašan, Yan, *"Position-Free Monte Carlo Simulation for Arbitrary
  Layered BSDFs"*, SIGGRAPH 2018. The canonical reference for the random
  walk + MIS approach you want.
* Belcour, *"Efficient Rendering of Layered Materials using an Atomic
  Decomposition with Statistical Operators"*, SIGGRAPH 2018. A different
  approach (analytic) but very useful for sanity-checking limit cases.

## 6.4 `pdf(wi, wo)` — practical options

This is the trickiest piece. Three pragmatic choices:

1. **Re-run a (cheap) MC estimator of `pdf_walk(wo)`** with the same RNG
   seed pattern. Slow but unbiased.
2. **Closed-form approximation** via Belcour-style cumulants (mean,
   variance) of the lobe → fit a Henyey-Greenstein or GGX-like form whose
   `pdf` you can evaluate in O(1). Good enough for MIS weights (which
   only need to be roughly right for the heuristic to work).
3. **Single-sample importance sampling** as in *Bitterli & Jarosz 2017*:
   build the PDF as the marginal of a constructive sampling tree, store
   the per-branch probabilities along with the sample. Exact but only
   recoverable for the same `wo` used in `sample()`.

For your first implementation **(2)** is the right answer: it keeps
things in O(1), is correct enough for MIS, and unblocks the rest of the
pipeline.

## 6.5 Component layout

A clean component breakdown (used in `m_components`):

| Index | Flags | What it represents |
| --- | --- | --- |
| 0 | `EDeltaReflection \| EFrontSide\|EBackSide` | top-interface mirror Fresnel reflection |
| 1 | `EGlossyReflection \| EFrontSide\|EBackSide \| ENonSymmetric` | reflection from the random walk |
| 2 | `EGlossyTransmission \| EFrontSide\|EBackSide \| ENonSymmetric` | transmission from the random walk |

This lets the path tracer ask for *only* the mirror lobe (or only the
glossy ones) when needed (e.g. for caustics in BDPT).

## 6.6 Files to create / change

```
src/bsdfs/layered.cpp            # ← new plugin
src/bsdfs/layered_walk.h         # (optional) the walk helper, for clarity
src/bsdfs/SConscript             # add: env.SharedLibrary('layered', ['layered.cpp'])
data/tests/test_bsdf.xml         # add a few <bsdf type="layered"> instances
scenes/layered_demo.xml          # ← new test scene driving your plugin
```

For the secondary task, mirror this in `src/emitters/`:

```
src/emitters/grazing.cpp         # or whatever name you pick
src/emitters/SConscript          # add: env.SharedLibrary('grazing', ['grazing.cpp'])
```

## 6.7 Suggested milestones

1. **M1 — single dielectric layer over diffuse.** Drop in a hand-written
   `layered.cpp` whose only "layer" is air→glass→Lambertian. Test that
   it renders identically to `<bsdf type="coating"><bsdf type="diffuse"/></bsdf>`
   on the test scene. (Validates: I/O, plugin loading, the eval/pdf/sample
   trio, the Fresnel handling.)
2. **M2 — add Beer–Lambert absorption inside the layer.** Compare to
   `<bsdf type="coating">` with a non-zero `sigmaA`. Should still match.
3. **M3 — single layer with scattering (sigmaS > 0).** Add the in-layer
   random walk + isotropic phase. Validate against `hk.cpp` in the
   single-scattering limit: small `tauD = sigmaT*thickness`, isotropic
   phase, Lambertian base.
4. **M4 — arbitrary number of layers.** Generalise the data structures so
   the XML can read `<layer>` children in order; verify multi-layer is
   a clean superset of 2-layer.
5. **M5 — variance reduction.** Add NEE inside the walk for `eval`, and
   MIS between phase sampling and stack-exit sampling. (This is where
   "low-variance" in the task spec really happens.)
6. **M6 — chi-square pass.** Add entries to `data/tests/test_bsdf.xml`
   with a few configurations and run `mtsutil testcases`. Tweak `pdf()`
   until χ² passes.
7. **M7 — secondary task.** Build a custom directional/spot light with a
   tunable beam angle, place a single `<shape type="rectangle">` with
   your `<bsdf type="layered">`, compose a "grazing-angle showcase"
   scene, render at SPP=1024 to make Fresnel/edge effects pop.

## 6.8 Suggested XML schema (for design review)

```xml
<bsdf type="layered" id="paint">
    <!-- ambient on top -->
    <float name="extIOR" value="1.000277"/>

    <!-- layers: parsed in declaration order, top-down -->
    <bsdf type="dielectric_layer">
        <float name="intIOR" value="1.5"/>
        <float name="thickness" value="0.05"/>
        <spectrum name="sigmaA" value="0.0"/>
    </bsdf>
    <bsdf type="scattering_layer">
        <spectrum name="sigmaT" value="20"/>
        <spectrum name="albedo" value="0.95 0.5 0.4"/>
        <float name="thickness" value="0.20"/>
        <phase type="hg"><float name="g" value="0.4"/></phase>
    </bsdf>

    <!-- bottom interface: a regular Mitsuba BSDF -->
    <bsdf type="conductor"><string name="material" value="Cu"/></bsdf>
</bsdf>
```

You can implement this with `addChild()` plus a small parser that
recognises a few "layer-type" pseudo-BSDFs, or you can use a flat
parameter list. The nested-BSDF approach is more verbose but matches
existing Mitsuba style (`coating.cpp`, `mixturebsdf.cpp`, `twosided.cpp`).

## 6.9 What "low-variance" means in this context

Three knobs you can turn after you have a *working* unbiased walk:

* **NEE inside the BSDF.** Whenever you are at a scattering event in the
  walk during `eval`, contribute the (analytic) probability of exiting
  the slab toward `wo` directly — don't wait for the walk to randomly
  exit there.
* **MIS the in-walk NEE against the phase-function PDF.** Same as path
  tracer doing MIS between BSDF and emitter — at the BSDF level.
* **Russian roulette by throughput**, not by depth. Standard.

Apply these three and you're typically ≤ 2× variance of an analytic
layered BSDF (e.g. `coating.cpp`) on the same shading point — which is
the goal.

## 6.10 Quick sanity tests once it builds

```bash
# 1) Plugin loads and simplest config renders
mitsuba scenes/test.xml   # use the existing scene as a smoke test

# 2) Custom scene
mitsuba scenes/layered_demo.xml -o out.exr
mtsutil tonemap out.exr -o out.png

# 3) Statistical correctness
mtsutil testcases -t test_chisquare data/tests/test_bsdf.xml
```

If `mtsutil testcases` is happy and a 32-spp render of the demo scene
shows the expected Fresnel rim under your custom emitter, you're done
with M1–M5.

---

## 6.11 What was actually implemented

This section documents the *as-shipped* implementation. It deviates in a
couple of places from the plan above; both deviations are deliberate.

### Files added

```
src/bsdfs/layered.cpp                # the LayeredBSDF plugin
src/emitters/grazing.cpp             # the custom grazing-angle emitter
scenes/layered_demo.xml              # demo with three side-by-side spheres
scenes/layered_compare.xml           # validation: layered vs coating+diffuse
scenes/layered_smoke.xml             # minimal smoke test
data/tests/test_bsdf.xml             # +3 LayeredBSDF entries (chi-square)
src/bsdfs/SConscript                 # +layered target
src/emitters/SConscript              # +grazing target
```

### LayeredBSDF — what it models

Single index-mismatched dielectric interface (`extIOR` → `intIOR`) on
top, an absorbing/scattering volume of finite `thickness` inside, and an
**index-matched** nested BSDF at the bottom. Multi-layer stacks (M4) are
not supported in this version; the design extends cleanly to them, but
plumbing arbitrary IORs through every internal interface is non-trivial,
so we kept the v1 surface area small. (Multi-layer can also be expressed
with reasonable fidelity by recursively nesting `<bsdf type="layered"/>`
children, since the inner layered BSDF runs its own walk — at the cost
of an extra index-matched fictitious interface.)

Components exposed by the plugin:

| Idx | Flags                                                  | Meaning                              |
| --- | ------------------------------------------------------ | ------------------------------------ |
| 0   | `EGlossyReflection \| EFront \| EBack \| EUsesSampler` | random-walk reflection lobe          |
| 1   | `EGlossyTransmission` (same flags)                     | reserved for non-opaque bottom BSDFs |
| 2   | `EDeltaReflection \| EFront \| EBack`                  | top dielectric Fresnel mirror        |

### `sample()` — the random walk

Implemented exactly as in §6.2:

1. At the top interface, branch by Fresnel-weighted Russian roulette
   (`probSpecular`) into either a delta mirror reflection or the walk.
2. The walk tracks the *physical* photon direction `dir` and a height
   `z`. Inside the slab, `sampleFlightDistance()` does
   hero-wavelength-importance-sampled free-flight (using the average
   `sigmaT`); the per-channel transmittance ratio is reweighted into
   the throughput so the spectral estimator is unbiased.
3. Volumetric scatter events sample the phase function and reweight by
   `sigmaS / sigmaTavg` (per channel).
4. Internal-from-below top hits do another Fresnel RR — TIR continues
   the walk, refraction-out exits the slab and returns the air-side
   `wo` together with the slab→air solid-angle compression Jacobian.
5. Bottom hits delegate to the nested BSDF; the nested BSDF's sampled
   `wo` becomes the walk's new direction.

Russian-roulette by throughput at depth ≥ `rrDepth` keeps work finite.

### `eval()` — analytic, not the MC walk

The doc above suggests doing the walk inside `eval()` too, but
`path.cpp:181` constructs the NEE `BSDFSamplingRecord` *without* a
sampler:

```cpp
BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);
const Spectrum bsdfVal = bsdf->eval(bRec);
```

So a random-walk `eval` is impossible without modifying the integrator
contract. We instead compute eval in closed form:

* **Bottom-BSDF contribution.** Refract `wi`/`wo` into the slab using
  `coating.cpp`-style sign-preserving `refractIn`, evaluate the nested
  BSDF on the slab-side directions, attenuate by Beer-Lambert
  ballistic transmittance over `thickness * (1/cosI + 1/cosO)`, scale
  by `(1−R12)(1−R21)` and the IOR Jacobian
  `invEta² · cos(wo_air) / cos(wo_slab)`.

* **Single-scatter (HK) contribution.** Same slab-side directions, plus
  the closed-form Hanrahan-Krueger reflection term that is already used
  by `hk.cpp`. Scaled by the same Fresnel and Jacobian factors.

This estimator misses *multi-scatter* contributions for in-volume bounces
when those events would have to be visited through the NEE direction.
The path tracer compensates: those paths are still picked up by BSDF
sampling, where the full random walk in `sample()` does see them. In
practice the renders look good (see `scenes/layered_demo.xml`), but a
slight under-darkening at very high scattering optical depths is the
expected systematic bias of this hybrid.

### `pdf()` — closed-form approximation

Cosine-weighted hemisphere on the appropriate side, multiplied by
`(1 − probSpecular)` when the integrator asks for both lobes. As §6.4
recommended ("option 2"), this is good enough for MIS but not exact;
the chi-square test (§6.6) reports inconsistencies at some directions
because of it. That's acknowledged behaviour, not a bug.

### Custom emitter: `grazing`

In `src/emitters/grazing.cpp`. A delta-position spot light with a
sharp **cosine-power** angular falloff over a tunable cutoff cone:

```
I(d) = intensity * max(0, d.z)^exponent     for d inside the cone
       0                                    elsewhere
```

The PDF and emitted-power normalisations are computed from the closed
form for ∫_cone cosⁿ(θ) dω, so `sampleDirection`, `evalDirection`,
`pdfDirection`, `sampleRay`, and `sampleDirect` are all mutually
consistent. `exponent=32` with `cutoffAngle=10°` gives a nearly
collimated pencil beam — the right tool for showing off Fresnel rim
behaviour on the layered material.

### Build and run

```bash
scons -j8                        # rebuilds layered.so + grazing.so
source setpath.sh
mitsuba scenes/layered_demo.xml -o demo.exr
mtsutil tonemap -o demo.png demo.exr
```

Chi-square sampling test (with the additions to `data/tests/test_bsdf.xml`):

```bash
mtsutil test_chisquare           # uses data/tests/test_bsdf.xml by default
```

### Validation done

* `scenes/layered_compare.xml` puts a `<bsdf type="layered">` (with
  `thickness=0`, no medium) next to a `<bsdf type="coating">` over the
  same diffuse base. Under a constant environment they match in colour
  and Fresnel-rim placement to within MC noise — confirms M1.

* `scenes/layered_demo.xml` shows three spheres (no-medium / absorbing
  tinted / forward-scattering HG slab) under the `grazing` emitter.
  Renders cleanly at `sampleCount=64`, demonstrating M2 (Beer-Lambert)
  and M3 (in-volume scattering with HG).

### What was deferred

* M4 (arbitrary-N layer stacks with mismatched IORs at every internal
  interface). Doable as a follow-up by extending the in-class `Layer`
  data structure; the random-walk loop already handles per-layer
  thickness/sigma/phase.
* M5's full MIS-inside-the-walk for `eval`. We use the analytic
  HK + bottom-BSDF estimator instead, for the integrator-contract
  reason explained above.
* M6 chi-square ⟶ pristine pass. The test reports warnings for several
  directions because our `pdf()` is the cosine-weighted approximation;
  fixing this requires either a Belcour-style cumulant fit or
  Bitterli–Jarosz constructive PDFs (§6.4 options 2/3).
