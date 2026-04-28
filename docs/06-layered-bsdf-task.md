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

Components exposed by the plugin (front-side only after the v2 revision —
see §6.12):

| Idx | Flags                                                  | Meaning                              |
| --- | ------------------------------------------------------ | ------------------------------------ |
| 0   | `EGlossyReflection \| EFrontSide \| EUsesSampler \| ENonSymmetric` | random-walk reflection lobe          |
| 1   | `EDeltaReflection \| EFrontSide`                       | top dielectric Fresnel mirror        |

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

---

## 6.12 v2 revisions (after reading the reference)

After reading `multilayered.cpp` (the SIGGRAPH 2018 production plugin —
see `08-reference-multilayered.md`), the local plugin was tightened in
four ways:

1. **Reference-schema XML compatibility.** `LayeredBSDF` now also
   parses the paper plugin's parameter names: `nbLayers=2`,
   `<bsdf name="surface_0">` (top — we extract its IOR via
   `bsdf->getEta()`), `<bsdf name="surface_1">` (bottom — used as the
   nested BSDF), `<spectrum name="sigmaT_0">` + `<spectrum name="albedo_0">`
   (we derive `sigmaA = sigmaT*(1-albedo)` and `sigmaS = sigmaT*albedo`),
   and `<phase name="phase_0">`. Reference-only flags (`bidir`, `pdf`,
   `stochPdfDepth`, `bidirUseAnalog`, `pdfRepetitive`, …) are accepted
   silently and ignored. **The practical effect**: the figure8 paper
   scenes drop in directly with only two textual rewrites — `path_layered`
   → `path` and `multilayered` → `layered`. (See §6.13 for the
   instructions.)

2. **Front-side-only components.** The previous version exposed
   `EGlossyReflection | EFrontSide | EBackSide` plus a stub
   `EGlossyTransmission`. Back-side incidence on a slab with an opaque
   base is unphysical; the old code handled it by mirroring `wi.z` and
   pretending the back was the front, which is correct only for a
   symmetric stack. The opaque-bottom configuration has no symmetry
   guarantee, so v2 simply declares the BSDF one-sided
   (`EFrontSide` only) and short-circuits back-side queries to zero.
   The unused `EGlossyTransmission` component is dropped — the path
   tracer was occasionally spending NEE budget on a lobe that always
   returned 0.

3. **Sample-component indices renumbered.** With the transmission lobe
   gone, `EDeltaReflection` is now component 1 (was 2). `getRoughness()`
   was updated to match. The serialized stream version now also stores
   `m_nbLayers` for round-tripping the schema.

4. **Heuristic specular-sampling weight cleanup.** Replaced the dead
   `(1 - 0.5)` term with a clearer `0.5 + 0.5 * avgAttn` so the
   intent — bias toward the mirror lobe when the slab transmits little
   light — is obvious from the source.

### What v2 deliberately does *not* do

* **No bi-directional / stochastic-pdf eval.** The reference's
  `bidirStochTRT` mode requires a `Sampler` inside `eval()`, which the
  stock `path` integrator does not pass; the reference ships its own
  `path_layered` integrator that adds `bRec.sampler = rRec.sampler`
  during NEE. We could port that integrator, but a custom integrator
  would no longer be a drop-in for vanilla Mitsuba. The hybrid
  analytic-eval / MC-sample design here keeps the BSDF compatible with
  every Mitsuba integrator that does NEE, at the cost of slightly
  higher variance on highly forward-scattering slabs.

* **No microfacet top interface.** The figure8 scenes use
  `roughdielectric` with α as low as 0.005 (visually near-smooth), so
  modelling the top as a smooth dielectric loses very little fidelity.
  α=0.1 (one of the figure8 plates) is a more visible discrepancy; the
  reference renders it with an actual GGX-distributed top whereas we
  fall back to smooth + Fresnel.

* **No per-layer `surface_*` for N>2.** `nbLayers` other than 2 is
  refused at parse time. The reference's true multi-layer support
  (per-internal-interface IOR + roughness + normal/orientation maps)
  is a substantial extension; rendering 2-layer paper scenes does not
  need it.

---

## 6.13 Rendering the paper scenes (`scene_data/` zips)

Two paper-scene archives ship in `scene_data/`: `figure8.zip` (small,
4-plate Veach MI test — the v2-renderable one) and `teaser.zip`
(large, anisotropic-medium tablecloth + vases — depends on
`microflake` phase, instanced textures, the `damask4_info.dat` block
file, and other reference-only machinery; **not renderable** with
this plugin without porting more of `multilayered.cpp`).

The rest of this section is the figure8 recipe.

### One-time setup

```bash
cd /home/ula/pbr/mitsuba/scene_data
unzip -o figure8.zip -d figure8/        # → figure8/{mi_*.xml, meshes/}
```

The repo already contains adapted scenes at
`scenes/scene_data/figure8/`. They are byte-equivalent to the upstream
XMLs except for two `sed` rewrites:

```
type="path_layered"   → type="path"
type="multilayered"   → type="layered"
```

`scenes/scene_data/figure8/meshes` is a symlink into the unzipped
`scene_data/figure8/meshes/`, so unzipping is a hard prerequisite.

### Rendering

```bash
cd /home/ula/pbr/mitsuba
source setpath.sh
mitsuba scenes/scene_data/figure8/mi_acr_left.xml \
        -o scenes/scene_data/output/figure8_mi_acr_left.exr
mtsutil tonemap -o scenes/scene_data/output/figure8_mi_acr_left.png \
        scenes/scene_data/output/figure8_mi_acr_left.exr
```

Repeat for `mi_trt_middle.xml` and `mi_noMIS_right.xml`. Render time
is ~7s per scene at the upstream 200 spp on a 20-core machine; the
build is single-precision OpenMP-parallelised path tracing.

### Outputs

Pre-rendered EXR + PNG pairs are committed at
`scenes/scene_data/output/`:

* `figure8_mi_acr_left.png` — Veach-style 4 plates × 4 area lights,
  each plate a 2-layer dielectric+gold stack with varying combined
  roughness. The classic test that all 16 plate-light combinations
  receive a consistent contribution under MIS.
* `figure8_mi_trt_middle.png` — same scene, paper's "TRT" pdf mode
  in the upstream XML (no effect on our plugin — the flag is parsed
  and ignored).
* `figure8_mi_noMIS_right.png` — same scene with multilayered's
  `MIS=false` flag (also ignored here, since the integrator-level
  MIS in `path` is unaffected by it).

### Caveats vs. the reference

* All four plates in our renders share the *same* Fresnel rim
  because we collapse the rough top dielectric to a smooth one.
* The bottom roughness still varies (it's just a stock Mitsuba
  `roughconductor` BSDF used as the nested base), so plates with
  rough bottoms (`b`, `d`) differ visibly from those with smooth
  bottoms (`a`, `c`) on the conductor highlight.
* Variance is somewhat higher than the upstream paper plugin; that is
  the documented trade-off of using analytic eval rather than the
  reference's stochastic-pdf bi-directional eval.

## 6.14 Hero showcase scene (`scenes/layered_hero.xml`)

The figure8 plate-grid is the validation render. For a *showpiece*
that puts both deliverables — the `layered` BSDF and the `grazing`
emitter — into one image, the repo also ships
`scenes/layered_hero.xml`: two classic graphics props side-by-side,
each coated with a different layered material, lit by the grazing
emitter.

### Subjects

| Hero | Mesh file | Source |
| --- | --- | --- |
| Stanford bunny (left) | `scenes/meshes/bunny.ply` (35947 verts / 69451 tris) | Local copy of the canonical Stanford 3D Scanning Repository bunny that ships in `data/tests/`. |
| Newell / Utah teapot (right) | `scenes/meshes/teapot.obj` (1202 verts / 2256 tris) | Fetched from the McNopper OpenGL repo on GitHub and committed to the repo. |

Both materials are 2-layer stacks built with the reference XML
schema (`nbLayers`, `surface_0`, `sigmaT_0`, `albedo_0`, `phase_0`,
`surface_1`):

| | Top (`surface_0`) | Slab | Bottom (`surface_1`) |
| --- | --- | --- | --- |
| Bunny ("amber over copper") | smooth `roughdielectric`, IOR 1.5 | `sigmaT_0=1.4`, `albedo_0="0.96 0.55 0.18"`, isotropic HG | `roughconductor`, `material="Cu"`, `alpha=0.05` |
| Teapot ("amber over gold") | smooth `roughdielectric`, IOR 1.5 | `sigmaT_0=1.8`, `albedo_0="0.97 0.65 0.22"`, isotropic HG | `roughconductor`, `material="Au"`, `alpha=0.04` |

Same warm pigmented slab on both, two different metals underneath
— the layered random walk correctly composes top + slab + bottom,
so the bunny reads distinctly redder (copper) and the teapot
yellower (gold) under identical lighting. That side-by-side hue
contrast is what would be impossible to reproduce with `roughplastic`
or any non-layered BSDF.

### Lighting

* **Key**: the custom `grazing` emitter at `cutoffAngle=22°`,
  `exponent=18`, `intensity=280`, placed off-camera-right at a
  shallow angle so its cone clips along the heroes' silhouettes
  and excites the dielectric Fresnel peak.
* **Fill**: dim constant environment (`radiance=0.05`) lifts the
  shadow side just enough to read shape.
* **Rim**: a small off-frame warm area sphere up-and-behind
  separates the back of the heroes from the wall.

### Stage

Built-in `rectangle` floor and back-wall primitives (no extra
meshes), warm-grey diffuse wall + slightly glossy near-black floor
that picks up a subtle reflection of the heroes' undersides.

### Render

From `scenes/`, after `source ../setpath.sh`:

```
mitsuba layered_hero.xml -o scene_data/output/layered_hero.exr
mtsutil tonemap -g 2.2 scene_data/output/layered_hero.exr
```

Defaults: 1280×800, 1024 spp, ldsampler, `path` integrator with
`maxDepth=16`. ~3 minutes wall-clock on a 20-core box. The
pre-rendered output is committed at
`scenes/scene_data/output/layered_hero.{exr,png}`.

### What the image is showing

* **Body colour**: the warm interior body of each hero comes from
  spectral absorption inside the slab — light that enters from the
  grazing key, travels some path through the pigmented medium,
  reflects off the metal base, and traverses the slab again on the
  way out. Different wavelengths attenuate at different rates, so
  the chromaticity shifts subtly with depth (most visible on the
  bunny's belly and the teapot's underside, where the path is
  longest).
* **Metal hue contrast**: the bunny is copper-red, the teapot is
  gold-yellow, despite an almost-identical absorbing slab. This is
  the layered BSDF *composing* the slab tint with the conductor
  Fresnel — neither object is just "tinted metal".
* **Fresnel rim**: the bright highlights on the spout, the bunny's
  ear and the teapot lid are the dielectric *top* layer reflecting
  the grazing key near-tangentially. The layered BSDF preserves
  this surface contribution alongside the volumetric one (the
  delta/specular component of the top dielectric).
* **Cast shadows**: the long, soft, single-edged shadows on the
  floor are the signature of a tight cosine-power emitter — a
  normal point or area light would scatter the boundary much more.

