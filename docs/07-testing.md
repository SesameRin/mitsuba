# Testing the LayeredBSDF + Grazing Emitter

This page collects every command you need to validate the two plugins
shipped in `src/bsdfs/layered.cpp` and `src/emitters/grazing.cpp`,
together with a description of how each test scene is structured and
what it is actually checking.

The tests fall into three layers:

1. **Visual rendering tests** — three hand-crafted scenes that exercise
   different parts of the implementation and let you eyeball the result.
2. **Reference comparison** — render LayeredBSDF in its degenerate
   regime side-by-side with the existing `coating` plugin; they must
   match.
3. **Statistical (chi-square) test** — `mtsutil testcase` draws ~100k
   samples per incident direction and compares the empirical histogram
   against the analytic `pdf()`.

## 0. Setup

You need to do this once per shell:

```bash
cd /home/ula/pbr/mitsuba
source setpath.sh        # puts dist/ on PATH, sets MITSUBA_PLUGIN_PATH
scons -j8                # incremental build; should report "up to date"
```

Confirm the two plugin shared libraries exist:

```bash
ls -la dist/plugins/layered.so dist/plugins/grazing.so
```

If either is missing, the build silently failed for that plugin —
re-run `scons -j8` and read its stderr. If `mitsuba` reports a
"plugin not found" error at render time, you forgot to source
`setpath.sh` (it sets `MITSUBA_PLUGIN_PATH=$PWD/dist/plugins`).

## 1. Visual tests

All three scenes live in `scenes/`. Each is small (256–640 px wide,
32–64 spp) and renders in seconds on a modern CPU.

### 1a. `scenes/layered_compare.xml` — reference comparison

**What it tests.** When `thickness=0`, `sigmaA=0`, `sigmaS=0`, the
LayeredBSDF reduces to "smooth dielectric layer over diffuse," which
is exactly what `coating + diffuse` produces. The scene renders a
sphere of each side-by-side under a constant environment.

**Why this is the most important test.** It validates the entire
top-interface logic — Fresnel reflection/transmission probabilities,
Snell refraction sign convention, the solid-angle compression Jacobian
for refracted exit rays — against a plugin (`coating`) that already
ships in Mitsuba and is known to be correct.

**How to run:**

```bash
mitsuba scenes/layered_compare.xml -o tmp/compare.exr
mtsutil tonemap -o tmp/compare.png tmp/compare.exr
```

**Expected result.** Both spheres render essentially identically
(within Monte Carlo noise). If they don't:

- A bias in the Fresnel branch will show as one sphere being uniformly
  brighter / darker than the other.
- A wrong refraction sign will make the LayeredBSDF sphere look
  noticeably less saturated than `coating` (because the diffuse base is
  receiving the wrong amount of incoming light).
- A wrong solid-angle Jacobian usually shows up as a narrow rim
  intensity mismatch at grazing angles.

### 1b. `scenes/layered_smoke.xml` — single-material smoke test

**What it tests.** A single layered sphere
(`thickness=0.5`, `sigmaA=(0.2, 0.5, 0.8)`, no scattering) lit by the
new `grazing` emitter aimed nearly horizontally. This exercises:

- the absorbing-medium path in the random walk,
- Beer–Lambert transmittance along the in-layer ray,
- the grazing emitter's cosine-power angular falloff and direct-sample
  contribution to NEE.

**How to run:**

```bash
mitsuba scenes/layered_smoke.xml -o tmp/smoke.exr
mtsutil tonemap -o tmp/smoke.png tmp/smoke.exr
```

**Expected result.**

- The base diffuse colour is `(0.8, 0.2, 0.2)` (red).
- `sigmaA = (0.2, 0.5, 0.8)` absorbs green and blue more than red, so
  the layer transmits red preferentially — the lit side should look
  warmer than the base reflectance alone would suggest, with a faint
  cyan-shifted shadow boundary.
- A bright Fresnel rim appears on the side facing the grazing light.

### 1c. `scenes/layered_demo.xml` — three-sphere showcase

**What it tests.** All three regimes simultaneously, side-by-side:

| Sphere | Parameters | What it demonstrates |
|---|---|---|
| Left  | `thickness=0`, no sigma                                                  | Dielectric Fresnel on red diffuse — should match the `coating` analogue. |
| Mid   | `thickness=0.6`, `sigmaA=(0.5, 1.0, 1.5)`                                | Pure absorption — grazing rays travel further through the layer and tint more strongly. |
| Right | `thickness=0.4`, `sigmaA=(0.05,0.1,0.2)`, `sigmaS=3.0`, HG `g=0.5`       | Volumetric scattering inside the layer — whitens the grazing rim. |

The committed `scenes/layered_demo.png` is the reference output.
Re-render and compare:

```bash
mitsuba scenes/layered_demo.xml -o tmp/demo.exr
mtsutil tonemap -o tmp/demo.png tmp/demo.exr
```

The third sphere should clearly look "milkier" near the rim than the
middle one — that is the forward-scattering HG phase function brightening
the silhouette. If the rim looks identical to the absorbing sphere,
either `sigmaS` is being ignored or the phase-function sample is not
being applied during the walk.

## 2. Statistical test (sample/PDF consistency)

The chi-square test is configured via `data/tests/test_bsdf.xml`.
That file lists every BSDF that should be tested for sampling/PDF
consistency; lines 219–247 contain three `<bsdf type="layered">`
entries that mirror the three demo scene regimes.

**What it does.** For each BSDF, `mtsutil testcase`:

1. picks a sweep of incident directions,
2. calls `sample()` ~100k times per direction and bins the outgoing
   directions over the sphere,
3. computes the expected per-bin frequencies from `pdf()` and the
   weight returned by `sample()`,
4. runs Pearson's chi-square test and reports a per-direction
   p-value.

**How to run:**

```bash
cd tmp                                            # so failure_*.m dumps stay out of the repo
mtsutil testcase ~/pbr/mitsuba/data/tests/test_bsdf.xml
```

**Expected result.**

- The first two LayeredBSDF entries (no-scattering and absorbing) pass
  cleanly.
- The scattering entry produces "Potential inconsistency" warnings.
  This is **known and documented** — see
  [`06-layered-bsdf-task.md` §6.11](06-layered-bsdf-task.md#611-what-was-actually-implemented).
  The shipped `pdf()` is a cosine-weighted approximation; matching it
  exactly to the random walk's true distribution is M6, deferred.

Whenever a histogram bin disagrees with `pdf()` past the tolerance,
`mtsutil testcase` writes an Octave/MATLAB file `failure_<N>.m` into the
current working directory containing the empirical and reference
tables. These are diagnostic dumps only — delete them when you're done.

## 3. Quick failure-mode reference

| Symptom | Likely cause |
|---|---|
| Sphere renders pure black | `layered.so` not loaded — re-source `setpath.sh`, check `dist/plugins/`. |
| `compare.png` halves don't match | Top-interface regression: Fresnel branching, refraction sign in `refractIn`, or the solid-angle Jacobian on exit. |
| `demo.png` middle sphere looks gray instead of tinted | Absorption / Beer–Lambert path is broken — the walk is exiting without accumulating `exp(-sigmaT * t)`. |
| `demo.png` scattering sphere identical to absorbing one | `sigmaS` is being ignored at sampling time, or the phase function is never invoked. |
| Grazing emitter renders flat / unfocused | `cutoffAngle` was passed in radians — the XML expects degrees. |
| Renderer aborts with `bRec.sampler != NULL` | Something is calling the random-walk path from inside `eval()`. The shipped `eval()` is analytic precisely to avoid this; do not regress it. |
| Many `failure_*.m` files in cwd | Chi-square reported inconsistencies — expected for the scattering layer (M6). Delete the files. |

## 4. Total runtime

On a recent laptop the full suite (build + three renders + chi-square)
takes well under a minute:

| Step | Approx. wall time |
|---|---|
| `scons -j8` (incremental) | < 5 s |
| `layered_compare.xml`     | ~5 s |
| `layered_smoke.xml`       | ~5 s |
| `layered_demo.xml`        | ~15 s |
| `mtsutil testcase`        | ~30 s |
