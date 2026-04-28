# 8. Reference: position-free MC layered BSDF (Guo et al. 2018)

The "official" layered BSDF this project compares against is the
SIGGRAPH Asia 2018 paper

> **Position-Free Monte Carlo Simulation for Arbitrary Layered BSDFs**
> Yu Guo, Miloš Hašan, Shuang Zhao — ACM TOG (SIGGRAPH Asia 2018).

Its source lives in the sibling repo at `../layeredbsdf/` (a fork of
Mitsuba 0.6.0). This doc points at the files that actually implement the
algorithm so you can read them alongside our `src/bsdfs/layered.cpp`.

## 8.1 Core files

| File | Lines | Role |
| --- | --- | --- |
| `layeredbsdf/src/bsdfs/multilayered.cpp` | 1672 | The algorithm. Registered as the `multilayered` BSDF plugin. Implements the 1D position-free random walk (uni-dir + bi-dir), MIS, and the stochastic pdf modes. |
| `layeredbsdf/src/bsdfs/rtrans.h` | 448 | Rough-transmittance lookup at rough dielectric/conductor boundaries (tricubic interp into a precomputed 3D table). |
| `layeredbsdf/src/bsdfs/microfacet.h` | 725 | Beckmann / GGX / Phong distributions used by the layered sampler and evaluator. |
| `layeredbsdf/src/integrators/path/path_layered.cpp` | 338 | Companion integrator, registered as `path_layered`. The README recommends it over plain `path` for `multilayered`. |

`multilayered.cpp` is the entry point — start there. The constructor
parses per-layer `sigmaT_i`, `density_i`, `albedo_i`, `orientation_i`,
`normal_i`, plus the toggles `MIS`, `bidir`, `bidirUseAnalog`,
`pdfMode` (default `bidirStochTRT`), `stochPdfDepth`, `pdfRepetitive`,
`maxSurvivalProb`. Both the uni-directional and bi-directional methods
from the paper live in this single file.

## 8.2 Supporting data

| Path | Used by |
| --- | --- |
| `layeredbsdf/data/microfacet/{beckmann,ggx,phong}.dat` | Loaded by `rtrans.h` for rough-transmittance lookups. |
| `layeredbsdf/data/ior/*.spd` | Measured spectral IORs for conductors / dielectrics. |

## 8.3 Example scenes — not in the repo

The repo ships **no example `.xml` scenes**. The only `.xml` under
`layeredbsdf/` is `src/mtsgui/resources/docs.xml` (GUI help — unrelated).

The paper's figure scenes are hosted externally and linked from
`layeredbsdf/README.md`:

* teaser, figure11, figure13 — Dropbox zips
* figure2, figure3, figure8, figure12t, figure12b, figure14, figure15 —
  in the `tflsguoyu/layeredbsdf_suppl` GitHub repo under
  `github/scenes/*.zip`

Workflow once a zip is downloaded:

```
cd layeredbsdf
mv config_linux.py config.py
scons -j$(nproc)
source setpath.sh
mitsuba path/to/figureXX.xml
```

## 8.4 How a `multilayered` scene is wired

From the README — the three lines that distinguish a paper scene from a
stock Mitsuba 0.6 scene:

```xml
<scene version="0.6.0">           <!-- 0.6.0, not 0.5.0 -->
  <integrator type="path_layered"/>
  <bsdf type="multilayered"> ... per-layer sigmaT/albedo/normal/... </bsdf>
</scene>
```

Caveats called out by the upstream README:

* `conductor` and `dielectric` are unsupported as interface BSDFs —
  substitute `roughconductor` / `roughdielectric` with `alpha=0.001`.
* If `single` precision crashes or floods warnings, switch `config.py`
  to `double`.

## 8.5 Relationship to this project

Our `src/bsdfs/layered.cpp` implements a much narrower slice:
horizontally homogeneous stack, 1D random walk, uni-directional only,
designed to be readable and to exercise `path` + MIS + NEE in stock
Mitsuba 0.6. The reference `multilayered` is the production-grade
counterpart — bi-dir estimator, stochastic-pdf MIS variants,
anisotropic media, per-layer normal/orientation maps. When validating
our plugin, the reference is the ground truth to beat (or, more
honestly, to converge toward).

After reading this reference, the local plugin was tightened to accept
the reference's XML schema (`nbLayers`, `surface_0/1`, `sigmaT_0`,
`albedo_0`, `phase_0`) so paper scenes drop in with only the
integrator/BSDF type names rewritten — see
[`06-layered-bsdf-task.md §6.12`](06-layered-bsdf-task.md#612-v2-revisions-after-reading-the-reference)
for the change-list and [`§6.13`](06-layered-bsdf-task.md#613-rendering-the-paper-scenes-scene_data-zips)
for the rendering recipe. The figure8 paper scene is renderable today;
teaser is not, because it depends on `microflake` phase + instanced
textures + the `damask4_info.dat` block file from the reference repo.
