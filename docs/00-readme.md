# Mitsuba 0.6 Onboarding Notes

These docs are a personal reading guide to the Mitsuba 0.6 codebase, written
to prepare for two implementation tasks:

1. **Core task** — a new `LayeredBSDF` plugin: arbitrary stacks of dielectric /
   conductor / scattering layers, evaluated and sampled with a 1D Monte Carlo
   random walk, fully compatible with MIS and NEE in the path tracer.
2. **Secondary task** — a custom light source designed to highlight the
   layered material (e.g. a refined spot/directional light for grazing-angle
   Fresnel showcasing).

> **Status:** both tasks are implemented. See
> [`06-layered-bsdf-task.md` §6.11](06-layered-bsdf-task.md#611-what-was-actually-implemented)
> for what shipped versus what was deferred. The plugin sources are
> `src/bsdfs/layered.cpp` and `src/emitters/grazing.cpp`; demo scenes
> are in `scenes/`.

The order I suggest reading them:

| File | Purpose |
| --- | --- |
| `01-overview.md` | Bird's-eye view of the repo: directory layout, libraries, plugin model, build outputs. |
| `02-bsdf-system.md` | The `BSDF` interface in detail (`sample` / `eval` / `pdf`, conventions, flags). Walks through `diffuse.cpp` and `coating.cpp`. |
| `03-emitter-system.md` | The `Emitter` interface, `samplePosition` / `sampleDirection` / `sampleDirect`, walked through using `spot.cpp` and `directional.cpp`. |
| `04-integrators-mis-nee.md` | How the `path` integrator combines BSDF sampling and next-event estimation with multiple importance sampling. |
| `05-build-and-plugin-system.md` | SCons, `MTS_EXPORT_PLUGIN`, `setpath.sh`, where built `.so` files end up, how to add a new plugin. |
| `06-layered-bsdf-task.md` | Concrete plan for the LayeredBSDF + custom emitter tasks: theory, files to create, milestones, how to test. |
| `07-testing.md` | How to validate the shipped LayeredBSDF + grazing emitter: setup, the three test scenes, the chi-square test, and a failure-mode cheat sheet. |

> Mitsuba 0.6 is the same code-base as Mitsuba 0.5/0.4 in spirit — most of
> the academic literature on Mitsuba (e.g. Wenzel Jakob's thesis, the
> "Arbitrarily Layered Microfacet" paper, the position-free Monte Carlo paper)
> drops in directly. Anything in `src/bsdfs/` is the closest reference.
