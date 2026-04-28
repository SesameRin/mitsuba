# 1. Mitsuba 0.6 — Architectural Overview

## What Mitsuba is

Mitsuba is a **research-oriented physically based renderer**. The code is
organised as a small set of C++ shared libraries plus a large number of
**plugins**. Each plugin is a `.so` (or `.dll` / `.dylib`) implementing one
"thing" — a BSDF, an emitter, an integrator, a shape, a sampler, etc. The
renderer loads the plugins it needs at run-time based on the `type=…`
attributes in your scene XML.

This is exactly the extension model you want: to add a new BSDF you write
**one** `.cpp` file in `src/bsdfs/`, register it in the `SConscript`,
recompile, and you can immediately use `<bsdf type="myplugin"/>` in any scene.

## Top-level layout

```
mitsuba/
├── SConstruct              # Top-level SCons build orchestration
├── build/                  # SCons configuration / install scripts
├── config.py               # Compiler/Boost/Qt paths chosen during ./configure
├── data/                   # Test scenes, IOR tables, fonts, ...
├── dist/                   # **Build output** (binaries + plugins/)
├── doc/                    # The original LaTeX manual
├── include/mitsuba/        # Public headers — the API you code against
│   ├── core/               # Math, Spectrum, Properties, Stream, ...
│   ├── render/             # Scene/BSDF/Emitter/Integrator/Medium APIs
│   ├── bidir/              # Helpers for bidirectional methods (BDPT/MLT)
│   └── hw/                 # GPU/realtime preview support
├── setpath.sh              # Source this AFTER building to put dist/ on PATH
└── src/
    ├── libcore/            # libmitsuba-core.so
    ├── librender/          # libmitsuba-render.so
    ├── libbidir/           # libmitsuba-bidir.so
    ├── libhw/              # libmitsuba-hw.so (OpenGL preview)
    ├── libpython/          # Python bindings
    ├── mitsuba/            # The `mitsuba` CLI binary + `mtssrv`, `mtsutil`
    ├── mtsgui/             # The Qt GUI
    ├── bsdfs/              # ★ BSDF plugins   (diffuse, dielectric, coating, hk, …)
    ├── emitters/           # ★ Light source plugins (point, spot, area, envmap, sun…)
    ├── integrators/        # path/, bdpt/, mlt/, photonmapper/, vpl/, …
    ├── medium/             # homogeneous, heterogeneous
    ├── phase/              # hg, isotropic, rayleigh, microflake, mixturephase
    ├── samplers/           # ldsampler, halton, sobol, independent
    ├── shapes/             # sphere, rectangle, ply, obj, instance, …
    ├── sensors/            # perspective, orthographic, thinlens, …
    ├── films/              # hdrfilm, ldrfilm
    ├── rfilters/           # box, gaussian, mitchell, catmullrom, …
    ├── subsurface/         # dipole, …
    ├── textures/           # bitmap, checkerboard, gridtexture, …
    ├── volume/             # constvolume, gridvolume, …
    └── tests/              # Unit tests + chi-square sampling validation
```

## Library layering

```
   libmitsuba-core.so  ──► libmitsuba-render.so ──► libmitsuba-bidir.so
                                       │
                                       └─► all plugins (.so) link here
```

* **`libcore`** — no dependency on rendering. Provides:
  * `Spectrum`, `Vector`, `Point`, `Normal`, `Frame`, `AABB`
  * `Properties` (the typed key–value store filled from XML attributes)
  * `Stream`, `InstanceManager` (binary serialization for network rendering)
  * `Object`/`ref<>` reference-counted base
  * Plugin loader, file resolver, RNG, math helpers (`fresnelDielectricExt`,
    `warp::squareToCosineHemisphere`, …)

* **`librender`** — abstract interfaces:
  * `Scene`, `Sensor`, `Sampler`, `Film`, `Integrator`
  * `BSDF`, `Emitter`, `Texture`, `Medium`, `PhaseFunction`, `Shape`, `Subsurface`
  * Sampling records: `BSDFSamplingRecord`, `DirectSamplingRecord`,
    `PositionSamplingRecord`, `DirectionSamplingRecord`,
    `MediumSamplingRecord`, `PhaseFunctionSamplingRecord`,
    `RadianceQueryRecord`.

* **`libbidir`** — only used by BDPT / MLT / PSSMLT. Vertex / edge graph
  abstraction for path-space rendering. You can ignore it for the core task.

## Plugin model in 30 seconds

Every plugin has the same shape:

```cpp
class MyThing : public BSDF /* or Emitter, Integrator, … */ {
public:
    MyThing(const Properties &props) : BSDF(props) { /* read XML params */ }
    MyThing(Stream *s, InstanceManager *m) : BSDF(s, m) { /* deserialize */ }

    /* override the abstract methods of the base class */

    MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(MyThing, /*abstract=*/false, BSDF)
MTS_EXPORT_PLUGIN(MyThing, "Human-readable name");
```

`MTS_EXPORT_PLUGIN` (defined in `include/mitsuba/core/cobject.h`) expands to
two C-linkage symbols, `CreateInstance` and `GetDescription`, which is what
`Plugin::loadPlugin` looks up via `dlsym` when it sees `<bsdf type="mything"/>`.

## Build → run loop

1. `./configure.py` (once) — picks a `config-*.py` template into `config.py`.
2. `scons -j8` — compiles everything; outputs land in `dist/`.
3. `source setpath.sh` — prepends `dist/` to `PATH` and `LD_LIBRARY_PATH`,
   sets `MITSUBA_DIR`, etc.
4. `mitsuba scene.xml` — the CLI loads `dist/plugins/<type>.so` for every
   `<bsdf type="…"/>`, `<emitter type="…"/>`, `<integrator type="…"/>` etc.

When you add or modify a plugin, just `scons -j8` again — it only rebuilds
the changed `.so` and reinstalls it in `dist/plugins/`.

## What renders a frame, in one paragraph

The **scene** holds **sensors**, **emitters**, **shapes** (each shape carries
a **BSDF**, optionally a **Subsurface** integrator, optionally an interior
**Medium**), and one **integrator**. The integrator drives the render loop:
for every pixel it asks the **sampler** for random numbers, traces rays
through the scene, evaluates BSDFs and emitters along the path, and writes
sample contributions into the **film** through the **reconstruction filter**.
For `path`/`volpath` integrators, this is a recursive random walk with
next-event estimation (NEE) and multiple importance sampling (MIS) — that's
exactly what your LayeredBSDF needs to be friends with.

## Files you'll be reading most for the task

* `include/mitsuba/render/bsdf.h` — the BSDF interface contract.
* `include/mitsuba/render/emitter.h` — the Emitter interface contract.
* `src/bsdfs/diffuse.cpp` — the simplest working BSDF.
* `src/bsdfs/coating.cpp` — a non-trivial BSDF that **wraps** a nested BSDF
  with refraction and Beer–Lambert absorption. Closest existing analogue to a
  layered material.
* `src/bsdfs/hk.cpp` — Hanrahan–Krueger *single-scattering* layer (no random
  walk, but uses sigmaT / albedo / phase function — same building blocks).
* `src/integrators/path/path.cpp` — canonical NEE+MIS path tracer. Your
  BSDF must satisfy the contract this code assumes.
* `src/emitters/spot.cpp`, `src/emitters/directional.cpp` — closest
  templates for the secondary task.
