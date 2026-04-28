# 5. Build and plugin system

## 5.1 SCons in 60 seconds

Mitsuba uses **SCons** (a Python-based build system). The relevant files:

* `SConstruct` (project root) — orchestrates everything.
* `build/SConscript.configure` — applies the chosen `config.py` to set
  compiler flags, includes, linker options.
* `build/SConscript.install` — copies built artifacts into `dist/`.
* `src/<module>/SConscript` — declares what to build for each module.
* `config.py` — written from one of the templates in `build/` by your
  initial `./configure`.

For a plugin, the `SConscript` is a one-liner. For example,
`src/bsdfs/SConscript`:

```python
plugins += env.SharedLibrary('diffuse', ['diffuse.cpp'])
plugins += env.SharedLibrary('coating', ['coating.cpp'])
…
```

`SharedLibrary('foo', ['foo.cpp'])` produces `dist/plugins/foo.so` (the
`build/SConscript.install` step copies it after linking).

## 5.2 What `MTS_EXPORT_PLUGIN` actually does

From `include/mitsuba/core/cobject.h`:

```cpp
#define MTS_EXPORT_PLUGIN(name, descr) \
    extern "C" { \
        void MTS_EXPORT *CreateInstance(const Properties &props) { \
            return new name(props); \
        } \
        const char MTS_EXPORT *GetDescription() { \
            return descr; \
        } \
    }
```

So every plugin `.so` exposes exactly two `extern "C"` symbols. When the
scene XML reads `<bsdf type="diffuse">`, Mitsuba's `Plugin` loader
`dlopen`s `dist/plugins/diffuse.so`, calls `dlsym("CreateInstance")`, and
hands it the `Properties` parsed from XML. The instance gets RTTI-checked
against `MTS_CLASS(BSDF)` (registered by `MTS_IMPLEMENT_CLASS_S`), and
that's it.

This means: **the file name is the plugin name in XML.** If you create
`src/bsdfs/layered.cpp` and add it as `env.SharedLibrary('layered',
['layered.cpp'])`, your scene uses `<bsdf type="layered"/>`.

`MTS_IMPLEMENT_CLASS` vs `MTS_IMPLEMENT_CLASS_S`:
* `_S` = also serialisable (provides `Stream`-based deserialisation
  constructor). Use this for anything that goes across the network for
  distributed renders, which is essentially everything you'd write.

## 5.3 Run-time loading: `setpath.sh`

`setpath.sh` (sourced, not executed) sets:
* `PATH` += `dist/`  (so `mitsuba`, `mtsutil`, `mtssrv`, `mtsgui` are found)
* `LD_LIBRARY_PATH` += `dist/`  (so the four `libmitsuba-*.so`s resolve)
* `MITSUBA_DIR` = `dist/` (used by the file resolver to find
  `dist/plugins/foo.so`)

If `mitsuba` says "Encountered an error: Plugin 'layered' could not be
loaded", 99% of the time the culprit is forgetting to re-source
`setpath.sh` or building into a different `dist/`.

## 5.4 Adding a new plugin — checklist

1. Drop `layered.cpp` into `src/bsdfs/`.
2. Edit `src/bsdfs/SConscript` and append:
   ```python
   plugins += env.SharedLibrary('layered', ['layered.cpp'])
   ```
3. Build: `scons -j8`.
4. Verify the new file: `ls dist/plugins/layered.so`.
5. Reference it from a scene: `<bsdf type="layered" …>`.
6. Render: `mitsuba scene.xml -o out.exr`.

For a new emitter you do the same in `src/emitters/SConscript`. The
plugin name conventions match what's already there (`spot`, `directional`,
`area`, `point`, …).

## 5.5 Useful CLI flags

* `mitsuba -v` (or `-vv`) — verbose / very-verbose log output. Most BSDFs
  print themselves to the log via `toString()`; very useful when debugging.
* `mitsuba -s scene.xml` — render then exit. (Default behavior, but worth
  knowing about.)
* `mitsuba -p N` — N rendering threads (override sampler choice).
* `mitsuba -D key=value` — define a parameter substitutable in scene XML
  via `$key`. Handy for parameter sweeps.
* `mtsutil` — runner for utilities (`mtsutil testcases`, `mtsutil
  blackbody`, `mtsutil joinrgb`, …). Test your plugin with chi-square here.

## 5.6 IDE / editor tips

* The build emits `compile_commands.json`-friendly output if you set
  `CXX = "clang++"` in `config.py`. Without that, use `bear -- scons -j8`
  or VS Code's C/C++ extension with manual include paths from `config.py`.
* Includes you'll touch most:
  ```
  -Iinclude
  -Iinclude/mitsuba                       (rare; use mitsuba/...)
  -I/usr/include/eigen3                   (system)
  -I<boost prefix>/include                (system)
  ```

## 5.7 Recompiling fast

When iterating on `layered.cpp`, the only file that gets recompiled is
that one + the link of `dist/plugins/layered.so`. A full `scons -j8` after
a one-line change should take <5 s. If you change a public header in
`include/mitsuba/render/`, the world rebuilds — try to keep edits inside
your plugin file and a private header `src/bsdfs/layered.h` until you're
sure you need a public surface.
