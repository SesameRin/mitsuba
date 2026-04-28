# 2. The BSDF system

This is the most important chapter for your task — your `LayeredBSDF` will
plug into this same interface.

Reference files: `include/mitsuba/render/bsdf.h`, `src/bsdfs/diffuse.cpp`,
`src/bsdfs/coating.cpp`, `src/bsdfs/hk.cpp`.

## 2.1 What a BSDF must implement

A subclass of `BSDF` overrides three pure virtual methods:

```cpp
virtual Spectrum sample(BSDFSamplingRecord &bRec,
                        const Point2 &sample) const = 0;

virtual Spectrum sample(BSDFSamplingRecord &bRec,
                        Float &pdf,
                        const Point2 &sample) const = 0;

virtual Spectrum eval(const BSDFSamplingRecord &bRec,
                      EMeasure measure = ESolidAngle) const = 0;

virtual Float    pdf (const BSDFSamplingRecord &bRec,
                      EMeasure measure = ESolidAngle) const = 0;
```

In addition, in `configure()` (called after construction) you populate
`m_components` (one entry per "lobe", each a bitmask of `EBSDFType` flags)
and let the base class compute `m_combinedType`.

### `BSDFSamplingRecord` (the query record)

Contains:
* `wi`, `wo` — the incident and outgoing **direction in the local shading
  frame**, both pointing **away** from the scattering event. (This is
  Mitsuba's convention; if you instead think of `wi` as "incoming light",
  remember that `cosTheta(wi) > 0` means light is hitting the front side.)
* `its` — the surface intersection (gives you UVs, world-space frame, etc.).
* `eta` — relative IOR along the sampled direction (will be set by sample()).
* `mode` — `ERadiance` or `EImportance` — needed for non-symmetric BSDFs
  (transmission across a refractive boundary scales differently for radiance
  vs. importance transport — this is the "non-symmetric scattering" issue
  in Veach's thesis).
* `typeMask`, `component` — let the integrator (or another wrapping BSDF)
  request only specific components.
* `sampledType`, `sampledComponent` — the BSDF must set these on success
  inside `sample()` (e.g. `EDiffuseReflection`, `EDeltaTransmission`).

### `EMeasure`

* `ESolidAngle` — smooth lobes (most common).
* `EDiscrete` — Dirac-delta lobes (mirrors, ideal refraction). For these,
  `eval` returns the spectral attenuation (no cosine, no PDF density), and
  `pdf` returns a *probability* (not a density).
* `ELength` — 1D delta lobes (rare, used by hair-like models).

### Component flags (`EBSDFType` — `bsdf.h:224`)

```
ENull, EDiffuseReflection, EDiffuseTransmission,
EGlossyReflection, EGlossyTransmission,
EDeltaReflection, EDeltaTransmission,
EDelta1DReflection, EDelta1DTransmission,
+ attribute flags: EAnisotropic, ESpatiallyVarying, ENonSymmetric,
                   EFrontSide, EBackSide, EUsesSampler.
```

For your layered BSDF you will probably want
`EGlossyReflection | EGlossyTransmission | EFrontSide | EBackSide |
ENonSymmetric | EUsesSampler` — the last bit is **important**: it tells
the renderer that you need extra random numbers from `bRec.sampler` beyond
the two passed to `sample()`. Random walks through layers absolutely need
this.

## 2.2 Conventions you will hear from the path tracer

The path integrator calls a BSDF roughly like this (see `path.cpp`):

```cpp
// (1) NEE: pick a direction toward a sampled emitter, then ask the BSDF
//     what its value (× cos) and PDF are.
BSDFSamplingRecord bRec(its, its.toLocal(emitterDir), ERadiance);
Spectrum f = bsdf->eval(bRec);            // returns f(wi,wo) * cos(theta_o)
Float    p = bsdf->pdf (bRec);            // density w.r.t. solid angle

// (2) BSDF sampling: ask the BSDF to choose a direction.
BSDFSamplingRecord bRec(its, sampler, ERadiance);
Float pdf;
Spectrum weight = bsdf->sample(bRec, pdf, sample2D);
// weight already equals f(wi,wo) * cos(theta_o) / pdf.
```

So the contract is:

* **`eval()` returns `f * cos(theta_o)` for smooth lobes.** Forgetting the
  cosine is the #1 mistake; the integrator does *not* multiply by it.
* **`sample()` returns `f * cos(theta_o) / pdf`.** This lets a clever BSDF
  exploit cancellation (e.g. in the diffuse case, `f * cos / pdf = albedo`).
* **`pdf()` is a solid-angle density**, except for delta lobes (then it's a
  probability and the corresponding `eval` already factored in the delta).

These rules apply to *every* lobe of *every* BSDF, including yours.

## 2.3 Reading `diffuse.cpp` — the simplest example

```cpp
class SmoothDiffuse : public BSDF {
    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure m) const {
        if (!(bRec.typeMask & EDiffuseReflection) || m != ESolidAngle ||
            cosTheta(bRec.wi) <= 0 || cosTheta(bRec.wo) <= 0)
            return Spectrum(0.0f);
        return reflectance(bRec.its) * (INV_PI * cosTheta(bRec.wo));
    }
    Float pdf(const BSDFSamplingRecord &bRec, EMeasure m) const {
        // returns cos(theta_o) / pi
        return warp::squareToCosineHemispherePdf(bRec.wo);
    }
    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &s) const {
        bRec.wo = warp::squareToCosineHemisphere(s);
        bRec.eta = 1.0f;
        bRec.sampledType = EDiffuseReflection;
        pdf = warp::squareToCosineHemispherePdf(bRec.wo);
        return reflectance(bRec.its);   // = (rho/pi)*cos / (cos/pi)
    }
};
```

Things to note:

1. The `(typeMask & ...)` early-out: if the integrator asked only for a
   different lobe, return zero.
2. The `cosTheta(wi) <= 0` guard means "ignore back-facing geometry" — this
   BSDF is one-sided. For your layered material, you'll want it
   **two-sided**, with the layer stack symmetric across `z = 0`.
3. `bRec.eta = 1.0f` because nothing refracts here. For a transmission
   event you'd write `bRec.eta = (entering ? eta : 1/eta)`.

## 2.4 Reading `coating.cpp` — most relevant analogue

`coating.cpp` is a **wrapping BSDF**: it adds a smooth dielectric layer on
top of a nested BSDF. This is morally a 2-layer LayeredBSDF, hand-coded.
The structure is exactly what you'll generalize.

Key ideas from `coating.cpp` you should steal:

* **Component fan-out.** It exposes the nested BSDF's components *plus* one
  extra `EDeltaReflection` component for the top-surface mirror reflection
  (`coating.cpp:165-171`). Your LayeredBSDF will expose at least one
  `EGlossy*` component (the layered scattering) plus the top-interface
  Fresnel reflection.

* **Refract in / refract out / Beer–Lambert.**
  ```cpp
  Vector wiPrime = refractIn(bRec.wi, R12);          // cos sign flips
  Vector woPrime = refractIn(bRec.wo, R21);
  result = nested->eval(...) * (1-R12)*(1-R21);
  result *= exp(-sigmaA * thickness *
                (1/|cos(wiPrime)| + 1/|cos(woPrime)|));
  result *= invEta*invEta * cosTheta(wo)/cosTheta(woPrime);  // solid-angle compression
  ```
  The 1/cos factors are the path-length-through-slab factors — your random
  walk will produce these naturally instead of analytically.

* **Importance sampling weights between top and nested.**
  `coating.cpp:175-182` computes `m_specularSamplingWeight` from the average
  absorption, then in `sample()` does a Russian-roulette-like split:
  ```cpp
  Float probSpecular = R12 * w / (R12*w + (1-R12)*(1-w));
  if (sample.x < probSpecular) /* specular branch */ else /* nested branch */
  ```
  You will do something similar at the **top interface** of your stack: with
  probability `R12` reflect specularly, with probability `1-R12` enter the
  random walk. (Inside the walk, you do the same at every layer interface.)

* **MIS-friendliness.** `coating.cpp` always returns a **valid PDF** for both
  branches it could have produced — even the specular one (a discrete prob.).
  The path integrator queries `pdf()` separately to compute MIS weights when
  the BSDF-sampled ray happens to hit a light. You must support this:
  given an arbitrary `(wi, wo)` pair, return the density that your sampling
  routine *would* have put there. For a random-walk sampler this is
  delicate (see §6 — the trick is to account for both reflection and
  transmission branches, not the full tree of internal scatter events).

## 2.5 Reading `hk.cpp` — closest physics analogue

`hk.cpp` (Hanrahan–Krueger) handles a single-scattering, finite-thickness
slab on top of a Lambertian. It already deals with `sigmaT`, `albedo`,
`thickness`, and a `PhaseFunction` child. You'll see:

```cpp
sigmaT = sigmaA + sigmaS;
tauD   = sigmaT * thickness;       // optical depth through the slab
albedo = sigmaS / sigmaT;
```

The big difference from your task: HK uses **closed-form single-scattering**
+ a diffusion term, not a Monte Carlo random walk. So this file is useful
to (a) see how a slab BSDF takes a `<phase>` child and `<bsdf>` child, and
(b) to validate your random walk: in the optically-thin limit your numerical
estimator should reproduce HK's analytic single-scattering term.

## 2.6 The `Frame` and `warp` helpers

* `Frame::cosTheta(v)`, `sinTheta`, `tanTheta`, `cosPhi`, `sinPhi` — all
  cheap because vectors are in the local frame `(s, t, n)`.
* `bRec.its.toLocal(worldDir)` and `its.toWorld(localDir)` convert between
  shading frame and world space.
* `warp::squareToCosineHemisphere(u)` (cosine-weighted), `squareToUniformCone`,
  `squareToUniformSphere`, `squareToBeckmann`, `squareToGGX` — all in
  `include/mitsuba/core/warp.h`.

## 2.7 Energy conservation, two-sidedness, validation

* The base class has `ensureEnergyConservation` helpers that clamp/scale
  textures whose maximum exceeds 1.
* For two-sided BSDFs, the convention in `bsdf.h:258-260` is to set both
  `EFrontSide` and `EBackSide` flags. The integrator passes whichever side
  the ray actually hit.
* **Chi-square test.** `src/tests/test_chisquare.cpp` plus
  `data/tests/test_bsdf.xml` runs any BSDF through a statistical test that
  fits a histogram of `sample()` outputs against `pdf()`. This is the
  cleanest way to convince yourself that your `pdf` matches your `sample`.
  Add an entry for your plugin to `test_bsdf.xml`.

## 2.8 Skeleton you'll start from

```cpp
class LayeredBSDF : public BSDF {
public:
    LayeredBSDF(const Properties &props) : BSDF(props) {
        // read top/bottom IOR, list of layers (etas, sigmaA, sigmaS, phase, thickness),
        // optional bottom BSDF, optional max walk depth, etc.
    }

    void configure() {
        m_components.clear();
        m_components.push_back(EGlossyReflection   | EFrontSide|EBackSide);
        m_components.push_back(EGlossyTransmission | EFrontSide|EBackSide);
        // Optional: an explicit specular component for the top interface
        m_components.push_back(EDeltaReflection    | EFrontSide|EBackSide);
        BSDF::configure();
    }

    Spectrum eval (const BSDFSamplingRecord &bRec, EMeasure m) const override;
    Float    pdf  (const BSDFSamplingRecord &bRec, EMeasure m) const override;
    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &s) const override;
    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &s) const override;

    void addChild(const std::string &name, ConfigurableObject *child) override;
    void serialize(Stream*, InstanceManager*) const override;

    MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(LayeredBSDF, false, BSDF)
MTS_EXPORT_PLUGIN(LayeredBSDF, "Layered BSDF (1D random walk)");
```

The interesting work is concentrated in `eval`, `pdf`, and `sample`. See
`06-layered-bsdf-task.md` for how to fill those in.
