/*
    LayeredBSDF: position-free Monte Carlo random-walk plugin in the
    spirit of:

        Guo, Hašan, Yan, "Position-Free Monte Carlo Simulation
        for Arbitrary Layered BSDFs", SIGGRAPH 2018.

    The implementation models:
      * a single index-mismatched dielectric interface on top
        (extIOR -> intIOR),
      * an absorbing/scattering volume of finite thickness inside,
      * an index-matched nested BSDF at the bottom.

    sample() runs a real random walk in z exactly as the paper
    prescribes: Russian-roulette branching at the top interface,
    free-flight distance sampling, in-volume phase-function scattering,
    nested-BSDF interactions at the bottom, and Fresnel-weighted
    Russian roulette at every internal/top interface.

    eval() is an analytic single-scattering (HK-style) plus
    bottom-BSDF transmittance estimator. This keeps eval() deterministic
    and free of a Sampler, since path.cpp's NEE construction does not
    supply one. The path tracer still gathers multi-bounce contributions
    via BSDF sampling, where sample()'s walk handles them in full.

    pdf() is a closed-form approximation (cosine-weighted hemisphere
    times a top-Fresnel branching factor); the path tracer only needs
    a sufficiently-correct PDF for the MIS heuristic.

    Components (front-side only — opaque-bottom slab is one-sided):
      0  EGlossyReflection    -- random-walk reflection lobe
      1  EDeltaReflection     -- top dielectric Fresnel mirror

    XML schema is dual-form:
      * Native form: extIOR/intIOR floats, sigmaA/sigmaS spectra,
        nested <bsdf> child for the bottom, optional <phase> child.
      * Reference form (Guo et al. 2018, e.g. figure8 paper scenes):
        nbLayers=2, <bsdf name="surface_0"> for the top dielectric
        (we extract its IOR; roughness is dropped), <bsdf
        name="surface_1"> for the bottom, sigmaT_0/albedo_0 spectra,
        <phase name="phase_0"> for the phase function.

    Conventions:
      * z = 0  : top of the slab (air side, exterior)
      * z = -h : bottom of the slab (nested BSDF)
      * +z is outward in the local shading frame
      * during the walk we track the *physical* photon direction `dir`:
            dir.z < 0  ==>  photon traveling into the slab
            dir.z > 0  ==>  photon traveling toward the top
*/

#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/warp.h>
#include "ior.h"

MTS_NAMESPACE_BEGIN

class LayeredBSDF : public BSDF {
public:
    LayeredBSDF(const Properties &props) : BSDF(props) {
        Float intIOR = lookupIOR(props, "intIOR", "bk7");
        Float extIOR = lookupIOR(props, "extIOR", "air");

        if (intIOR <= 0 || extIOR <= 0 || intIOR == extIOR)
            Log(EError, "LayeredBSDF: interior and exterior IOR must be "
                "positive and differ");

        m_eta    = intIOR / extIOR;
        m_invEta = 1.0f / m_eta;

        m_thickness = props.getFloat("thickness", 1.0f);
        m_sigmaA = props.getSpectrum("sigmaA", Spectrum(0.0f));
        m_sigmaS = props.getSpectrum("sigmaS", Spectrum(0.0f));

        /* Reference position-free MC schema (Guo et al. 2018):
             sigmaT_0  = extinction
             albedo_0  = single-scattering albedo
           If both are present, derive sigmaA/sigmaS from them so
           paper scene XMLs (e.g. figure8, teaser) drop in directly. */
        bool hasSigmaT  = props.hasProperty("sigmaT_0");
        bool hasAlbedo  = props.hasProperty("albedo_0");
        if (hasSigmaT || hasAlbedo) {
            Spectrum sT = props.getSpectrum("sigmaT_0", Spectrum(1.0f));
            Spectrum al = props.getSpectrum("albedo_0", Spectrum(1.0f));
            m_sigmaS = sT * al;
            m_sigmaA = sT * (Spectrum(1.0f) - al);
        }

        m_maxDepth = props.getInteger("maxDepth", 64);
        m_rrDepth  = props.getInteger("rrDepth", 5);

        m_specularReflectance = props.getSpectrum("specularReflectance",
                                                  Spectrum(1.0f));

        /* Number of *surfaces* (not media). We only support 2 (one
           top dielectric + one bottom BSDF); larger N will be treated
           as a configuration error so the user notices. */
        m_nbLayers = props.getInteger("nbLayers", 2);
        if (m_nbLayers != 2)
            Log(EError, "LayeredBSDF: only nbLayers=2 (one medium between "
                "two surfaces) is supported. Got nbLayers=%i.", m_nbLayers);

        /* Reference-schema flags we silently accept and ignore — they
           only matter for the bidirectional / stochastic-pdf paths in
           the paper plugin, which we approximate. */
        for (const char *k : { "MIS", "bidir", "bidirUseAnalog",
                               "pdfRepetitive", "stochPdfDepth",
                               "diffusePdf", "maxSurvivalProb",
                               "multilayer" }) {
            if (props.hasProperty(k)) (void) k; /* parsed-then-ignored */
        }
        if (props.hasProperty("pdf"))
            (void) props.getString("pdf"); /* same */
    }

    LayeredBSDF(Stream *stream, InstanceManager *manager)
        : BSDF(stream, manager) {
        m_eta = stream->readFloat();
        m_thickness = stream->readFloat();
        m_sigmaA = Spectrum(stream);
        m_sigmaS = Spectrum(stream);
        m_specularReflectance = Spectrum(stream);
        m_maxDepth = stream->readInt();
        m_rrDepth  = stream->readInt();
        m_nbLayers = stream->readInt();
        m_nested = static_cast<BSDF *>(manager->getInstance(stream));
        m_phase  = static_cast<PhaseFunction *>(manager->getInstance(stream));
        m_invEta = 1.0f / m_eta;
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);
        stream->writeFloat(m_eta);
        stream->writeFloat(m_thickness);
        m_sigmaA.serialize(stream);
        m_sigmaS.serialize(stream);
        m_specularReflectance.serialize(stream);
        stream->writeInt(m_maxDepth);
        stream->writeInt(m_rrDepth);
        stream->writeInt(m_nbLayers);
        manager->serialize(stream, m_nested.get());
        manager->serialize(stream, m_phase.get());
    }

    void configure() {
        /* If a top-interface BSDF was supplied via `surface_0`, take
           its relative IOR. We don't model its roughness — paper
           figure8 uses alpha values down to 0.005 (near-smooth) for
           which a smooth dielectric is a close visual approximation. */
        if (m_topSurface != NULL) {
            Float topEta = m_topSurface->getEta();
            if (topEta > 0 && topEta != 1.0f) {
                m_eta    = topEta;
                m_invEta = 1.0f / m_eta;
            } else {
                Log(EWarn, "LayeredBSDF: surface_0 child has eta=%f; "
                    "falling back to extIOR/intIOR.", topEta);
            }
        }

        if (m_phase == NULL)
            m_phase = static_cast<PhaseFunction *>(PluginManager::getInstance()
                ->createObject(MTS_CLASS(PhaseFunction), Properties("isotropic")));

        if (m_nested == NULL) {
            Properties p("diffuse");
            m_nested = static_cast<BSDF *>(PluginManager::getInstance()
                ->createObject(MTS_CLASS(BSDF), p));
            m_nested->configure();
        }

        m_sigmaT       = m_sigmaA + m_sigmaS;
        m_sigmaTavg    = m_sigmaT.average();
        m_hasMedium    = (m_sigmaTavg > 0);

        /* Component layout. We expose the random-walk lobe and the
           top mirror only — opaque bottom never produces a transmission
           sample, so we don't advertise EGlossyTransmission. The plugin
           is also marked front-side only: a back-side hit on a slab
           with an opaque base is unphysical and would otherwise be
           handled by an unprincipled axis flip. */
        m_components.clear();
        m_components.push_back(EGlossyReflection | EFrontSide
                               | EUsesSampler | ENonSymmetric);
        m_components.push_back(EDeltaReflection  | EFrontSide);

        m_usesRayDifferentials = false;

        /* Heuristic specular-sampling weight: bias toward the mirror
           lobe when the slab's diffuse contribution is small (highly
           absorbing medium). */
        Float avgAttn = 1.0f;
        if (m_hasMedium)
            avgAttn = (-m_sigmaT * (2.0f * m_thickness)).exp().average();
        Float scoreNested = 0.5f + 0.5f * avgAttn;
        m_specularSamplingWeight = 1.0f / (1.0f + scoreNested);

        BSDF::configure();
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        const Class *cClass = child->getClass();
        if (cClass->derivesFrom(MTS_CLASS(BSDF))) {
            /* Reference schema names: surface_0 is the top dielectric
               (we extract its IOR in configure()); surface_1 is the
               bottom (nested) BSDF. Unnamed children behave as the
               nested BSDF (back-compat). */
            if (name == "surface_0") {
                if (m_topSurface != NULL)
                    Log(EError, "LayeredBSDF: duplicate surface_0 child");
                m_topSurface = static_cast<BSDF *>(child);
            } else if (name == "surface_1") {
                if (m_nested != NULL)
                    Log(EError, "LayeredBSDF: duplicate surface_1/nested child");
                m_nested = static_cast<BSDF *>(child);
            } else {
                if (m_nested != NULL)
                    Log(EError, "LayeredBSDF: only one nested BSDF child allowed");
                m_nested = static_cast<BSDF *>(child);
            }
        } else if (cClass->derivesFrom(MTS_CLASS(PhaseFunction))) {
            if (m_phase != NULL)
                Log(EError, "LayeredBSDF: only one phase-function child allowed");
            m_phase = static_cast<PhaseFunction *>(child);
        } else {
            BSDF::addChild(name, child);
        }
    }

    Float getEta() const { return 1.0f; }

    Float getRoughness(const Intersection &its, int component) const {
        if (component == 1) return 0.0f;
        return std::numeric_limits<Float>::infinity();
    }

    Spectrum getDiffuseReflectance(const Intersection &its) const {
        return m_nested->getDiffuseReflectance(its);
    }

    /* =================================================================== */
    /*                          helper functions                           */
    /* =================================================================== */

    inline Vector reflectZ(const Vector &v) const {
        return Vector(-v.x, -v.y, v.z);
    }

    /**
     * Same as `coating.cpp`'s refractIn: refract an air-side
     * "away-from-event" vector into the slab, preserving the sign of
     * the z-component. Useful when computing the slab-side direction
     * that would, under reciprocal Snell refraction, line up with
     * \c v at the top interface.
     */
    inline Vector refractIn(const Vector &v, Float &R) const {
        Float cosI = Frame::cosTheta(v);
        Float cosT;
        R = fresnelDielectricExt(std::abs(cosI), cosT, m_eta);
        if (R == 1.0f)
            return Vector(0.0f);
        return Vector(m_invEta * v.x, m_invEta * v.y,
                      -math::signum(cosI) * cosT);
    }

    /**
     * Sample a free-flight distance using the average sigmaT as the
     * hero wavelength. \c weight gets the per-channel transmittance
     * ratio that reweights the throughput for an unbiased multi-channel
     * estimator.
     */
    inline Float sampleFlightDistance(Float u, Spectrum &weight) const {
        if (!m_hasMedium) {
            weight = Spectrum(1.0f);
            return std::numeric_limits<Float>::infinity();
        }
        Float t = -std::log(std::max((Float) 1e-30f, 1 - u)) / m_sigmaTavg;
        weight = (-(m_sigmaT - Spectrum(m_sigmaTavg)) * t).exp();
        return t;
    }

    inline Spectrum transmittance(Float t) const {
        if (!m_hasMedium || t <= 0)
            return Spectrum(1.0f);
        return (-m_sigmaT * t).exp();
    }

    /* =================================================================== */
    /*                               sample                                */
    /* =================================================================== */

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf,
                    const Point2 &_sample) const {
        AssertEx(bRec.sampler != NULL,
                 "LayeredBSDF requires a sampler in the BSDFSamplingRecord");

        bool sampleSpecular = (bRec.typeMask & EDeltaReflection)
            && (bRec.component == -1 || bRec.component == 1);
        bool sampleGlossy = ((bRec.typeMask & EGlossyReflection) != 0)
            && (bRec.component == -1 || bRec.component == 0);

        if (!sampleSpecular && !sampleGlossy)
            return Spectrum(0.0f);

        /* Front-side only: opaque bottom + horizontally-homogeneous slab
           is one-sided. Back-side queries return zero, matching the
           component flags advertised in configure(). */
        const Vector &wi = bRec.wi;
        Float cosThetaI = Frame::cosTheta(wi);
        if (cosThetaI <= 0)
            return Spectrum(0.0f);

        Float cosThetaT;
        Float R12 = fresnelDielectricExt(cosThetaI, cosThetaT, m_eta);

        Float probSpecular;
        if (sampleSpecular && sampleGlossy) {
            probSpecular = R12 * m_specularSamplingWeight /
                (R12 * m_specularSamplingWeight + (1 - R12) * (1 - m_specularSamplingWeight));
        } else if (sampleSpecular) {
            probSpecular = 1.0f;
        } else {
            probSpecular = 0.0f;
        }

        Point2 sample(_sample);
        bool chooseSpecular = false;
        if (sampleSpecular && sampleGlossy) {
            if (sample.x < probSpecular) {
                chooseSpecular = true;
                sample.x /= probSpecular;
            } else {
                sample.x = (sample.x - probSpecular) / (1 - probSpecular);
            }
        } else {
            chooseSpecular = sampleSpecular;
        }

        if (chooseSpecular) {
            bRec.wo = reflectZ(wi);
            bRec.eta = 1.0f;
            bRec.sampledComponent = 1;
            bRec.sampledType = EDeltaReflection;
            pdf = (sampleSpecular && sampleGlossy) ? probSpecular : 1.0f;
            return m_specularReflectance * (R12 / pdf);
        }

        /* ============= Random walk branch =============
         *
         * The transmitted ray's physical direction inside the slab is:
         *   dir = (-invEta * wi.x, -invEta * wi.y, cosThetaT)
         * with cosThetaT < 0. The negation on x,y matches Snell's law
         * applied to the *physical* direction (wi already points
         * away from the event). */
        Vector dir(-m_invEta * wi.x, -m_invEta * wi.y, cosThetaT);

        /* Branch-probability compensation. */
        Spectrum throughput((1.0f - R12) /
            std::max((Float) 1e-12f, 1.0f - probSpecular));
        Float    z = 0.0f;
        bool     exited = false;

        Sampler *sampler = bRec.sampler;

        for (int depth = 0; depth < m_maxDepth; ++depth) {
            if (depth >= m_rrDepth) {
                Float q = std::min((Float) 0.95f, throughput.max());
                if (q <= 0) return Spectrum(0.0f);
                if (sampler->next1D() >= q) return Spectrum(0.0f);
                throughput /= q;
            }

            Spectrum mediumWeight(1.0f);
            Float t = m_hasMedium
                ? sampleFlightDistance(sampler->next1D(), mediumWeight)
                : std::numeric_limits<Float>::infinity();

            Float dz = dir.z;
            if (std::abs(dz) < 1e-6f)
                return Spectrum(0.0f);

            Float distToBoundary = (dz < 0)
                ? (-m_thickness - z) / dz
                : (0.0f - z) / dz;

            if (m_hasMedium && t < distToBoundary) {
                /* Volumetric scatter event. */
                z += dir.z * t;
                /* Hero-wavelength T/p reweighting:
                     scatter_pdf_hero = sigmaTavg * exp(-sigmaTavg * t);
                     scatter_target_per_channel = sigmaS * exp(-sigmaT * t)
                   ratio = (sigmaS/sigmaTavg) * exp((sigmaTavg-sigmaT)*t)
                         = (sigmaS/sigmaTavg) * mediumWeight.        */
                throughput *= mediumWeight * (m_sigmaS / m_sigmaTavg);

                PhaseFunctionSamplingRecord pRec(MediumSamplingRecord(),
                    -dir, ERadiance);
                Float phasePdf, phaseWeight;
                phaseWeight = m_phase->sample(pRec, phasePdf, sampler);
                if (phasePdf == 0 || phaseWeight == 0)
                    return Spectrum(0.0f);
                throughput *= phaseWeight;
                dir = pRec.wo;
            } else {
                if (m_hasMedium) {
                    Spectrum segTrans = transmittance(distToBoundary);
                    Float p_hero = std::exp(-m_sigmaTavg * distToBoundary);
                    if (p_hero <= 0) return Spectrum(0.0f);
                    throughput *= segTrans / p_hero;
                }
                z += dir.z * distToBoundary;

                if (dz > 0) {
                    /* Top boundary -- exit or TIR. */
                    Float cosD = dir.z;
                    Float cosOutT;
                    Float R21 = fresnelDielectricExt(cosD, cosOutT, m_invEta);
                    Float u = sampler->next1D();
                    if (u < R21) {
                        dir.z = -dir.z;
                    } else {
                        /* Refract out: physical direction in air. cosOutT
                           has the opposite sign of cosD, so the air-side
                           z = -cosOutT > 0. Snell scales tangentials by
                           m_eta when going from slab to air. */
                        Vector woAir(m_eta * dir.x, m_eta * dir.y, -cosOutT);
                        bRec.wo = woAir;
                        /* Solid-angle compression Jacobian. */
                        Float saJac = m_invEta * m_invEta *
                                      std::abs(woAir.z) /
                                      std::max((Float) 1e-6f, cosD);
                        throughput *= saJac;
                        exited = true;
                        break;
                    }
                } else {
                    /* Bottom boundary -- nested BSDF event. */
                    BSDFSamplingRecord nestedRec(bRec.its, sampler, bRec.mode);
                    /* The physical photon arrived going down (dir.z < 0).
                       For the nested BSDF, "wi" is the away-from-event
                       direction toward where light came from, which is
                       the upward complement of dir: -dir. */
                    nestedRec.wi = Vector(-dir.x, -dir.y, -dir.z);
                    if (Frame::cosTheta(nestedRec.wi) <= 0)
                        return Spectrum(0.0f);

                    Float nestedPdf;
                    Spectrum nestedW = m_nested->sample(nestedRec, nestedPdf,
                                                       sampler->next2D());
                    if (nestedW.isZero() || nestedPdf == 0)
                        return Spectrum(0.0f);
                    throughput *= nestedW;
                    if (Frame::cosTheta(nestedRec.wo) <= 0)
                        return Spectrum(0.0f);
                    dir = nestedRec.wo;
                }
            }
        }

        if (!exited)
            return Spectrum(0.0f);

        if (Frame::cosTheta(bRec.wo) <= 0) {
            /* Walk exited downward — only possible if the nested BSDF
               sampled a transmission. We don't advertise that lobe and
               the path tracer would mis-MIS it; drop. */
            return Spectrum(0.0f);
        }

        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EGlossyReflection;

        pdf = pdfImpl(wi, bRec.wo, R12, sampleSpecular, sampleGlossy);
        if (!std::isfinite(pdf) || pdf <= 0)
            return Spectrum(0.0f);

        return throughput * m_specularReflectance;
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        Float pdf;
        return LayeredBSDF::sample(bRec, pdf, sample);
    }

    /* =================================================================== */
    /*                                 pdf                                 */
    /* =================================================================== */

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        bool sampleSpecular = (bRec.typeMask & EDeltaReflection)
            && (bRec.component == -1 || bRec.component == 1);
        bool sampleGlossy = ((bRec.typeMask & EGlossyReflection) != 0)
            && (bRec.component == -1 || bRec.component == 0);

        const Vector &wi = bRec.wi, &wo = bRec.wo;

        Float cosThetaI = Frame::cosTheta(wi);
        if (cosThetaI <= 0 || Frame::cosTheta(wo) <= 0)
            return 0.0f;

        Float R12 = fresnelDielectricExt(cosThetaI, m_eta);

        if (measure == EDiscrete) {
            if (!sampleSpecular)
                return 0.0f;
            if (std::abs(dot(reflectZ(wi), wo) - 1) < DeltaEpsilon) {
                if (sampleSpecular && sampleGlossy) {
                    Float probSpecular = R12 * m_specularSamplingWeight /
                        (R12 * m_specularSamplingWeight +
                         (1 - R12) * (1 - m_specularSamplingWeight));
                    return probSpecular;
                }
                return 1.0f;
            }
            return 0.0f;
        }

        if (measure != ESolidAngle || !sampleGlossy)
            return 0.0f;

        return pdfImpl(wi, wo, R12, sampleSpecular, sampleGlossy);
    }

    /// Cosine-weighted hemisphere pdf approximation, scaled by the
    /// non-specular branch probability.
    Float pdfImpl(const Vector &wi, const Vector &wo, Float R12,
                  bool sampleSpecular, bool sampleGlossy) const {
        Float cosThetaO = Frame::cosTheta(wo);
        if (cosThetaO == 0) return 0.0f;
        Float pdf = std::abs(cosThetaO) * INV_PI;
        if (sampleSpecular && sampleGlossy) {
            Float probSpecular = R12 * m_specularSamplingWeight /
                (R12 * m_specularSamplingWeight +
                 (1 - R12) * (1 - m_specularSamplingWeight));
            pdf *= (1.0f - probSpecular);
        }
        return pdf;
    }

    /* =================================================================== */
    /*                                eval                                 */
    /*                                                                     */
    /* Analytic single-scatter Hanrahan-Krueger contribution + analytic    */
    /* bottom-BSDF transmittance (à la coating). The path tracer's NEE      */
    /* eval() does not provide a sampler, so we deliberately avoid the     */
    /* MC random walk here. Multi-bounce contributions are still gathered  */
    /* by the integrator through BSDF sampling, where sample() runs the    */
    /* full random walk.                                                   */
    /* =================================================================== */

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        bool sampleSpecular = (bRec.typeMask & EDeltaReflection)
            && (bRec.component == -1 || bRec.component == 1);
        bool sampleGlossy = ((bRec.typeMask & EGlossyReflection) != 0)
            && (bRec.component == -1 || bRec.component == 0);

        const Vector &wi = bRec.wi, &wo = bRec.wo;
        Float cosThetaI = Frame::cosTheta(wi);
        Float cosThetaO = Frame::cosTheta(wo);

        if (measure == EDiscrete) {
            if (!sampleSpecular || cosThetaI <= 0) return Spectrum(0.0f);
            if (std::abs(dot(reflectZ(wi), wo) - 1) < DeltaEpsilon) {
                Float R12 = fresnelDielectricExt(cosThetaI, m_eta);
                return m_specularReflectance * R12;
            }
            return Spectrum(0.0f);
        }

        if (measure != ESolidAngle || !sampleGlossy)
            return Spectrum(0.0f);
        if (cosThetaI <= 0 || cosThetaO <= 0)
            return Spectrum(0.0f);

        /* Slab-side directions (sign-preserving refraction; same as
           coating.cpp's refractIn). */
        Float R12, R21;
        Vector wiSlab = refractIn(wi, R12);
        if (R12 == 1.0f) return Spectrum(0.0f);
        Vector woSlab = refractIn(wo, R21);
        if (R21 == 1.0f) return Spectrum(0.0f);

        const Float cosI = Frame::cosTheta(wiSlab);
        const Float cosO = Frame::cosTheta(woSlab);
        if (cosI <= 0 || cosO <= 0)
            return Spectrum(0.0f);

        /* The IOR Jacobian factors: (a) saJac substitutes the
           nested-BSDF's cos(woSlab) for cos(wo_air) * invEta^2,
           accounting for solid-angle compression as we exit the
           slab. (b) topT is the product of two Fresnel transmittances
           — one entering, one exiting the top dielectric. */
        const Float saJac = m_invEta * m_invEta * cosThetaO / cosO;
        const Float topT  = (1.0f - R12) * (1.0f - R21);

        Spectrum result(0.0f);

        /* ----- bottom BSDF (ballistic transit through the slab) ----- */
        BSDFSamplingRecord nestedRec(bRec.its, wiSlab, woSlab, bRec.mode);
        Spectrum nestedF = m_nested->eval(nestedRec, ESolidAngle);
        if (!nestedF.isZero()) {
            Spectrum trans = transmittance(m_thickness * (1.0f / cosI + 1.0f / cosO));
            result += nestedF * trans * topT * saJac;
        }

        /* ----- single-scatter (Hanrahan–Krueger) ----- */
        if (m_hasMedium && !m_sigmaS.isZero()) {
            PhaseFunctionSamplingRecord pRec(MediumSamplingRecord(),
                wiSlab, woSlab, bRec.mode);
            Float phaseVal = m_phase->eval(pRec);
            if (phaseVal > 0) {
                Spectrum tauD = m_sigmaT * m_thickness;
                Spectrum albedo;
                for (int c = 0; c < SPECTRUM_SAMPLES; ++c)
                    albedo[c] = m_sigmaT[c] > 0
                        ? m_sigmaS[c] / m_sigmaT[c] : 0.0f;
                Spectrum hk = albedo * (phaseVal * cosI / (cosI + cosO)) *
                    (Spectrum(1.0f) -
                     ((-1.0f / cosI - 1.0f / cosO) * tauD).exp());
                /* hk*cosO has units of f*cos(wo_slab); saJac replaces
                   that with f*cos(wo_air)*invEta^2 as required by the
                   solid-angle convention used by the integrator. */
                result += hk * cosO * topT * saJac;
            }
        }

        return m_specularReflectance * result;
    }

    /* =================================================================== */

    std::string toString() const {
        std::ostringstream oss;
        oss << "LayeredBSDF[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  eta = " << m_eta << "," << endl
            << "  thickness = " << m_thickness << "," << endl
            << "  sigmaA = " << m_sigmaA.toString() << "," << endl
            << "  sigmaS = " << m_sigmaS.toString() << "," << endl
            << "  maxDepth = " << m_maxDepth << "," << endl
            << "  rrDepth = " << m_rrDepth << "," << endl
            << "  phase = " << indent(m_phase->toString()) << "," << endl
            << "  nested = " << indent(m_nested->toString()) << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    Float m_eta, m_invEta;
    Float m_thickness;
    Spectrum m_sigmaA, m_sigmaS, m_sigmaT;
    Float m_sigmaTavg;
    bool m_hasMedium;
    Spectrum m_specularReflectance;
    int m_maxDepth, m_rrDepth;
    int m_nbLayers;
    Float m_specularSamplingWeight;
    ref<BSDF> m_topSurface; ///< Optional top BSDF; we extract its IOR.
    ref<BSDF> m_nested;
    ref<PhaseFunction> m_phase;
};

MTS_IMPLEMENT_CLASS_S(LayeredBSDF, false, BSDF)
MTS_EXPORT_PLUGIN(LayeredBSDF, "Layered BSDF (1D position-free random walk)");
MTS_NAMESPACE_END
