/*
    Grazing beam emitter -- a delta-position spot light with a sharp
    cosine-power angular falloff, designed to highlight grazing-angle
    Fresnel response on layered materials.

    The emitter sits at the origin of its local frame and points along
    +z. The radiant intensity along a unit direction d (in local frame)
    is

        I(d) = intensity * max(0, cos(theta))^exponent,
              where cos(theta) = d.z

    capped to zero beyond the cutoff angle. exponent controls beam
    tightness; exponent=1 reduces to a cosine spot, while exponent>>1
    produces a near-collimated pencil. The plugin computes the
    correctly-normalised PDF (cosine-power within a cone) so that
    sampleRay / sampleDirection / pdfDirection are all consistent.

    Like the built-in spot emitter, this emitter is a delta-position
    light: sampleDirect returns a discrete probability of 1 and the
    pre-divided contribution intensity * falloff / dist^2.
*/

#include <mitsuba/render/scene.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

class GrazingEmitter : public Emitter {
public:
    GrazingEmitter(const Properties &props) : Emitter(props) {
        m_intensity = props.getSpectrum("intensity", Spectrum(1.0f));
        Float cutoffDeg = props.getFloat("cutoffAngle", 30.0f);
        m_cutoffAngle = degToRad(cutoffDeg);
        m_exponent = std::max((Float) 0.0f, props.getFloat("exponent", 32.0f));
        m_type = EDeltaPosition;
    }

    GrazingEmitter(Stream *stream, InstanceManager *manager)
        : Emitter(stream, manager) {
        m_intensity = Spectrum(stream);
        m_cutoffAngle = stream->readFloat();
        m_exponent = stream->readFloat();
        configure();
    }

    void configure() {
        Emitter::configure();
        m_cosCutoff = std::cos(m_cutoffAngle);

        /*
         * Normalisation: integrate cos^n(theta) sin(theta) dtheta dphi over
         * the cone theta in [0, cutoff]:
         *
         *      integral_0^{2pi} d phi
         *      integral_0^{cutoff} cos^n(theta) sin(theta) dtheta
         *      = 2 pi * (1 - cos^{n+1}(cutoff)) / (n + 1)
         *
         * The angular pdf inside the cone is
         *      pdf(d) = (n+1) / (2 pi (1 - cos^{n+1}(cutoff))) * cos^n(theta)
         */
        m_pdfNormalisation = (m_exponent + 1.0f) /
            (2.0f * M_PI *
             (1.0f - std::pow((Float) m_cosCutoff, (Float)(m_exponent + 1.0f))));

        /* Total integrated solid-angle weight, used for the radiant power
           returned by samplePosition. */
        m_omegaTotal = 2.0f * M_PI *
            (1.0f - std::pow((Float) m_cosCutoff, (Float)(m_exponent + 1.0f))) /
            (m_exponent + 1.0f);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Emitter::serialize(stream, manager);
        m_intensity.serialize(stream);
        stream->writeFloat(m_cutoffAngle);
        stream->writeFloat(m_exponent);
    }

    /* ------------------------------------------------------------------ */
    /*                       angular profile helpers                      */
    /* ------------------------------------------------------------------ */

    /// Local-frame angular intensity profile (no normalisation).
    inline Float falloff(const Vector &d) const {
        if (d.z <= m_cosCutoff)
            return 0.0f;
        return std::pow(d.z, m_exponent);
    }

    /// Sample a local-frame direction d with cos^n falloff inside the cone.
    inline Vector sampleConeDir(const Point2 &sample) const {
        Float u = sample.x;
        Float oneMinusBase = 1.0f - std::pow((Float) m_cosCutoff,
                                             (Float)(m_exponent + 1.0f));
        Float cosTheta = std::pow(1.0f - u * oneMinusBase,
                                  (Float)(1.0f / (m_exponent + 1.0f)));
        Float sinTheta = math::safe_sqrt(1.0f - cosTheta * cosTheta);
        Float sinPhi, cosPhi;
        math::sincos(2.0f * M_PI * sample.y, &sinPhi, &cosPhi);
        return Vector(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);
    }

    /// Solid-angle pdf of \ref sampleConeDir inside the cone.
    inline Float pdfConeDir(const Vector &d) const {
        if (d.z <= m_cosCutoff) return 0.0f;
        return m_pdfNormalisation * std::pow(d.z, m_exponent);
    }

    /* ------------------------------------------------------------------ */
    /*                        Emitter API methods                         */
    /* ------------------------------------------------------------------ */

    Spectrum samplePosition(PositionSamplingRecord &pRec, const Point2 &sample,
            const Point2 *extra) const {
        const Transform &trafo = m_worldTransform->eval(pRec.time);
        pRec.p = trafo.transformAffine(Point(0.0f));
        pRec.n = Normal(0.0f);
        pRec.pdf = 1.0f;
        pRec.measure = EDiscrete;
        /*
         * The emitted-power conversion mirrors what spot.cpp does (return
         * intensity * 4 pi as the position term). Here we instead return
         * intensity * Omega_total, since the angular profile is already
         * integrated against this weight in evalDirection() / sampleDirection().
         */
        return m_intensity * m_omegaTotal;
    }

    Spectrum evalPosition(const PositionSamplingRecord &pRec) const {
        return (pRec.measure == EDiscrete) ?
            (m_intensity * m_omegaTotal) : Spectrum(0.0f);
    }

    Float pdfPosition(const PositionSamplingRecord &pRec) const {
        return (pRec.measure == EDiscrete) ? 1.0f : 0.0f;
    }

    Spectrum sampleDirection(DirectionSamplingRecord &dRec,
            PositionSamplingRecord &pRec, const Point2 &sample,
            const Point2 *extra) const {
        const Transform &trafo = m_worldTransform->eval(pRec.time);
        Vector localDir = sampleConeDir(sample);
        dRec.d = trafo(localDir);
        dRec.pdf = pdfConeDir(localDir);
        dRec.measure = ESolidAngle;
        return Spectrum(1.0f / m_omegaTotal);
    }

    Float pdfDirection(const DirectionSamplingRecord &dRec,
            const PositionSamplingRecord &pRec) const {
        if (dRec.measure != ESolidAngle) return 0.0f;
        const Transform &trafo = m_worldTransform->eval(pRec.time);
        Vector local = trafo.inverse()(dRec.d);
        return pdfConeDir(local);
    }

    Spectrum evalDirection(const DirectionSamplingRecord &dRec,
            const PositionSamplingRecord &pRec) const {
        if (dRec.measure != ESolidAngle) return Spectrum(0.0f);
        const Transform &trafo = m_worldTransform->eval(pRec.time);
        Vector local = trafo.inverse()(dRec.d);
        return Spectrum(falloff(local) / m_omegaTotal);
    }

    Spectrum sampleRay(Ray &ray, const Point2 &spatialSample,
            const Point2 &directionalSample, Float time) const {
        const Transform &trafo = m_worldTransform->eval(time);
        Vector localDir = sampleConeDir(directionalSample);
        ray.setTime(time);
        ray.setOrigin(trafo.transformAffine(Point(0.0f)));
        ray.setDirection(trafo(localDir));
        Float dirPdf = pdfConeDir(localDir);
        if (dirPdf == 0)
            return Spectrum(0.0f);
        return m_intensity * (falloff(localDir) / dirPdf);
    }

    Spectrum sampleDirect(DirectSamplingRecord &dRec, const Point2 &sample) const {
        const Transform &trafo = m_worldTransform->eval(dRec.time);
        dRec.p = trafo.transformAffine(Point(0.0f));
        dRec.uv = Point2(0.5f);
        dRec.n = Normal(0.0f);
        dRec.d = dRec.p - dRec.ref;
        Float distSqr = dRec.d.lengthSquared();
        dRec.dist = std::sqrt(distSqr);
        dRec.d /= dRec.dist;
        dRec.pdf = 1.0f;
        dRec.measure = EDiscrete;
        Vector localDir = trafo.inverse()(-dRec.d);
        return m_intensity * falloff(localDir) / distSqr;
    }

    Float pdfDirect(const DirectSamplingRecord &dRec) const {
        return dRec.measure == EDiscrete ? 1.0f : 0.0f;
    }

    AABB getAABB() const {
        return m_worldTransform->getTranslationBounds();
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "GrazingEmitter[" << endl
            << "  intensity = " << m_intensity.toString() << "," << endl
            << "  cutoffAngle = " << (m_cutoffAngle * 180.0f / M_PI) << "," << endl
            << "  exponent = " << m_exponent << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    Spectrum m_intensity;
    Float m_cutoffAngle, m_cosCutoff;
    Float m_exponent;
    Float m_pdfNormalisation, m_omegaTotal;
};

MTS_IMPLEMENT_CLASS_S(GrazingEmitter, false, Emitter)
MTS_EXPORT_PLUGIN(GrazingEmitter, "Grazing-angle spot emitter "
                  "(cosine-power falloff)");
MTS_NAMESPACE_END
