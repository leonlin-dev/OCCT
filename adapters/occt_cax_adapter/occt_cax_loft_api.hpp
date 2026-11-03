// ===========================================================================
// adapters/occt_cax_adapter/occt_cax_loft_api.hpp
//
// High-level native wrapper for Open CASCADE Technology (OCCT):
// occt_cax::loft_row_coplanar(...) receives native OCCT curve handles
// (Handle(Geom_BSplineCurve) or occt_cax::OcctCurve) and returns the native
// OCCT surface handle (Handle(Geom_BSplineSurface)) produced by the Module-1..5
// pipeline loft::CoplanarClosedFormLoft<OcctCaxLoftContext>::build_loft_surface.
// ===========================================================================
#pragma once

#include "occt_cax_loft_context.hpp"

#include <loft/loft_context_concepts.hpp>
#include <loft/coplanar_closed_form_loft.hpp>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace occt_cax {

using CurveHandle   = OcctCaxLoftContext::NurbsEngine::CurvT;
using SurfaceHandle = OcctCaxLoftContext::NurbsEngine::SurfT;

// Request knobs (pipeline pass-throughs).
struct LoftRequest {
    bool is_closed_loft = false;
    loft::SectionMode mode = loft::SectionMode::RowCoplanarBlend;
    int degree_u = 3;
    int degree_v = 3;
    loft::TolerancePolicy<double> tol{};
    std::size_t min_num_cp_v = 16;
    std::vector<double> anchors_u;
};

// Surface + audit payload for callers that need the invariant report.
struct LoftRowResult {
    SurfaceHandle surface;
    std::size_t num_cp_u = 0;
    std::size_t num_cp_v = 0;
    loft::LoftAuditReport<double> audit;

    double max_contact_deviation = 0.0;
    double lebesgue_v = 0.0;
    bool structural_ok = true;
    bool any_escalated = false;
};

namespace detail {

inline bool is_unit_domain(const CurveHandle& c) noexcept {
    if (c.knots.empty()) return false;
    const double lo = c.knots.front();
    const double hi = c.knots.back();
    return std::fabs(lo) <= 1e-14 && std::fabs(hi - 1.0) <= 1e-14;
}

inline void require_clamped(const CurveHandle& c, const char* what) {
    const int p = c.degree;
    if (static_cast<int>(c.knots.size()) <
        static_cast<int>(c.cp_homo.size()) + p + 1) {
        throw std::invalid_argument(
            std::string(what) +
            ": knot-vector/control-point size mismatch (n = num_kt - p - 1)");
    }
    const std::size_t p_sz = static_cast<std::size_t>(p);
    for (std::size_t i = 0; i <= p_sz; ++i) {
        if (c.knots[i] != c.knots.front() ||
            c.knots[c.knots.size() - 1 - i] != c.knots.back()) {
            throw std::invalid_argument(
                std::string(what) +
                ": knot vector must be clamped ((p+1)-fold end knots)");
        }
    }
}

inline void require_clamped_unit(const CurveHandle& c, const char* what) {
    require_clamped(c, what);
    if (!is_unit_domain(c)) {
        throw std::invalid_argument(
            std::string(what) + ": normalized knot domain must be [0, 1]");
    }
}

inline CurveHandle normalize_to_unit(const CurveHandle& c) {
    return is_unit_domain(c) ? c
                             : OcctCaxLoftContext::NurbsEngine::reparameterize_linear(
                                   c, 0.0, 1.0);
}

} // namespace detail

inline LoftRowResult loft_row_coplanar(
    const std::vector<CurveHandle>& sections,
    const std::vector<CurveHandle>& guides,
    const CurveHandle* spine = nullptr,
    const LoftRequest& request = LoftRequest{}) {
    if (sections.size() < 2) {
        throw std::invalid_argument(
            "occt_cax::loft_row_coplanar: need at least 2 profile sections");
    }

    std::vector<loft::ProfileStationInput<CurveHandle>> stations;
    stations.reserve(sections.size());
    for (std::size_t ell = 0; ell < sections.size(); ++ell) {
        detail::require_clamped(sections[ell], "section");
        loft::ProfileStationInput<CurveHandle> st;
        st.base_curve = detail::normalize_to_unit(sections[ell]);
        detail::require_clamped_unit(st.base_curve, "normalized section");
        stations.push_back(std::move(st));
    }

    std::vector<CurveHandle> unit_guides;
    unit_guides.reserve(guides.size());
    for (std::size_t k = 0; k < guides.size(); ++k) {
        detail::require_clamped(guides[k], "guide");
        unit_guides.push_back(detail::normalize_to_unit(guides[k]));
        detail::require_clamped_unit(unit_guides.back(), "normalized guide");
    }

    CurveHandle unit_spine;
    const CurveHandle* spine_ptr = nullptr;
    if (spine != nullptr) {
        detail::require_clamped(*spine, "spine");
        unit_spine = detail::normalize_to_unit(*spine);
        detail::require_clamped_unit(unit_spine, "normalized spine");
        spine_ptr = &unit_spine;
    }

    const auto out =
        loft::CoplanarClosedFormLoft<OcctCaxLoftContext>::build_loft_surface(
            loft::Span<const loft::ProfileStationInput<CurveHandle>>(
                stations.data(), stations.size()),
            loft::Span<const CurveHandle>(unit_guides.data(),
                                          unit_guides.size()),
            spine_ptr,
            request.is_closed_loft,
            request.mode,
            request.degree_u,
            request.degree_v,
            request.tol);

    LoftRowResult res;
    res.surface = out.surface;
    res.num_cp_u = out.num_cp_u;
    res.num_cp_v = out.num_cp_v;
    res.audit = out.audit;
    res.structural_ok = true;
    res.any_escalated = false;
    return res;
}

inline LoftRowResult loft_row_coplanar(
    const std::vector<Handle(Geom_BSplineCurve)>& sections,
    const std::vector<Handle(Geom_BSplineCurve)>& guides,
    const Handle(Geom_BSplineCurve)* spine = nullptr,
    const LoftRequest& request = LoftRequest{}) {
    std::vector<CurveHandle> sec_handles;
    sec_handles.reserve(sections.size());
    for (const auto& s : sections) sec_handles.emplace_back(s);

    std::vector<CurveHandle> gd_handles;
    gd_handles.reserve(guides.size());
    for (const auto& g : guides) gd_handles.emplace_back(g);

    if (spine != nullptr && !spine->IsNull()) {
        CurveHandle sp_handle(*spine);
        return loft_row_coplanar(sec_handles, gd_handles, &sp_handle, request);
    }
    return loft_row_coplanar(sec_handles, gd_handles, nullptr, request);
}

} // namespace occt_cax
