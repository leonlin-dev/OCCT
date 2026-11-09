// ===========================================================================
// adapters/occt_cax_adapter/occt_cax_loft_context.hpp
//
// Four-Pillar LoftingContext adapter binding Open CASCADE Technology (OCCT 7.8.1)
// native geometry (Geom_BSplineCurve, Geom_BSplineSurface, gp_Pnt, gp_Vec)
// to the Module-1 contract (include/loft/loft_context_concepts.hpp).
// ===========================================================================
#pragma once

#include <loft/loft_context_concepts.hpp>
#include <loft/coplanar_closed_form_loft.hpp>

#include <Standard_Handle.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <TColStd_Array2OfReal.hxx>

#include <loft_adapter_eigen_algebra.hpp>
#include <loft_adapter_point_projection.hpp>
#include <loft_adapter_surface_knot_insert.hpp>
#include <loft_adapter_curve_reverse.hpp>
#include <loft_adapter_curve_topology.hpp>
#include <loft_adapter_curve_concat.hpp>
#include <loft_adapter_curve_extremum.hpp>
#include <loft_adapter_guide_candidate_selection.hpp>
#include <loft_adapter_unified_basis.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace occt_cax {

using loft::Span;
using loft::MatrixView2D;

namespace detail {

inline double clampd(double v, double lo, double hi) noexcept {
    return v < lo ? lo : (v > hi ? hi : v);
}

inline void compress_knots(const std::vector<double>& flat_knots,
                           std::vector<double>& unique_knots,
                           std::vector<int>& mults,
                           double tol = 1e-12) {
    unique_knots.clear();
    mults.clear();
    if (flat_knots.empty()) return;
    unique_knots.push_back(flat_knots[0]);
    mults.push_back(1);
    for (std::size_t i = 1; i < flat_knots.size(); ++i) {
        if (std::fabs(flat_knots[i] - unique_knots.back()) <= tol) {
            mults.back()++;
        } else {
            unique_knots.push_back(flat_knots[i]);
            mults.push_back(1);
        }
    }
}

// A2.1: largest s with kv[s] <= t, clamped to [p, size-1-p-1]; t at the
// right end closes the last domain span (end-point convention).
inline std::size_t find_knot_span(int p, const std::vector<double>& kv,
                                  double t) noexcept {
    assert(kv.size() >= static_cast<std::size_t>(p) + 2);
    const std::size_t last = kv.size() - 1;
    const std::size_t hi0 = last - static_cast<std::size_t>(p) - 1;
    if (t >= kv[last]) return hi0;
    std::size_t lo = static_cast<std::size_t>(p);
    std::size_t hi = last;
    while (hi - lo > 1) {
        const std::size_t mid = (lo + hi) / 2;
        if (t < kv[mid]) hi = mid; else lo = mid;
    }
    return lo;
}

// de Boor value triangle + derivative recursion (The NURBS Book A2.3):
// ders_out is (max_d + 1) x (p + 1), row = derivative order.
inline void eval_basis_ders(int p, const std::vector<double>& kv,
                            std::size_t span, double t, int max_d,
                            std::vector<double>& ders_out) {
    const std::size_t pp = static_cast<std::size_t>(p) + 1;
    ders_out.assign((static_cast<std::size_t>(max_d) + 1) * pp, 0.0);
    if (t < kv[static_cast<std::size_t>(p)])
        t = kv[static_cast<std::size_t>(p)];
    const double t_end = kv[kv.size() - 1];
    if (t > t_end) t = t_end;
    std::vector<double> val(pp * pp, 0.0);
    const bool at_right_end =
        (span + static_cast<std::size_t>(p) + 1 == kv.size() - 1);
    val[0] = (t >= kv[span] &&
              (t < kv[span + 1] || (at_right_end && t == kv[span + 1])))
                 ? 1.0 : 0.0;
    for (int d = 1; d <= p; ++d) {
        const std::size_t dd = static_cast<std::size_t>(d);
        for (std::size_t j = 0; j <= dd; ++j) {
            const std::size_t g = span - dd + j;
            double acc = 0.0;
            if (j >= 1) {
                const double den = kv[g + dd] - kv[g];
                if (den > 0.0)
                    acc += (t - kv[g]) / den * val[(dd - 1) * pp + (j - 1)];
            }
            if (j <= dd - 1) {
                const double den = kv[g + 1 + dd] - kv[g + 1];
                if (den > 0.0)
                    acc += (kv[g + 1 + dd] - t) / den * val[(dd - 1) * pp + j];
            }
            val[dd * pp + j] = acc;
        }
    }
    for (std::size_t j = 0; j < pp; ++j) ders_out[j] = val[pp * pp - pp + j];

    std::vector<double> prev(val), next(pp * pp, 0.0);
    for (int k = 1; k <= max_d; ++k) {
        for (std::size_t idx = 0; idx < pp * pp; ++idx) next[idx] = 0.0;
        for (int d = 1; d <= p; ++d) {
            const std::size_t dd = static_cast<std::size_t>(d);
            for (std::size_t j = 0; j <= dd; ++j) {
                const std::size_t g = span - dd + j;
                double acc = 0.0;
                if (j >= 1) {
                    const double den = kv[g + dd] - kv[g];
                    if (den > 0.0) acc += prev[(dd - 1) * pp + (j - 1)] / den;
                }
                {
                    const double den = kv[g + 1 + dd] - kv[g + 1];
                    if (den > 0.0) acc -= prev[(dd - 1) * pp + j] / den;
                }
                next[dd * pp + j] = static_cast<double>(d) * acc;
            }
        }
        prev = next;
        for (std::size_t j = 0; j < pp; ++j)
            ders_out[static_cast<std::size_t>(k) * pp + j] =
                next[pp * pp - pp + j];
    }
}

inline void boehm_insert(int p, std::vector<double>& kv,
                         std::vector<std::array<double, 4>>& c, double x) {
    const std::size_t span = find_knot_span(p, kv, x);
    std::vector<std::array<double, 4>> q(c.size() + 1);
    for (std::size_t j = 0; j <= span - static_cast<std::size_t>(p); ++j)
        q[j] = c[j];
    for (std::size_t j = span - static_cast<std::size_t>(p) + 1; j <= span;
         ++j) {
        const double den = kv[j + static_cast<std::size_t>(p)] - kv[j];
        const double alpha = (den > 0.0) ? (x - kv[j]) / den : 0.0;
        for (int w = 0; w < 4; ++w) {
            q[j][static_cast<std::size_t>(w)] =
                alpha * c[j][static_cast<std::size_t>(w)] +
                (1.0 - alpha) * c[j - 1][static_cast<std::size_t>(w)];
        }
    }
    for (std::size_t j = span + 1; j < q.size(); ++j) q[j] = c[j - 1];
    kv.insert(kv.begin() + static_cast<std::ptrdiff_t>(span) + 1, x);
    c.swap(q);
}

inline void refine_knots(int p, std::vector<double>& kv,
                         std::vector<std::array<double, 4>>& c,
                         const std::vector<double>& kv_target) {
    std::size_t ti = 0;
    while (ti < kv_target.size()) {
        std::size_t tj = ti;
        while (tj < kv_target.size() && kv_target[tj] == kv_target[ti]) ++tj;
        const double x = kv_target[ti];
        std::size_t mult = 0;
        for (const double k : kv) { if (k == x) ++mult; }
        for (std::size_t rep = mult; rep < tj - ti; ++rep)
            boehm_insert(p, kv, c, x);
        ti = tj;
    }
}

} // namespace detail

// ===========================================================================
// OcctCurve — value-semantic OCCT curve wrapper.
//
// Embeds Handle(Geom_BSplineCurve) and mirrors degree, knots, cp_homo
// duck-typed by the loft module layer.
// ===========================================================================
struct OcctCurve {
    Handle(Geom_BSplineCurve) native;

    int degree = 3;
    std::vector<double> knots;
    std::vector<std::array<double, 4>> cp_homo;

    OcctCurve() = default;

    OcctCurve(const Handle(Geom_BSplineCurve)& c) : native(c) {
        sync_views();
    }

    static OcctCurve build(int p,
                           const std::vector<double>& kv,
                           const std::vector<std::array<double, 4>>& cp) {
        std::vector<double> uknots;
        std::vector<int> umults;
        detail::compress_knots(kv, uknots, umults);

        Standard_Integer nb_poles = static_cast<Standard_Integer>(cp.size());
        Standard_Integer nb_knots = static_cast<Standard_Integer>(uknots.size());

        TColgp_Array1OfPnt poles(1, nb_poles);
        TColStd_Array1OfReal weights(1, nb_poles);
        bool is_rational = false;

        for (Standard_Integer i = 1; i <= nb_poles; ++i) {
            const auto& h = cp[static_cast<std::size_t>(i - 1)];
            const Standard_Real wt = h[3];
            if (std::fabs(wt - 1.0) > 1e-12) is_rational = true;
            weights.SetValue(i, wt);
            if (std::fabs(wt) > 1e-300) {
                poles.SetValue(i, gp_Pnt(h[0] / wt, h[1] / wt, h[2] / wt));
            } else {
                poles.SetValue(i, gp_Pnt(h[0], h[1], h[2]));
            }
        }

        TColStd_Array1OfReal knots_arr(1, nb_knots);
        TColStd_Array1OfInteger mults_arr(1, nb_knots);
        for (Standard_Integer i = 1; i <= nb_knots; ++i) {
            knots_arr.SetValue(i, uknots[static_cast<std::size_t>(i - 1)]);
            mults_arr.SetValue(i, umults[static_cast<std::size_t>(i - 1)]);
        }

        Handle(Geom_BSplineCurve) handle;
        if (is_rational) {
            handle = new Geom_BSplineCurve(poles, weights, knots_arr, mults_arr, p);
        } else {
            handle = new Geom_BSplineCurve(poles, knots_arr, mults_arr, p);
        }
        return OcctCurve(handle);
    }

    void sync_views() {
        if (native.IsNull()) {
            degree = 0;
            knots.clear();
            cp_homo.clear();
            return;
        }

        degree = native->Degree();

        // Decompress knots
        const TColStd_Array1OfReal& uk = native->Knots();
        const TColStd_Array1OfInteger& um = native->Multiplicities();
        knots.clear();
        for (Standard_Integer i = uk.Lower(); i <= uk.Upper(); ++i) {
            const Standard_Real k = uk.Value(i);
            const Standard_Integer m = um.Value(i);
            for (Standard_Integer r = 0; r < m; ++r) {
                knots.push_back(k);
            }
        }

        // Homogeneous CPs
        const Standard_Integer nb_poles = native->NbPoles();
        cp_homo.clear();
        cp_homo.reserve(static_cast<std::size_t>(nb_poles));
        const TColStd_Array1OfReal* w = native->Weights();
        for (Standard_Integer i = 1; i <= nb_poles; ++i) {
            const gp_Pnt p = native->Pole(i);
            const Standard_Real wt = (w != nullptr) ? w->Value(i) : 1.0;
            cp_homo.push_back({{wt * p.X(), wt * p.Y(), wt * p.Z(), wt}});
        }
    }

    Handle(Geom_BSplineCurve) to_native() const {
        return native;
    }
};

struct OcctCaxLoftContext {
    // -----------------------------------------------------------------------
    // Pillar 1: Kernel — scalars & basic Euclidean/homogeneous primitives.
    // -----------------------------------------------------------------------
    struct Kernel {
        using FT       = double;
        using Vec2T    = std::array<FT, 2>;
        using Pt3T     = std::array<FT, 3>;
        using Vec3T    = std::array<FT, 3>;
        using PtHomo3T = std::array<FT, 4>;  // [w*x, w*y, w*z, w]

        static FT dot2(const Vec2T& a, const Vec2T& b) noexcept {
            return a[0] * b[0] + a[1] * b[1];
        }
        static FT dot3(const Vec3T& a, const Vec3T& b) noexcept {
            return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        }
        static Vec3T cross3(const Vec3T& a, const Vec3T& b) noexcept {
            return {{ a[1] * b[2] - a[2] * b[1],
                      a[2] * b[0] - a[0] * b[2],
                      a[0] * b[1] - a[1] * b[0] }};
        }
        static FT norm3(const Vec3T& v) noexcept {
            return std::sqrt(dot3(v, v));
        }
        static Vec3T normalize3(const Vec3T& v) noexcept {
            const FT n = norm3(v);
            if (n == FT(0)) return {{0, 0, 0}};
            return {{ v[0] / n, v[1] / n, v[2] / n }};
        }
        static Vec3T diff3(const Pt3T& a, const Pt3T& b) noexcept {
            return {{ a[0] - b[0], a[1] - b[1], a[2] - b[2] }};
        }
        static Pt3T add_scaled3(const Pt3T& p, const Vec3T& v, FT s) noexcept {
            return {{ p[0] + s * v[0], p[1] + s * v[1], p[2] + s * v[2] }};
        }
    };

    // -----------------------------------------------------------------------
    // Pillar 2: NurbsEngine — native OCCT handles + shared unified basis math.
    // -----------------------------------------------------------------------
    struct NurbsEngine {
        using FT    = double;
        using CurvT = OcctCurve;
        using SurfT = Handle(Geom_BSplineSurface);

        static int degree(const CurvT& c) noexcept { return c.degree; }

        static std::vector<FT> knots(const CurvT& c) { return c.knots; }

        static std::vector<Kernel::PtHomo3T> control_points_homo(const CurvT& c) {
            return c.cp_homo;
        }

        static std::size_t find_knot_span(int p, Span<const FT> kv, FT t) {
            const std::vector<FT> kvs(kv.begin(), kv.end());
            return detail::find_knot_span(p, kvs, t);
        }

        static void eval_basis_ders(int p, Span<const FT> kv, std::size_t span,
                                    FT t, int max_d,
                                    MatrixView2D<FT> ders_out) {
            const std::vector<FT> kvs(kv.begin(), kv.end());
            std::vector<FT> flat;
            detail::eval_basis_ders(p, kvs, span, t, max_d, flat);
            const std::size_t pp = static_cast<std::size_t>(p) + 1;
            for (std::size_t k = 0; k <= static_cast<std::size_t>(max_d); ++k)
                for (std::size_t j = 0; j < pp; ++j)
                    ders_out(k, j) = flat[k * pp + j];
        }

        static void compute_oslo_matrix(int p, Span<const FT> kv_old,
                                        Span<const FT> kv_new,
                                        MatrixView2D<FT> M_out) {
            const std::vector<FT> ko(kv_old.begin(), kv_old.end());
            const std::vector<FT> kn(kv_new.begin(), kv_new.end());
            const std::size_t n_old = ko.size() - static_cast<std::size_t>(p) - 1;
            const std::size_t n_new = kn.size() - static_cast<std::size_t>(p) - 1;
            for (std::size_t a = 0; a < n_old; ++a) {
                std::vector<double> row_kv = ko;
                std::vector<std::array<double, 4>> row(
                    n_old, std::array<double, 4>{{0, 0, 0, 0}});
                row[a] = {{1, 1, 1, 1}};
                detail::refine_knots(p, row_kv, row, kn);
                for (std::size_t b = 0; b < n_new; ++b)
                    M_out(a, b) = row[b][0];
            }
        }

        static CurvT create_nurbs_curve(int p, Span<const FT> kv,
                                        Span<const Kernel::PtHomo3T> cp_homo) {
            std::vector<std::array<double, 4>> cp;
            cp.reserve(cp_homo.size());
            for (const auto& h : cp_homo) cp.push_back(h);
            return OcctCurve::build(p, std::vector<double>(kv.begin(), kv.end()), cp);
        }

        static CurvT reparameterize_linear(const CurvT& c, FT a, FT b) {
            const std::size_t m = c.knots.size();
            const FT tmin = c.knots.front();
            const FT tmax = c.knots.back();
            std::vector<double> new_knots(m);
            std::vector<std::array<double, 4>> new_cps(c.cp_homo.size());

            if (b >= a) {
                for (std::size_t i = 0; i < m; ++i)
                    new_knots[i] = a + (b - a) * (c.knots[i] - tmin) / (tmax - tmin);
                new_cps = c.cp_homo;
            } else {
                for (std::size_t i = 0; i < m; ++i)
                    new_knots[i] = a + (b - a) * (c.knots[m - 1 - i] - tmin) / (tmax - tmin);
                const std::size_t ncp = c.cp_homo.size();
                for (std::size_t i = 0; i < ncp; ++i)
                    new_cps[i] = c.cp_homo[ncp - 1 - i];
            }
            return OcctCurve::build(c.degree, new_knots, new_cps);
        }

        // Shared adapter-side reversal body (operator-sinking review,
        // Task 2.2): the former microkernel copy, hosted by the lab.
        static CurvT reverse_curve(const CurvT& c) {
            return loft_adapter::reverse_curve_via<NurbsEngine, Kernel::PtHomo3T, FT>(c);
        }

        // Shared adapter-side topology bodies (operator-sinking review,
        // Task 2.2 remainder), hosted by the algorithm lab.
        static CurvT clip_curve(const CurvT& c, FT t0, FT t1, FT tol) {
            if (t0 > t1) return reverse_curve(clip_curve(c, t1, t0, tol));
            std::vector<double> kv = NurbsEngine::knots(c);
            std::vector<Kernel::PtHomo3T> cp_h = NurbsEngine::control_points_homo(c);
            if (kv.size() < 2 * (static_cast<std::size_t>(NurbsEngine::degree(c)) + 1))
                return c;
            std::vector<double> kv_out;
            std::vector<std::array<double, 4>> cp_out;
            loft_adapter::topo_clip_arrays(NurbsEngine::degree(c), std::move(kv),
                                           std::move(cp_h), static_cast<double>(t0),
                                           static_cast<double>(t1), static_cast<double>(tol),
                                           kv_out, cp_out);
            std::vector<Kernel::PtHomo3T> homo(cp_out.begin(), cp_out.end());
            return create_nurbs_curve(NurbsEngine::degree(c),
                                      Span<const FT>(kv_out.data(), kv_out.size()),
                                      Span<const Kernel::PtHomo3T>(homo.data(), homo.size()));
        }
        static CurvT convert_periodic_to_clamped(int p, Span<const FT> break_knots,
                                                 Span<const Kernel::PtHomo3T> cp_homo_per) {
            std::vector<double> kv_out;
            std::vector<std::array<double, 4>> cp_out;
            std::vector<std::array<double, 4>> per(cp_homo_per.begin(), cp_homo_per.end());
            loft_adapter::topo_periodic_to_clamped_arrays(
                p, std::vector<double>(break_knots.begin(), break_knots.end()), per, kv_out, cp_out);
            std::vector<Kernel::PtHomo3T> homo(cp_out.begin(), cp_out.end());
            return create_nurbs_curve(p, Span<const FT>(kv_out.data(), kv_out.size()),
                                      Span<const Kernel::PtHomo3T>(homo.data(), homo.size()));
        }

        static CurvT elevate_degree(const CurvT& c, int target_degree) {
            const int p = c.degree;
            if (target_degree <= p) return c;
            const auto res = loft_adapter::unified_degree_elevate(
                p, target_degree, c.knots, c.cp_homo);
            return OcctCurve::build(target_degree, res.first, res.second);
        }

        static SurfT create_nurbs_surface(int p_u, int p_v,
                                          Span<const FT> kv_u,
                                          Span<const FT> kv_v,
                                          Span<const Kernel::PtHomo3T> cp_grid_u_major) {
            const std::size_t n_u = kv_u.size() - static_cast<std::size_t>(p_u) - 1;
            const std::size_t n_v = kv_v.size() - static_cast<std::size_t>(p_v) - 1;
            if (cp_grid_u_major.size() != n_u * n_v) {
                throw std::invalid_argument(
                    "occt_cax::NurbsEngine::create_nurbs_surface: cp grid size mismatch");
            }

            std::vector<double> uknots, vknots;
            std::vector<int> umults, vmults;
            detail::compress_knots(std::vector<double>(kv_u.begin(), kv_u.end()), uknots, umults);
            detail::compress_knots(std::vector<double>(kv_v.begin(), kv_v.end()), vknots, vmults);

            Standard_Integer nb_u_knots = static_cast<Standard_Integer>(uknots.size());
            Standard_Integer nb_v_knots = static_cast<Standard_Integer>(vknots.size());
            TColStd_Array1OfReal u_knots_arr(1, nb_u_knots);
            TColStd_Array1OfInteger u_mults_arr(1, nb_u_knots);
            for (Standard_Integer i = 1; i <= nb_u_knots; ++i) {
                u_knots_arr.SetValue(i, uknots[static_cast<std::size_t>(i - 1)]);
                u_mults_arr.SetValue(i, umults[static_cast<std::size_t>(i - 1)]);
            }

            TColStd_Array1OfReal v_knots_arr(1, nb_v_knots);
            TColStd_Array1OfInteger v_mults_arr(1, nb_v_knots);
            for (Standard_Integer i = 1; i <= nb_v_knots; ++i) {
                v_knots_arr.SetValue(i, vknots[static_cast<std::size_t>(i - 1)]);
                v_mults_arr.SetValue(i, vmults[static_cast<std::size_t>(i - 1)]);
            }

            TColgp_Array2OfPnt poles(1, static_cast<Standard_Integer>(n_u),
                                     1, static_cast<Standard_Integer>(n_v));
            TColStd_Array2OfReal weights(1, static_cast<Standard_Integer>(n_u),
                                         1, static_cast<Standard_Integer>(n_v));
            bool is_rational = false;

            for (std::size_t j = 0; j < n_u; ++j) {
                for (std::size_t i = 0; i < n_v; ++i) {
                    const auto& h = cp_grid_u_major[j * n_v + i];
                    const Standard_Real wt = h[3];
                    if (std::fabs(wt - 1.0) > 1e-12) is_rational = true;
                    weights.SetValue(static_cast<Standard_Integer>(j + 1),
                                     static_cast<Standard_Integer>(i + 1), wt);
                    if (std::fabs(wt) > 1e-300) {
                        poles.SetValue(static_cast<Standard_Integer>(j + 1),
                                       static_cast<Standard_Integer>(i + 1),
                                       gp_Pnt(h[0] / wt, h[1] / wt, h[2] / wt));
                    } else {
                        poles.SetValue(static_cast<Standard_Integer>(j + 1),
                                       static_cast<Standard_Integer>(i + 1),
                                       gp_Pnt(h[0], h[1], h[2]));
                    }
                }
            }

            Handle(Geom_BSplineSurface) surf;
            if (is_rational) {
                surf = new Geom_BSplineSurface(poles, weights, u_knots_arr, v_knots_arr,
                                               u_mults_arr, v_mults_arr, p_u, p_v);
            } else {
                surf = new Geom_BSplineSurface(poles, u_knots_arr, v_knots_arr,
                                               u_mults_arr, v_mults_arr, p_u, p_v);
            }
            return surf;
        }

                // Shared adapter-side Boehm insertion body (operator-sinking
        // review, Task 2.1) through the OCCT representation: expand the
        // compressed knot sequences, de-homogenize poles x weights, insert,
        // and rebuild via the OCCT bridge.
        static void insert_surface_knots(SurfT& s, int uv_dir, FT t, int mult) {
            if (mult <= 0 || (uv_dir != 1 && uv_dir != 2)) return;
            if (s.IsNull()) return;
            const Geom_BSplineSurface* sf = s.get();
            const int p_u = sf->UDegree();
            const int p_v = sf->VDegree();
            std::vector<double> ku, kv;
            for (Standard_Integer i = sf->UKnots().Lower(); i <= sf->UKnots().Upper(); ++i)
                for (Standard_Integer m = 0; m < sf->UMultiplicity(i); ++m)
                    ku.push_back(sf->UKnot(i));
            for (Standard_Integer i = sf->VKnots().Lower(); i <= sf->VKnots().Upper(); ++i)
                for (Standard_Integer m = 0; m < sf->VMultiplicity(i); ++m)
                    kv.push_back(sf->VKnot(i));
            const std::size_t n_u = ku.size() - static_cast<std::size_t>(p_u) - 1;
            const std::size_t n_v = kv.size() - static_cast<std::size_t>(p_v) - 1;
            std::vector<double> cp(n_u * n_v * 4);
            for (Standard_Integer j = sf->Poles().LowerCol(); j <= sf->Poles().UpperCol(); ++j)
                for (Standard_Integer i = sf->Poles().LowerRow(); i <= sf->Poles().UpperRow(); ++i) {
                    const gp_Pnt& pnt = sf->Pole(i, j);
                    const double w = sf->Weight(i, j);
                    const std::size_t ju = static_cast<std::size_t>(i - sf->Poles().LowerRow());
                    const std::size_t iv = static_cast<std::size_t>(j - sf->Poles().LowerCol());
                    const std::size_t o = (ju * n_v + iv) * 4;
                    cp[o + 0] = pnt.X() * w; cp[o + 1] = pnt.Y() * w;
                    cp[o + 2] = pnt.Z() * w; cp[o + 3] = w;
                }
            loft_adapter::insert_surface_knots(ku, kv, cp, p_u, p_v, 4, uv_dir,
                                               static_cast<double>(t), mult);
            const std::size_t nn_v = kv.size() - static_cast<std::size_t>(p_v) - 1;
            const std::size_t nn_u = ku.size() - static_cast<std::size_t>(p_u) - 1;
            std::vector<Kernel::PtHomo3T> homo(nn_u * nn_v);
            for (std::size_t i = 0; i < nn_u * nn_v; ++i)
                for (int w = 0; w < 4; ++w)
                    homo[i][static_cast<std::size_t>(w)] = cp[i * 4 + static_cast<std::size_t>(w)];
            s = create_nurbs_surface(p_u, p_v,
                                     Span<const FT>(ku.data(), ku.size()),
                                     Span<const FT>(kv.data(), kv.size()),
                                     Span<const Kernel::PtHomo3T>(homo.data(), homo.size()));
        }
    };

    // -----------------------------------------------------------------------
    // Pillar 3: GeomEngine — single-entity differential geometry queries.
    // -----------------------------------------------------------------------
    struct GeomEngine {
        using FT    = double;
        using Pt3T  = std::array<FT, 3>;
        using Vec3T = std::array<FT, 3>;
        using CurvT = OcctCurve;

        static void eval_curve_ders(const CurvT& c, FT t, int max_d,
                                    Span<Vec3T> ders_out) {
            assert(max_d >= 0);
            assert(ders_out.size() >= static_cast<std::size_t>(max_d) + 1);
            if (c.native.IsNull()) return;

            const FT lo = c.native->FirstParameter();
            const FT hi = c.native->LastParameter();
            const FT tc = detail::clampd(t, lo, hi);

            if (max_d == 0) {
                gp_Pnt p;
                c.native->D0(tc, p);
                ders_out[0] = {{p.X(), p.Y(), p.Z()}};
            } else if (max_d == 1) {
                gp_Pnt p;
                gp_Vec v1;
                c.native->D1(tc, p, v1);
                ders_out[0] = {{p.X(), p.Y(), p.Z()}};
                ders_out[1] = {{v1.X(), v1.Y(), v1.Z()}};
            } else if (max_d == 2) {
                gp_Pnt p;
                gp_Vec v1, v2;
                c.native->D2(tc, p, v1, v2);
                ders_out[0] = {{p.X(), p.Y(), p.Z()}};
                ders_out[1] = {{v1.X(), v1.Y(), v1.Z()}};
                ders_out[2] = {{v2.X(), v2.Y(), v2.Z()}};
            } else {
                gp_Pnt p;
                gp_Vec v1, v2;
                c.native->D2(tc, p, v1, v2);
                ders_out[0] = {{p.X(), p.Y(), p.Z()}};
                ders_out[1] = {{v1.X(), v1.Y(), v1.Z()}};
                ders_out[2] = {{v2.X(), v2.Y(), v2.Z()}};
                for (int d = 3; d <= max_d; ++d) {
                    gp_Vec vd = c.native->DN(tc, d);
                    ders_out[static_cast<std::size_t>(d)] = {{vd.X(), vd.Y(), vd.Z()}};
                }
            }
        }

        // Shared adapter-side projection body (operator-sinking review,
        // Task 1.3): a t_hint inside the domain is a seeded local search,
        // anything else (incl. NaN) selects the global screening mode.
        static FT project_point_to_curve(const CurvT& c, const Pt3T& p, FT t_hint, FT tol) {
            const std::vector<FT> kv = NurbsEngine::knots(c);
            const FT lo = kv[static_cast<std::size_t>(NurbsEngine::degree(c))];
            const FT hi = kv.back();
            const bool has_hint = (t_hint >= lo && t_hint <= hi);
            const double pt[3] = {p[0], p[1], p[2]};
            const FT tol_use = (tol > FT(0)) ? tol : FT(1e-12);
            return static_cast<FT>(loft_adapter::project_point_to_curve_arrays(
                NurbsEngine::degree(c), kv, NurbsEngine::control_points_homo(c), pt, has_hint,
                has_hint ? static_cast<double>(t_hint) : 0.0, static_cast<double>(tol_use)));
        }

        static bool intersect_curve_plane(const CurvT& c, const Pt3T& plane_q,
                                          const Vec3T& plane_n, FT t_hint,
                                          FT tol, FT& t_out) {
            const FT lo = c.native->FirstParameter();
            const FT hi = c.native->LastParameter();
            auto f = [&](FT t) {
                std::array<Vec3T, 1> d{};
                eval_curve_ders(c, t, 0, Span<Vec3T>(d.data(), 1));
                return Kernel::dot3(Kernel::diff3(d[0], plane_q), plane_n);
            };

            const int steps = 512;
            const FT h = (hi - lo) / static_cast<FT>(steps);
            FT best_t = lo;
            bool found = false;
            FT best_dist = FT(0);
            FT fa = f(lo);
            if (fa == FT(0)) { t_out = lo; return true; }
            FT ta = lo;

            for (int i = 0; i < steps; ++i) {
                const FT tb = (i + 1 == steps) ? hi : lo + static_cast<FT>(i + 1) * h;
                const FT fb = f(tb);
                if (fa * fb < FT(0)) {
                    FT a = ta, b = tb;
                    const FT fa_sign = fa;
                    for (int it = 0; it < 80 && (b - a) > tol; ++it) {
                        const FT m = FT(0.5) * (a + b);
                        const FT fm = f(m);
                        if (fm == FT(0)) { a = m; b = m; break; }
                        if (fa_sign * fm < FT(0)) b = m; else a = m;
                    }
                    const FT cand = FT(0.5) * (a + b);
                    const FT dist = std::fabs(cand - detail::clampd(t_hint, lo, hi));
                    if (!found || dist < best_dist) {
                        found = true;
                        best_t = cand;
                        best_dist = dist;
                    }
                }
                if (fb == FT(0)) {
                    const FT dist = std::fabs(tb - detail::clampd(t_hint, lo, hi));
                    if (!found || dist < best_dist) {
                        found = true;
                        best_t = tb;
                        best_dist = dist;
                    }
                }
                fa = fb;
                ta = tb;
            }
            if (found) t_out = best_t;
            return found;
        }

        static bool extract_supporting_plane(const CurvT& c, FT tol,
                                             Pt3T& q_out, Vec3T& n_out,
                                             FT& max_planarity_err_out) {
            const int ns = 64;
            const FT lo = c.native->FirstParameter();
            const FT hi = c.native->LastParameter();
            std::vector<Pt3T> pts(static_cast<std::size_t>(ns));
            std::array<Vec3T, 1> d{};
            for (int i = 0; i < ns; ++i) {
                const FT t = lo + (hi - lo) * static_cast<FT>(i) / static_cast<FT>(ns - 1);
                eval_curve_ders(c, t, 0, Span<Vec3T>(d.data(), 1));
                pts[static_cast<std::size_t>(i)] = d[0];
            }
            Pt3T q{{0, 0, 0}};
            for (const auto& pt : pts)
                for (int w = 0; w < 3; ++w)
                    q[static_cast<std::size_t>(w)] += pt[static_cast<std::size_t>(w)];
            for (int w = 0; w < 3; ++w)
                q[static_cast<std::size_t>(w)] /= static_cast<FT>(ns);

            Vec3T n{{0, 0, 0}};
            for (std::size_t i = 0; i < pts.size(); ++i) {
                const Pt3T& a = pts[i];
                const Pt3T& b = pts[(i + 1) % pts.size()];
                const Vec3T e1 = Kernel::diff3(b, a);
                const Vec3T e2{{
                    a[0] + b[0] - FT(2) * q[0],
                    a[1] + b[1] - FT(2) * q[1],
                    a[2] + b[2] - FT(2) * q[2] }};
                const Vec3T cr = Kernel::cross3(e1, e2);
                for (int w = 0; w < 3; ++w)
                    n[static_cast<std::size_t>(w)] += cr[static_cast<std::size_t>(w)];
            }
            const FT nlen = Kernel::norm3(n);
            if (nlen <= FT(1e-300)) {
                q_out = q;
                n_out = {{0, 0, 0}};
                max_planarity_err_out = FT(0);
                return false;
            }
            for (int w = 0; w < 3; ++w)
                n[static_cast<std::size_t>(w)] /= nlen;

            FT err = FT(0);
            for (const auto& pt : pts) {
                const FT dist = std::fabs(Kernel::dot3(Kernel::diff3(pt, q), n));
                if (dist > err) err = dist;
            }
            q_out = q;
            n_out = n;
            max_planarity_err_out = err;
            return err <= tol;
        }

        static FT compute_arc_length(const CurvT& c, FT t0, FT t1) {
            static const FT xg[2] = {FT(0.3399810435848563),
                                     FT(0.8611363115940526)};
            static const FT wg[2] = {FT(0.6521451548625461),
                                     FT(0.3478548451374538)};
            const FT lo = c.native->FirstParameter();
            const FT hi = c.native->LastParameter();
            const FT a = detail::clampd(std::fmin(t0, t1), lo, hi);
            const FT b = detail::clampd(std::fmax(t0, t1), lo, hi);
            const int panels = 64;
            const FT hp = (b - a) / static_cast<FT>(panels);
            FT total = FT(0);
            std::array<Vec3T, 2> d{};
            for (int i = 0; i < panels; ++i) {
                const FT am = a + static_cast<FT>(i) * hp;
                for (int g = 0; g < 2; ++g) {
                    for (int sgn = 0; sgn < 2; ++sgn) {
                        const FT xi = (sgn == 0) ? xg[g] : -xg[g];
                        const FT tm = am + FT(0.5) * hp * (FT(1) + xi);
                        eval_curve_ders(c, tm, 1, Span<Vec3T>(d.data(), 2));
                        total += wg[g] * Kernel::norm3(d[1]);
                    }
                }
            }
            total *= FT(0.5) * hp;
            return (t1 >= t0) ? total : -total;
        }

        // Neutral adaptive chord-height polyline (Pillar 3): deterministic
        // midpoint-stack chord-height subdivision with depth cap 8 and a
        // 256-sample ceiling — the same refinement contract as the OD
        // host's VxCrvAproxPolyPar oracle, evaluated over this kernel's
        // own curve evaluation so every four-pillar host shares one
        // longitudinal-enrichment semantics.
        static void approximate_curve_polyline(const CurvT& c, FT tol,
                                               std::vector<loft::CurvePolylineSample<FT, Pt3T>>& out) {
            out.clear();
            const int p = c.degree;
            const FT lo = c.native->FirstParameter();
            const FT hi = c.native->LastParameter();
            if (!(hi > lo)) {
                std::array<Vec3T, 1> d0{};
                eval_curve_ders(c, lo, 0, Span<Vec3T>(d0.data(), d0.size()));
                out.push_back({d0[0], lo});
                return;
            }
            const FT use_tol = (tol > FT(0)) ? tol : FT(1e-4);
            auto eval_pt = [&](FT t, Pt3T& pnt) {
                std::array<Vec3T, 1> d{};
                eval_curve_ders(c, t, 0, Span<Vec3T>(d.data(), d.size()));
                pnt = d[0];
            };
            if (p == 1) {
                const std::size_t num_cp =
                    c.knots.size() - static_cast<std::size_t>(p) - 1;
                out.reserve(num_cp);
                for (std::size_t i = 0; i < num_cp; ++i) {
                    const FT t = c.knots[i + 1];
                    Pt3T pnt{};
                    eval_pt(t, pnt);
                    out.push_back({pnt, t});
                }
                return;
            }
            auto dist_to_seg = [](const Pt3T& q, const Pt3T& a, const Pt3T& b) -> FT {
                const Vec3T ab{{b[0] - a[0], b[1] - a[1], b[2] - a[2]}};
                const FT L2 = Kernel::dot3(ab, ab);
                if (L2 <= FT(1e-14)) {
                    return Kernel::norm3(Kernel::diff3(q, a));
                }
                const FT s = (std::max)(FT(0), (std::min)(FT(1),
                    Kernel::dot3(Kernel::diff3(q, a), ab) / L2));
                const Vec3T pr{{a[0] + s * ab[0], a[1] + s * ab[1], a[2] + s * ab[2]}};
                return Kernel::norm3(Kernel::diff3(q, pr));
            };
            struct Seg { FT ta, tb; Pt3T pa, pb; int depth; };
            std::vector<Seg> stack;
            Pt3T p_start{}, p_end{};
            eval_pt(lo, p_start);
            eval_pt(hi, p_end);
            out.push_back({p_start, lo});
            stack.push_back({lo, hi, p_start, p_end, 0});
            while (!stack.empty() && out.size() < 256) {
                const Seg s = stack.back();
                stack.pop_back();
                const FT tmid = FT(0.5) * (s.ta + s.tb);
                Pt3T pmid{};
                eval_pt(tmid, pmid);
                const FT chord_ht = dist_to_seg(pmid, s.pa, s.pb);
                if (chord_ht > use_tol && s.depth < 8) {
                    stack.push_back({tmid, s.tb, pmid, s.pb, s.depth + 1});
                    stack.push_back({s.ta, tmid, s.pa, pmid, s.depth + 1});
                } else {
                    out.push_back({s.pb, s.tb});
                }
            }
        }

        // Shared adapter-side chain assembly (operator-sinking review,
        // Task 2.3): the former LoftPrepIntent::merge_profile_chain body.
        static CurvT concatenate_curves(Span<const CurvT> segments, Span<const FT> slot_bounds,
                                        std::vector<FT> &joints_out) {
            return loft_adapter::concatenate_curves_via<NurbsEngine, GeomEngine>(segments, slot_bounds,
                                                                                 joints_out);
        }

        // Shared adapter-side extremum oracle (operator-sinking review,
        // Step 1): the former closest_guide_profile_u screening grid.
        static bool compute_curve_curve_extremum(const CurvT& c1, const CurvT& c2, FT& t1_out, FT& t2_out,
                                                 FT& min_dist_out) {
            const std::vector<FT> kv1 = NurbsEngine::knots(c1);
            const std::vector<FT> kv2 = NurbsEngine::knots(c2);
            std::vector<std::array<double, 4>> cp1(NurbsEngine::control_points_homo(c1));
            std::vector<std::array<double, 4>> cp2(NurbsEngine::control_points_homo(c2));
            double t1 = 0.0, t2 = 0.0, dmin = 0.0;
            if (!loft_adapter::curve_curve_min_dist_arrays(NurbsEngine::degree(c1), kv1, cp1,
                                                           NurbsEngine::degree(c2), kv2, cp2, t1, t2, dmin))
                return false;
            t1_out = static_cast<FT>(t1);
            t2_out = static_cast<FT>(t2);
            min_dist_out = static_cast<FT>(dmin);
            return true;
        }
    };

    // -----------------------------------------------------------------------
    // Pillar 4: AlgebraEngine — shared Eigen body.
    // -----------------------------------------------------------------------
    struct AlgebraEngine : loft_adapter::EigenAlgebraEngineBase {
        using FT = double;
    };
};

} // namespace occt_cax

// Contract gate: assert OcctCaxLoftContext satisfies the Four-Pillar contract
static_assert(loft::is_lofting_context<occt_cax::OcctCaxLoftContext>::value,
              "OcctCaxLoftContext must satisfy the Four-Pillar LoftingContext contract");

// Unqualified alias
using OcctCaxLoftContext = occt_cax::OcctCaxLoftContext;

using ::loft::Span;
using ::loft::MatrixView2D;
