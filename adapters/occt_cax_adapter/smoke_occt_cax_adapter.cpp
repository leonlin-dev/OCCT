// ===========================================================================
// adapters/occt_cax_adapter/smoke_occt_cax_adapter.cpp
//
// Smoke witness for the OCCT CAX adapter. Checks:
//   [C1] OcctCaxLoftContext satisfies the Module-1 Four-Pillar contract
//        (compile-time; loft::assert_valid_lofting_context).
//   [C2] NurbsEngine dialect operators behave: span/basis/Oslo row-sum,
//        reparameterize_linear forward + reversal, elevate_degree
//        geometry preservation (rational quarter-circle stays on the arc).
//   [C2b] Oslo refinement identity B_old = M * B_new on an interior-knot vector.
//   [C3] GeomEngine: plane intersection, supporting-plane extraction, arc length.
//   [C4] AlgebraEngine (Eigen): KKT solve of a small saddle system and a
//        DFT cyclic cross-correlation peak.
//   [C4b] Cyclic cross-correlation against the direct O(N^2) reference.
//   [C5] End-to-end occt_cax::loft_row_coplanar over native OCCT handles
//        (4 tilted planar sections + spine) -> native Geom_BSplineSurface
//        with passing audit invariants.
// ===========================================================================
#include "occt_cax_loft_api.hpp"

#include <adapter_test_guard.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

int g_checks = 0;
int g_failures = 0;

void check(bool cond, const char* what) {
    ++g_checks;
    if (!cond) {
        ++g_failures;
        std::printf("  FAIL: %s\n", what);
    }
}

using Ctx    = occt_cax::OcctCaxLoftContext;
using Nurbs  = Ctx::NurbsEngine;
using Geom   = Ctx::GeomEngine;
using Curve  = Nurbs::CurvT;

Curve make_section(double z, double tilt_xy, double scale) {
    TColgp_Array1OfPnt poles(1, 3);
    poles.SetValue(1, gp_Pnt(-scale + tilt_xy, -scale, z));
    poles.SetValue(2, gp_Pnt(scale + tilt_xy, -scale, z));
    poles.SetValue(3, gp_Pnt(tilt_xy, scale, z));
    TColStd_Array1OfReal knots(1, 2);
    knots.SetValue(1, 0.0);
    knots.SetValue(2, 1.0);
    TColStd_Array1OfInteger mults(1, 2);
    mults.SetValue(1, 3);
    mults.SetValue(2, 3);
    Handle(Geom_BSplineCurve) c = new Geom_BSplineCurve(poles, knots, mults, 2);
    return Curve(c);
}

Curve make_spine() {
    TColgp_Array1OfPnt poles(1, 2);
    poles.SetValue(1, gp_Pnt(0.0, 0.0, 0.0));
    poles.SetValue(2, gp_Pnt(0.0, 0.0, 3.0));
    TColStd_Array1OfReal knots(1, 2);
    knots.SetValue(1, 0.0);
    knots.SetValue(2, 1.0);
    TColStd_Array1OfInteger mults(1, 2);
    mults.SetValue(1, 2);
    mults.SetValue(2, 2);
    Handle(Geom_BSplineCurve) g = new Geom_BSplineCurve(poles, knots, mults, 1);
    return Curve(g);
}

Curve make_rational_quarter_circle(double r = 1.0) {
    // 90-degree circular arc: degree 2, 3 CPs, middle weight w1 = 1/sqrt(2).
    const double w1 = 1.0 / std::sqrt(2.0);
    TColgp_Array1OfPnt poles(1, 3);
    poles.SetValue(1, gp_Pnt(r, 0.0, 0.0));
    poles.SetValue(2, gp_Pnt(r, r, 0.0));
    poles.SetValue(3, gp_Pnt(0.0, r, 0.0));
    TColStd_Array1OfReal weights(1, 3);
    weights.SetValue(1, 1.0);
    weights.SetValue(2, w1);
    weights.SetValue(3, 1.0);
    TColStd_Array1OfReal knots(1, 2);
    knots.SetValue(1, 0.0);
    knots.SetValue(2, 1.0);
    TColStd_Array1OfInteger mults(1, 2);
    mults.SetValue(1, 3);
    mults.SetValue(2, 3);
    Handle(Geom_BSplineCurve) c = new Geom_BSplineCurve(poles, weights, knots, mults, 2);
    return Curve(c);
}

void run_engine_checks() {
    // [C2] span + basis partition of unity on the section knot vector.
    const Curve sec = make_section(0.0, 0.0, 1.0);
    const std::size_t span = Nurbs::find_knot_span(2, sec.knots, 0.5);
    check(span == 2, "find_knot_span mid-domain span == p");
    std::vector<double> buf(2 * 3, 0.0);
    loft::MatrixView2D<double> ders(buf.data(), 2, 3);
    Nurbs::eval_basis_ders(2, sec.knots, span, 0.5, 1, ders);
    double sum0 = 0.0;
    for (std::size_t j = 0; j < 3; ++j) sum0 += ders(0, j);
    check(std::fabs(sum0 - 1.0) <= 1e-14, "basis row 0 partition of unity");
    check(std::fabs(ders(1, 0) + ders(1, 2)) <= 1e-12,
          "basis row 1 symmetry (symmetric hat)");

    // [C2] Oslo matrix: column sum == 1
    std::vector<double> kv_old = {0, 0, 0, 0, 1, 1, 1, 1};
    std::vector<double> kv_new = {0, 0, 0, 0, 0.5, 1, 1, 1, 1};
    std::vector<double> m_buf(4 * 5, 0.0);
    loft::MatrixView2D<double> M(m_buf.data(), 4, 5);
    Nurbs::compute_oslo_matrix(3, kv_old, kv_new, M);
    double worst_col = 0.0;
    for (std::size_t b = 0; b < 5; ++b) {
        double col_sum = 0.0;
        for (std::size_t a = 0; a < 4; ++a) col_sum += M(a, b);
        const double err = std::fabs(col_sum - 1.0);
        if (err > worst_col) worst_col = err;
    }
    check(worst_col <= 1e-13, "Oslo matrix column-sum == 1 partition of unity");

    // [C2] reparameterize_linear forward [0,1] -> [2, 5]
    const Curve fwd = Nurbs::reparameterize_linear(sec, 2.0, 5.0);
    check(std::fabs(fwd.knots.front() - 2.0) <= 1e-14 &&
          std::fabs(fwd.knots.back() - 5.0) <= 1e-14,
          "reparameterize_linear forward endpoints");
    std::array<std::array<double, 3>, 1> d_sec{}, d_fwd{};
    Geom::eval_curve_ders(sec, 0.5, 0, loft::Span<std::array<double, 3>>(d_sec.data(), 1));
    Geom::eval_curve_ders(fwd, 3.5, 0, loft::Span<std::array<double, 3>>(d_fwd.data(), 1));
    check(std::fabs(d_sec[0][0] - d_fwd[0][0]) <= 1e-13 &&
          std::fabs(d_sec[0][1] - d_fwd[0][1]) <= 1e-13 &&
          std::fabs(d_sec[0][2] - d_fwd[0][2]) <= 1e-13,
          "reparameterize_linear forward point preservation");

    // [C2] reparameterize_linear reversal [0, 1] -> [1, 0]
    const Curve rev = Nurbs::reparameterize_linear(sec, 1.0, 0.0);
    check(std::fabs(rev.knots.front() - 0.0) <= 1e-14 &&
          std::fabs(rev.knots.back() - 1.0) <= 1e-14,
          "reparameterize_linear reversal knot endpoints normalized");
    std::array<std::array<double, 3>, 1> d_rev{};
    Geom::eval_curve_ders(rev, 0.25, 0, loft::Span<std::array<double, 3>>(d_rev.data(), 1));
    Geom::eval_curve_ders(sec, 0.75, 0, loft::Span<std::array<double, 3>>(d_sec.data(), 1));
    check(std::fabs(d_rev[0][0] - d_sec[0][0]) <= 1e-13 &&
          std::fabs(d_rev[0][1] - d_sec[0][1]) <= 1e-13 &&
          std::fabs(d_rev[0][2] - d_sec[0][2]) <= 1e-13,
          "reparameterize_linear reversal evaluates backwards");

    // [C2] elevate_degree: quarter-circle p=2 -> p=4 stays on arc
    const Curve circ2 = make_rational_quarter_circle(2.0);
    const Curve circ4 = Nurbs::elevate_degree(circ2, 4);
    check(circ4.degree == 4, "elevate_degree target degree == 4");
    double worst_radius_err = 0.0;
    for (int i = 0; i <= 20; ++i) {
        const double t = static_cast<double>(i) / 20.0;
        std::array<std::array<double, 3>, 1> d{};
        Geom::eval_curve_ders(circ4, t, 0, loft::Span<std::array<double, 3>>(d.data(), 1));
        const double r = std::sqrt(d[0][0] * d[0][0] + d[0][1] * d[0][1]);
        const double err = std::fabs(r - 2.0);
        if (err > worst_radius_err) worst_radius_err = err;
    }
    check(worst_radius_err <= 1e-13, "elevate_degree preserves circle radius");

    // [C3] GeomEngine: plane intersection on non-unit-domain curve
    const Curve nonunit = Nurbs::reparameterize_linear(sec, 10.0, 20.0);
    double t_hit = 0.0;
    const std::array<double, 3> plane_q{{0.0, 0.0, 0.0}};
    const std::array<double, 3> plane_n{{1.0, 0.0, 0.0}}; // x = 0
    const bool hit = Geom::intersect_curve_plane(nonunit, plane_q, plane_n, 15.0, 1e-12, t_hit);
    check(hit, "intersect_curve_plane finds root on non-unit domain");
    std::array<std::array<double, 3>, 1> d_hit{};
    Geom::eval_curve_ders(nonunit, t_hit, 0, loft::Span<std::array<double, 3>>(d_hit.data(), 1));
    check(std::fabs(d_hit[0][0]) <= 1e-11, "intersect_curve_plane root x coord == 0");

    // [C3] supporting plane extraction on planar hat
    std::array<double, 3> sup_q{}, sup_n{};
    double max_err = 0.0;
    const bool pl_ok = Geom::extract_supporting_plane(sec, 1e-8, sup_q, sup_n, max_err);
    check(pl_ok && max_err <= 1e-12, "extract_supporting_plane planar section");
    check(std::fabs(std::fabs(sup_n[2]) - 1.0) <= 1e-12,
          "extract_supporting_plane normal aligned with z");

    // [C3] arc length on spine gamma(v) = (0, 0, 3v), [0, 1] -> 3.0
    const Curve spine = make_spine();
    const double len = Geom::compute_arc_length(spine, 0.0, 1.0);
    check(std::fabs(len - 3.0) <= 1e-12, "compute_arc_length straight spine len == 3.0");

    // [C4] AlgebraEngine: KKT solve
    std::vector<double> Q = {2.0};
    std::vector<double> A = {1.0};
    std::vector<double> rq = {2.0};
    std::vector<double> ra = {1.0};
    std::vector<double> x(1, 0.0);
    const bool kkt = Ctx::AlgebraEngine::solve_kkt(
        loft::MatrixView2D<const double>(Q.data(), 1, 1),
        loft::MatrixView2D<const double>(A.data(), 1, 1),
        rq, ra, x);
    check(kkt && std::fabs(x[0] - 1.0) <= 1e-12, "KKT saddle solve x = 1");

    // [C4] DFT cyclic cross-correlation
    std::vector<double> sig(8, 0.0), corr(8, 0.0);
    sig[2] = 1.0;
    Ctx::AlgebraEngine::circular_cross_correlation(sig, sig, corr);
    double peak = corr[0];
    bool peak_at_zero = true;
    for (std::size_t tau = 1; tau < 8; ++tau)
        if (corr[tau] > peak + 1e-12) peak_at_zero = false;
    check(peak_at_zero && std::fabs(peak - 1.0) <= 1e-12,
          "DFT cyclic cross-correlation self-peak at tau=0");
}

void run_oslo_identity_check() {
    const int p = 3;
    const std::vector<double> kv_old = {0, 0, 0, 0, 0.35, 0.7, 1, 1, 1, 1};
    const std::vector<double> kv_new = {0, 0, 0, 0, 0.35, 0.5, 0.7, 1, 1, 1, 1};
    const std::size_t n_old = kv_old.size() - 4;  // 6
    const std::size_t n_new = kv_new.size() - 4;  // 7
    std::vector<double> m_buf(n_old * n_new, 0.0);
    Nurbs::compute_oslo_matrix(
        p, kv_old, kv_new,
        loft::MatrixView2D<double>(m_buf.data(), n_old, n_new));

    const std::size_t pp_sz = static_cast<std::size_t>(p) + 1;
    auto basis_row = [&](const std::vector<double>& kv, double u) {
        std::vector<double> row(kv.size() - pp_sz, 0.0);
        const std::size_t span = Nurbs::find_knot_span(p, kv, u);
        std::vector<double> buf(pp_sz, 0.0);
        loft::MatrixView2D<double> ders(buf.data(), 1, pp_sz);
        Nurbs::eval_basis_ders(p, kv, span, u, 0, ders);
        for (std::size_t j = 0; j < pp_sz; ++j)
            row[span - static_cast<std::size_t>(p) + j] = buf[j];
        return row;
    };

    const std::vector<double> us = {0.0, 0.1, 0.27, 0.35, 0.5, 0.66, 0.7, 0.83, 1.0};
    double worst = 0.0;
    for (const double u : us) {
        const std::vector<double> n_old_row = basis_row(kv_old, u);
        const std::vector<double> n_new_row = basis_row(kv_new, u);
        for (std::size_t a = 0; a < n_old; ++a) {
            double recon = 0.0;
            for (std::size_t b = 0; b < n_new; ++b)
                recon += m_buf[a * n_new + b] * n_new_row[b];
            const double res = std::fabs(recon - n_old_row[a]);
            if (res > worst) worst = res;
        }
    }
    check(worst <= 1e-13, "Oslo identity B_old = M*B_new on interior knots");
    std::printf("  oslo identity worst residual=%.3e\n", worst);
}

void run_correlation_reference_check() {
    const std::size_t n = 16;
    std::vector<double> a(n), b(n);
    unsigned int seed = 12345u;
    for (std::size_t i = 0; i < n; ++i) {
        seed = seed * 1103515245u + 12345u;
        a[i] = static_cast<double>(static_cast<int>((seed >> 16) % 2000u) - 1000) / 1000.0;
        seed = seed * 1103515245u + 12345u;
        b[i] = static_cast<double>(static_cast<int>((seed >> 16) % 2000u) - 1000) / 1000.0;
    }
    std::vector<double> corr(n, 0.0);
    Ctx::AlgebraEngine::circular_cross_correlation(a, b, corr);
    double worst = 0.0;
    for (std::size_t tau = 0; tau < n; ++tau) {
        double ref = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            ref += a[i] * b[(i + tau) % n];
        const double res = std::fabs(corr[tau] - ref);
        if (res > worst) worst = res;
    }
    check(worst <= 1e-9, "cyclic correlation matches direct reference");
    std::printf("  correlation vs direct reference worst=%.3e\n", worst);
}

void run_loft_end_to_end() {
    // [C5] 4 tilted planar sections along the spine -> native OCCT Geom_BSplineSurface
    std::vector<Curve> sections;
    for (int i = 0; i < 4; ++i)
        sections.push_back(make_section(static_cast<double>(i), 0.1 * i, 1.0));
    const Curve spine = make_spine();

    occt_cax::LoftRequest req;
    req.degree_u = 2;
    req.degree_v = 3;

    const occt_cax::LoftRowResult r =
        occt_cax::loft_row_coplanar(sections, {}, &spine, req);

    check(!r.surface.IsNull(), "lofted native Geom_BSplineSurface is non-null");
    check(r.structural_ok, "restricted structural_ok is true");
    check(r.num_cp_u >= 3 && r.num_cp_v >= 2, "surface control-grid sizes");
    check(r.audit.all_passed, "audit invariants all passed");
    check(r.audit.max_row_coplanarity_error <= 1e-14, "row coplanarity <= 1e-14");
    check(r.audit.max_profile_algebraic_error <= 1e-13, "profile algebraic interpolation <= 1e-13");

    check(r.surface->UDegree() == 2 && r.surface->VDegree() == 3, "native surface degrees preserved");
    check(static_cast<std::size_t>(r.surface->NbUPoles()) == r.num_cp_u, "native NbUPoles matches num_cp_u");
    check(static_cast<std::size_t>(r.surface->NbVPoles()) == r.num_cp_v, "native NbVPoles matches num_cp_v");

    std::printf(
        "  loft: N_u=%zu N_v=%zu NbUPoles=%d NbVPoles=%d "
        "coplanarity=%.3e algebraic=%.3e guide_tol=%.3e isoline_bound=%.3e\n",
        r.num_cp_u, r.num_cp_v, r.surface->NbUPoles(), r.surface->NbVPoles(),
        r.audit.max_row_coplanarity_error, r.audit.max_profile_algebraic_error,
        r.audit.max_guide_tolerance_error, r.audit.max_isoline_planar_bound);
}

} // namespace

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    return adapter_guard::run("smoke_occt_cax_adapter", []() {
        loft::assert_valid_lofting_context<Ctx>();

        run_engine_checks();
        run_oslo_identity_check();
        run_correlation_reference_check();
        run_loft_end_to_end();

        std::printf("smoke_occt_cax_adapter: checks=%d failures=%d\n",
                    g_checks, g_failures);
        return (g_failures == 0) ? 0 : 1;
    });
}
