// LFP half-cell model: Li foil | separator | porous LiFePO4 cathode.
//
// A self-contained C++17 program: edit the input file, run, and it writes Time_Voltage.txt.
// It reads the same namelist input as the Fortran program and the Python package, and
// contains its own copy of the BAND block-tridiagonal solver (adapted from bandsolver
// v0.1.1, https://github.com/jcbernard87/bandsolver, BSD-3-Clause), so it needs no library.
//
//   build:  c++ -std=c++17 -O2 -ffp-contract=off lfp.cpp -o lfp_cpp
//   usage:  lfp_cpp [input.nml]          (default input file: lfp.nml)
//
// Equations: docs/model.md. Parameters: docs/parameters.md. mode = 'faithful' reproduces
// the original research code including the defects in docs/deviations.md; mode =
// 'corrected' applies the fixes.
//
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// ============================== BAND solver ==============================
// Solves A[j] dc[j-1] + B[j] dc[j] + D[j] dc[j+1] = G[j], j = 0..nj-1 (row-major n x n blocks).
namespace band {

enum class Pivot { partial, legacy };

// Solve Bm * S = R in place (Bm n x n, R n x m, row-major). Returns false if singular.
bool solve_partial(int n, int m, double* Bm, double* R) {
    const double eps = std::numeric_limits<double>::epsilon();
    double bscale = 0;
    for (int i = 0; i < n * n; ++i) bscale = std::max(bscale, std::abs(Bm[i]));
    if (bscale == 0) return false;
    for (int k = 0; k < n; ++k) {
        int p = k;
        for (int i = k + 1; i < n; ++i)
            if (std::abs(Bm[i * n + k]) > std::abs(Bm[p * n + k])) p = i;
        if (std::abs(Bm[p * n + k]) <= n * eps * bscale) return false;
        if (p != k) {
            std::swap_ranges(Bm + p * n, Bm + p * n + n, Bm + k * n);
            std::swap_ranges(R + p * m, R + p * m + m, R + k * m);
        }
        double f = 1.0 / Bm[k * n + k];
        for (int c = k; c < n; ++c) Bm[k * n + c] *= f;
        for (int c = 0; c < m; ++c) R[k * m + c] *= f;
        for (int i = 0; i < n; ++i) {
            if (i == k) continue;
            f = Bm[i * n + k];
            if (f == 0) continue;
            for (int c = k; c < n; ++c) Bm[i * n + c] -= f * Bm[k * n + c];
            for (int c = 0; c < m; ++c) R[i * m + c] -= f * R[k * m + c];
        }
    }
    return true;
}

// Pivot heuristic of the archival MATINV (Newman, Appendix C), with the operation order kept.
bool solve_legacy(int n, int m, double* Bm, double* R) {
    const double eps = std::numeric_limits<double>::epsilon();
    const double bmax0 = static_cast<double>(1.1f);
    double bscale = 0;
    for (int i = 0; i < n * n; ++i) bscale = std::max(bscale, std::abs(Bm[i]));
    std::vector<char> used(n, 0);
    int irow = 0, jcol = 0, jc = 0;
    for (int nn = 0; nn < n; ++nn) {
        double bmax = bmax0;
        bool found = false;
        for (int i = 0; i < n; ++i) {
            if (used[i]) continue;
            double bnext = 0, btry = 0;
            for (int j = 0; j < n; ++j) {
                if (used[j]) continue;
                const double a = std::abs(Bm[i * n + j]);
                if (a <= bnext) continue;
                bnext = a;
                if (bnext <= btry) continue;
                bnext = btry;
                btry = a;
                jc = j;
                found = true;
            }
            if (bnext >= bmax * btry) continue;
            bmax = bnext / btry;
            irow = i;
            jcol = jc;
        }
        if (!found) return false;
        used[jcol] = 1;
        if (jcol != irow) {
            std::swap_ranges(Bm + irow * n, Bm + irow * n + n, Bm + jcol * n);
            std::swap_ranges(R + irow * m, R + irow * m + m, R + jcol * m);
        }
        if (std::abs(Bm[jcol * n + jcol]) <= n * eps * bscale) return false;
        double f = 1.0 / Bm[jcol * n + jcol];
        for (int j = 0; j < n; ++j) Bm[jcol * n + j] *= f;
        for (int k = 0; k < m; ++k) R[jcol * m + k] *= f;
        for (int i = 0; i < n; ++i) {
            if (i == jcol) continue;
            f = Bm[i * n + jcol];
            for (int j = 0; j < n; ++j) Bm[i * n + j] -= f * Bm[jcol * n + j];
            for (int k = 0; k < m; ++k) R[i * m + k] -= f * R[jcol * m + k];
        }
    }
    return true;
}

// Block elimination and back substitution. Returns false on a singular or non-finite system.
bool solve(int n, int nj, const std::vector<double>& A, const std::vector<double>& B,
           const std::vector<double>& D, const std::vector<double>& G, std::vector<double>& dc, Pivot pivot) {
    const std::size_t nn = static_cast<std::size_t>(n) * n;
    dc.assign(static_cast<std::size_t>(n) * nj, 0.0);
    for (const auto* v : {&A, &B, &D, &G})
        for (double x : *v)
            if (!std::isfinite(x)) return false;
    const int np1 = n + 1;
    std::vector<double> E(static_cast<std::size_t>(nj) * n * np1), Am(nn), Bm(nn), R(static_cast<std::size_t>(n) * np1);
    auto Ej = [&](int j) { return E.data() + static_cast<std::size_t>(j) * n * np1; };
    auto blk = [&](const std::vector<double>& M, int j) { return M.data() + static_cast<std::size_t>(j) * nn; };
    auto block_solve = [&](int m) {
        return pivot == Pivot::legacy ? solve_legacy(n, m, Bm.data(), R.data()) : solve_partial(n, m, Bm.data(), R.data());
    };
    // node 0: B0 [E0 | e0] = [D0 | G0]
    std::copy(blk(B, 0), blk(B, 0) + nn, Bm.begin());
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) R[i * np1 + k] = blk(D, 0)[i * n + k];
        R[i * np1 + n] = G[i];
    }
    if (!block_solve(np1)) return false;
    for (int k = 0; k < n; ++k) {
        Ej(0)[k * np1 + n] = R[k * np1 + n];
        for (int l = 0; l < n; ++l) Ej(0)[k * np1 + l] = -R[k * np1 + l];
    }
    for (int j = 1; j < nj; ++j) {
        std::copy(blk(A, j), blk(A, j) + nn, Am.begin());
        std::copy(blk(B, j), blk(B, j) + nn, Bm.begin());
        for (int i = 0; i < n; ++i)
            for (int k = 0; k < n; ++k) R[i * np1 + k] = blk(D, j)[i * n + k];
        const double* E1 = Ej(j - 1);
        for (int i = 0; i < n; ++i) {
            R[i * np1 + n] = -G[static_cast<std::size_t>(j) * n + i];
            for (int l = 0; l < n; ++l) {
                R[i * np1 + n] += Am[i * n + l] * E1[l * np1 + n];
                for (int k = 0; k < n; ++k) Bm[i * n + k] += Am[i * n + l] * E1[l * np1 + k];
            }
        }
        if (!block_solve(np1)) return false;
        for (int k = 0; k < n; ++k)
            for (int q = 0; q < np1; ++q) Ej(j)[k * np1 + q] = -R[k * np1 + q];
    }
    for (int k = 0; k < n; ++k) dc[static_cast<std::size_t>(nj - 1) * n + k] = Ej(nj - 1)[k * np1 + n];
    for (int j = nj - 2; j >= 0; --j) {
        const double* Ec = Ej(j);
        double* x = dc.data() + static_cast<std::size_t>(j) * n;
        const double* xn = x + n;
        for (int k = 0; k < n; ++k) {
            x[k] = Ec[k * np1 + n];
            for (int l = 0; l < n; ++l) x[k] += Ec[k * np1 + l] * xn[l];
        }
    }
    for (double x : dc)
        if (!std::isfinite(x)) return false;
    return true;
}

}  // namespace band

// ============================== input ==============================
double r32(double x) { return static_cast<double>(static_cast<float>(x)); }

// ---- corrected mode: log variables and Scharfetter-Gummel fluxes (docs/model.md section 10) ----
double bern(double x) {  // B(x) = x/(e^x - 1), with its series near 0
    if (std::abs(x) < 1.0e-3) return 1.0 - x / 2.0 + x * x / 12.0;
    if (x > 700.0) return x * std::exp(-x);
    return x / (std::exp(x) - 1.0);
}
double bern_p(double x) {  // B'(x) = B(x) (1 - B(-x)) / x, with its series near 0
    if (std::abs(x) < 1.0e-3) return -0.5 + x / 6.0 - x * x * x / 180.0;
    return bern(x) * (1.0 - bern(-x)) / x;
}
double softplus(double x) { return std::max(x, 0.0) + std::log(1.0 + std::exp(-std::abs(x))); }
// 1/(1 + e^-x) without cancellation for either sign (use sigm(-x) for 1 - sigm(x))
double sigm(double x) {
    const double e = std::exp(-std::abs(x));
    return x >= 0.0 ? 1.0 / (1.0 + e) : e / (1.0 + e);
}
// x**n by repeated multiplication (the same rounding as the Fortran program)
double ipow(double x, int n) {
    double r = 1.0;
    for (int q = 0; q < n; ++q) r = r * x;
    return r;
}

struct Params {
    double L_cath_um = 24.0, L_sep = 25.0e-4;  // cathode thickness in um (L_cath = L_cath_um*1e-4 cm)
    int nj = 101, sep_node = 22;
    int nj_crystal = 21;                   // crystal model: nodes across a crystal (center and surface included)
    std::string particle_model = "uniform";  // corrected mode: 'uniform' particles or 'crystal' (solid diffusion)
    std::string crystal_shape = "sphere";    // crystal model: 'sphere', 'cylinder' or 'slab'
    double D_c = 8.0e-14;                  // crystal model: solid diffusivity [cm2/s]
    double eps = 0.5, eps_AM = 0.8, eps_sep = 0.39, tau_sep = 4.0, bruggeman = -0.5;
    double f_AM = 0.8;  // corrected mode: active fraction of the solid phase (D-9)
    double D = 2.0e-6, t_plus = 0.25, c_bulk = 1.0e-3, z_plus = 1.0, z_minus = -1.0;
    double sigma = 3.0e-3, M = 125.759, rho = 3.6, Q_th = 0.170, R_p = 200.0e-7;
    double k_rxn = -1.0;  // < 0: default for the chosen mode
    double alpha_a = 0.5, alpha_c = 0.5, k_Li = 1.0e-6, c_Li_ref = 1.0e-3;
    double R = 8.314, T = 298.0, F = 96485.0;
    double C_rate = 1.0, phi1_init = 3.6, phi2_init = 0.0, cs_init = 1.0e-5, t_max = 36000.0;
    int n_steps = 36000;
    double V_min = 2.5, V_max = 4.2;       // corrected-mode cutoffs (D-3)
    double fd_step = 1.0e-6;
    double newton_tol = 1.0e-10;           // corrected-mode Newton tolerance (D-7)
    int newton_max_iter = 25;
    std::string mode = "faithful";
    std::string steps;                     // corrected-mode protocol (docs/protocol.md)
    int cycles = 1;
    double write_interval = 18.0;          // [s], corrected mode
    double kappa_bg = 1.0e-8;              // corrected: background (solvent) ionic conductivity [S/cm]
    std::string file = "Time_Voltage.txt";
};

std::string lower(std::string s) {
    for (auto& ch : s) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    return s;
}

// Minimal Fortran-namelist reader: &group name = value, ... / ; '!' comments; d exponents.
Params read_input(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("input file not found: " + path);
    std::string text, line;
    while (std::getline(in, line)) {
        std::string out;
        char quote = 0;
        for (char ch : line) {
            if (quote) { if (ch == quote) quote = 0; }
            else if (ch == '\'' || ch == '"') quote = ch;
            else if (ch == '!') break;
            out += ch;
        }
        text += out + "\n";
    }
    Params p;
    std::map<std::string, double*> reals = {
        {"l_cath_um", &p.L_cath_um}, {"l_sep", &p.L_sep}, {"eps", &p.eps}, {"eps_am", &p.eps_AM}, {"f_am", &p.f_AM},
        {"eps_sep", &p.eps_sep}, {"tau_sep", &p.tau_sep}, {"bruggeman", &p.bruggeman}, {"d", &p.D},
        {"t_plus", &p.t_plus}, {"c_bulk", &p.c_bulk}, {"d_c", &p.D_c}, {"kappa_bg", &p.kappa_bg}, {"z_plus", &p.z_plus}, {"z_minus", &p.z_minus},
        {"sigma", &p.sigma}, {"m", &p.M}, {"rho", &p.rho}, {"q_th", &p.Q_th}, {"r_p", &p.R_p},
        {"k_rxn", &p.k_rxn}, {"alpha_a", &p.alpha_a}, {"alpha_c", &p.alpha_c}, {"k_li", &p.k_Li},
        {"c_li_ref", &p.c_Li_ref}, {"r", &p.R}, {"t", &p.T}, {"f", &p.F}, {"c_rate", &p.C_rate},
        {"phi1_init", &p.phi1_init}, {"phi2_init", &p.phi2_init}, {"cs_init", &p.cs_init},
        {"t_max", &p.t_max}, {"fd_step", &p.fd_step}, {"v_min", &p.V_min}, {"v_max", &p.V_max},
        {"newton_tol", &p.newton_tol}, {"write_interval", &p.write_interval}};
    std::map<std::string, int*> ints = {{"nj", &p.nj}, {"sep_node", &p.sep_node}, {"n_steps", &p.n_steps},
                                        {"newton_max_iter", &p.newton_max_iter}, {"cycles", &p.cycles},
                                        {"nj_crystal", &p.nj_crystal}};
    std::map<std::string, std::string*> strs = {{"mode", &p.mode}, {"file", &p.file}, {"steps", &p.steps},
                                                {"particle_model", &p.particle_model},
                                                {"crystal_shape", &p.crystal_shape}};
    const std::regex group(R"(&(\w+)([\s\S]*?)/)");
    const std::regex entry(R"((\w+)\s*=\s*('[^']*'|"[^"]*"|[^,\s/]+))");
    for (std::sregex_iterator g(text.begin(), text.end(), group), end; g != end; ++g) {
        const std::string body = (*g)[2];
        for (std::sregex_iterator e(body.begin(), body.end(), entry); e != end; ++e) {
            const std::string name = lower((*e)[1]);
            std::string val = (*e)[2];
            if (strs.count(name)) {
                *strs[name] = (val.front() == '\'' || val.front() == '"') ? val.substr(1, val.size() - 2) : val;
            } else {
                std::replace(val.begin(), val.end(), 'd', 'e');
                std::replace(val.begin(), val.end(), 'D', 'e');
                if (ints.count(name)) *ints[name] = std::stoi(val);
                else if (reals.count(name)) *reals[name] = std::stod(val);
                else throw std::runtime_error("unknown name in input: " + name);
            }
        }
    }
    return p;
}

// ============================== model ==============================
constexpr int NV = 4, IC = 0, IP1 = 1, IP2 = 2, ICS = 3;

// ============================== protocol ==============================
enum class Kind { cc, cv, rest };
struct Step {
    Kind kind;
    double C = 0.0, V = 0.0, t = -1.0, Vmin = 0.0, Vmax = 0.0, Imin = -1.0;  // t, Imin < 0: unset
};

std::vector<Step> parse_protocol(const std::string& text, double C_rate, double V_min, double V_max, int cycles) {
    std::vector<Step> steps;
    auto blank = [](const std::string& x) { return x.find_first_not_of(" \t") == std::string::npos; };
    if (blank(text)) {
        Step st{Kind::cc};
        st.C = C_rate; st.Vmin = V_min; st.Vmax = V_max;
        return {st};
    }
    std::stringstream parts(text);
    std::string part;
    int n = 0;
    while (std::getline(parts, part, ';')) {
        if (blank(part)) continue;
        ++n;
        auto fail = [&](const std::string& msg) { throw std::runtime_error("protocol step " + std::to_string(n) + ": " + msg); };
        std::istringstream words(part);
        std::string w;
        words >> w;
        Step st{Kind::cc};
        const std::string kind = lower(w);
        if (kind == "cc") st.kind = Kind::cc;
        else if (kind == "cv") st.kind = Kind::cv;
        else if (kind == "rest") st.kind = Kind::rest;
        else fail("unknown step type " + w);
        st.Vmin = V_min; st.Vmax = V_max;
        bool hasC = false, hasV = false;
        while (words >> w) {
            const auto eq = w.find('=');
            if (eq == std::string::npos) fail("expected key=value, got " + w);
            const std::string key = lower(w.substr(0, eq));
            std::string val = w.substr(eq + 1);
            std::replace(val.begin(), val.end(), 'd', 'e');
            std::replace(val.begin(), val.end(), 'D', 'e');
            double v;
            try { v = std::stod(val); } catch (...) { fail("bad number " + val); }
            if (st.kind == Kind::cc) {
                if (key == "c") { st.C = v; hasC = true; }
                else if (key == "t") st.t = v;
                else if (key == "vmin") st.Vmin = v;
                else if (key == "vmax") st.Vmax = v;
                else fail("cc does not take " + key);
            } else if (st.kind == Kind::cv) {
                if (key == "v") { st.V = v; hasV = true; }
                else if (key == "t") st.t = v;
                else if (key == "imin") st.Imin = v;
                else fail("cv does not take " + key);
            } else {
                if (key != "t") fail("rest does not take " + key);
                st.t = v;
            }
            if (key == "t" && v <= 0) fail("t must be positive");
        }
        if (st.kind == Kind::cc && !hasC) fail("cc needs C=");
        if (st.kind == Kind::cv && !hasV) fail("cv needs V=");
        if (st.kind == Kind::cv && st.t < 0 && st.Imin < 0) fail("cv needs t= or Imin= to end");
        if (st.kind == Kind::rest && st.t < 0) fail("rest needs t=");
        steps.push_back(st);
    }
    std::vector<Step> all;
    for (int k = 0; k < std::max(1, cycles); ++k) all.insert(all.end(), steps.begin(), steps.end());
    return all;
}
using Mat = double[NV][NV];

struct Model {
    Params p;
    bool faithful;
    double L_cath, lit36, ocp_c[7], spec_a, tortuosity, i_spec, eps_sep_face, phi1_sign, mass_area, i_1C;
    double vf_AM;  // active volume fraction: eps_AM (faithful) or f_AM*(1-eps) (corrected)
    double i_app;  // applied current density [A/cm2]; changes step by step in corrected mode
    bool full_current;
    double dcat_s, dan_s, ucat_s, uan_s, dcat_c, dan_c, ucat_c, uan_c;
    double dplus = 0.0, dminus = 0.0;  // ion diffusivities (corrected mode)
    int s;  // 0-based interface node
    std::vector<double> dx, aW, aE, bW, bE;
    // crystal model (corrected mode; docs/model.md section 12): one crystal per cathode volume,
    // electrode nodes s+1 .. nj-2; vertex-centred mesh with node volumes xtal_V and face areas xtal_A
    bool crystal = false;
    int nl = 0;
    double xtal_h = 0.0, xtal_AR = 0.0, a_x = 0.0;
    std::vector<double> xtal_V, xtal_A;

    explicit Model(Params in) : p(std::move(in)) {
        static const double OCP_FIT[7] = {3.114559, 4.438792, 71.7352, 70.85337, 4.240252, 68.5605, 67.730082};
        if (p.mode == "faithful") faithful = true;
        else if (p.mode == "corrected") faithful = false;
        else throw std::runtime_error("mode must be 'faithful' or 'corrected'");
        if (faithful) {  // single-precision literals of the original (D-4)
            p.R = r32(p.R); p.c_bulk = r32(p.c_bulk); p.Q_th = r32(p.Q_th); p.M = r32(p.M); p.rho = r32(p.rho);
            p.phi1_init = r32(p.phi1_init); p.eps_sep = r32(p.eps_sep); p.eps_AM = r32(p.eps_AM);
            p.C_rate = r32(p.C_rate);
            if (p.k_rxn < 0) p.k_rxn = 1.0e-8 * r32(std::pow(10.0, r32(0.966)));
            lit36 = r32(3.6);
            for (int k = 0; k < 7; ++k) ocp_c[k] = r32(OCP_FIT[k]);
            eps_sep_face = p.eps;   // D-2
            phi1_sign = -1.0;       // D-1
            full_current = false;   // D-11
        } else {
            if (p.k_rxn < 0) p.k_rxn = 1.0e-8 * std::pow(10.0, 0.966);
            lit36 = 3.6;
            for (int k = 0; k < 7; ++k) ocp_c[k] = OCP_FIT[k];
            eps_sep_face = p.eps_sep;
            phi1_sign = 1.0;
            full_current = true;
        }
        vf_AM = faithful ? p.eps_AM : p.f_AM * (1.0 - p.eps);
        L_cath = p.L_cath_um * 1.0e-4;  // as in the original, 24 * 1.0d-4
        spec_a = 3 * vf_AM / p.R_p;
        tortuosity = std::pow(p.eps, p.bruggeman);
        i_spec = p.Q_th * p.C_rate;
        i_app = i_spec * L_cath * vf_AM * p.rho;
        mass_area = L_cath * vf_AM * p.rho;
        i_1C = p.Q_th * mass_area;
        const double t_an = 1.0 - p.t_plus;
        const double d_cat = p.D * (1.0 + (t_an / p.t_plus)) / (2.0 * t_an / p.t_plus);
        const double d_an = d_cat * t_an / p.t_plus;
        dplus = d_cat; dminus = d_an;
        const double u_cat = d_cat / (p.R * p.T), u_an = d_an / (p.R * p.T);
        dcat_s = d_cat / p.tau_sep; dan_s = d_an / p.tau_sep; ucat_s = u_cat / p.tau_sep; uan_s = u_an / p.tau_sep;
        dcat_c = d_cat / tortuosity; dan_c = d_an / tortuosity; ucat_c = u_cat / tortuosity; uan_c = u_an / tortuosity;

        const int nj = p.nj;
        s = p.sep_node - 1;
        dx.assign(nj, 0.0); aW.assign(nj, 0.0); aE.assign(nj, 0.0); bW.assign(nj, 0.0); bE.assign(nj, 0.0);
        const double h_sep = p.L_sep / static_cast<double>(p.sep_node - 2);
        const double h_cat = L_cath / static_cast<double>(nj - p.sep_node - 1);
        for (int j = 1; j < s; ++j) dx[j] = h_sep;
        for (int j = s + 1; j < nj - 1; ++j) dx[j] = h_cat;
        for (int j = 1; j < nj; ++j) { aW[j] = dx[j - 1] / (dx[j - 1] + dx[j]); bW[j] = 2.0 / (dx[j - 1] + dx[j]); }
        for (int j = 0; j < nj - 1; ++j) { aE[j] = dx[j] / (dx[j + 1] + dx[j]); bE[j] = 2.0 / (dx[j] + dx[j + 1]); }
        crystal_setup();
    }

    // validate the particle-model inputs and build the crystal mesh (see crystal.py in Python)
    void crystal_setup() {
        if (p.particle_model != "uniform" && p.particle_model != "crystal")
            throw std::runtime_error("particle_model must be one of ('uniform', 'crystal'), got '" + p.particle_model + "'");
        int k;
        if (p.crystal_shape == "slab") k = 0;
        else if (p.crystal_shape == "cylinder") k = 1;
        else if (p.crystal_shape == "sphere") k = 2;
        else throw std::runtime_error("crystal_shape must be one of ('slab', 'cylinder', 'sphere'), got '" + p.crystal_shape + "'");
        crystal = p.particle_model == "crystal";
        nl = 0;
        if (!crystal) return;
        if (faithful) throw std::runtime_error("particle_model='crystal' needs mode='corrected'");
        if (p.nj_crystal < 4) throw std::runtime_error("nj_crystal must be at least 4, got " + std::to_string(p.nj_crystal));
        const int nc = p.nj_crystal;
        nl = p.nj - 2 - s;
        xtal_h = p.R_p / static_cast<double>(nc - 1);
        xtal_V.assign(nc, 0.0);
        xtal_A.assign(nc - 1, 0.0);
        for (int j = 0; j < nc; ++j) {
            double rr = xtal_h * static_cast<double>(j);
            if (j == nc - 1) rr = p.R_p;
            const double rW = std::max(rr - xtal_h / 2.0, 0.0), rE = std::min(rr + xtal_h / 2.0, p.R_p);
            xtal_V[j] = (ipow(rE, k + 1) - ipow(rW, k + 1)) / static_cast<double>(k + 1);
        }
        for (int j = 0; j < nc - 1; ++j) xtal_A[j] = ipow(xtal_h * (static_cast<double>(j) + 0.5), k);
        xtal_AR = ipow(p.R_p, k);
        a_x = static_cast<double>(k + 1) * vf_AM / p.R_p;
    }

    // ---- kinetics ----
    double cs_max() const { return (p.rho / p.M) * p.M * p.Q_th * 1000.0 * lit36 / p.F; }
    // open-circuit potential; corrected mode adds the Nernst term (RT/F) ln(c/c_bulk) (D-14)
    double ocp(double cs, double c) const {
        const double th = (cs / (p.rho / p.M)) / (p.M * p.Q_th * 1000.0 * lit36 / p.F);
        double u = ocp_c[0] + ocp_c[1] * std::atan(-(ocp_c[2] * th) + ocp_c[3]) - ocp_c[4] * std::atan(-(ocp_c[5] * th) + ocp_c[6]);
        if (!faithful) u += p.R * p.T / p.F * std::log(c / p.c_bulk);
        return u;
    }
    double rate(double c, double cs, double p1, double p2) const {
        const double eta = p1 - p2 - ocp(cs, c);
        double i0 = p.F * p.k_rxn * std::pow(c, p.alpha_a) * std::pow(cs_max() - cs, p.alpha_a) * std::pow(cs, p.alpha_c);
        if (faithful) i0 = r32(i0);  // D-4
        return i0 * (std::exp(p.alpha_a * p.F * eta / (p.R * p.T)) - std::exp(-(p.alpha_c * p.F * eta / (p.R * p.T))));
    }
    // rate and finite-difference derivatives w.r.t. (c, phi1, phi2, cs) (D-6)
    // rate and finite-difference derivatives (D-6); faithful mode only
    void rate_derivs(double c, double cs, double p1, double p2, double& i, double di[NV]) const {
        const double h = p.fd_step;
        i = rate(c, cs, p1, p2);
        di[IC] = c <= h ? (rate(c + h, cs, p1, p2) - i) / h
                        : (rate(c + h, cs, p1, p2) - rate(c - h, cs, p1, p2)) / (2.0 * h);
        di[ICS] = cs <= h ? (rate(c, cs + h, p1, p2) - i) / h
                          : (rate(c, cs + h, p1, p2) - rate(c, cs - h, p1, p2)) / (2.0 * h);
        di[IP1] = (rate(c, cs, p1 + h, p2) - rate(c, cs, p1 - h, p2)) / (2.0 * h);
        di[IP2] = (rate(c, cs, p1, p2 + h) - rate(c, cs, p1, p2 - h)) / (2.0 * h);
    }

    // ---- corrected-mode time step ----
    // storage coefficients: the time-derivative part of each row is T * (c - c_old)

    // ---- assembly ----
    void face(double e, double dcat, double ucat, double dan, double uan, double cf, double gphi2, Mat& dd, Mat& ff) const {
        const double k = p.z_plus * p.z_plus * ucat + p.z_minus * p.z_minus * uan;
        dd[IC][IC] = -(e * dcat);
        ff[IC][IC] = -(e * p.z_plus * ucat * p.F * gphi2);
        dd[IC][IP2] = -(e * p.z_plus * ucat * p.F * cf);
        dd[IP2][IC] = -(e * p.F * (p.z_plus * dcat + p.z_minus * dan));
        ff[IP2][IC] = -(e * (p.F * p.F) * k * gphi2);
        dd[IP2][IP2] = -(e * (p.F * p.F) * k * cf);
    }
    double current(const Mat& ff, const Mat& dd, const double* cf, const double* gf) const {
        double v = ff[IP2][IP2] * cf[IP2] + dd[IP2][IP2] * gf[IP2];
        if (full_current) v = v + dd[IP2][IC] * gf[IC];  // D-11
        return v;
    }
    void solid_row(double i, const double di[NV], double dt, Mat& rj, double* g) const {
        for (int k = 0; k < NV; ++k) rj[ICS][k] = -(spec_a * di[k] / p.F);
        rj[ICS][ICS] = -(spec_a * di[ICS] / p.F) - 1.0 * vf_AM / dt;
        g[ICS] = +(spec_a * i / p.F);
    }

    void assemble(const std::vector<double>& c, double dt, std::vector<double>& A, std::vector<double>& B,
                  std::vector<double>& Dm, std::vector<double>& G) const {
        const int nj = p.nj;
        const double F = p.F, a = spec_a, eps = p.eps, eps_sep = p.eps_sep, sig = p.sigma;
        std::fill(A.begin(), A.end(), 0.0); std::fill(B.begin(), B.end(), 0.0);
        std::fill(Dm.begin(), Dm.end(), 0.0); std::fill(G.begin(), G.end(), 0.0);
        auto C = [&](int j, int k) { return c[static_cast<std::size_t>(j) * NV + k]; };
        for (int j = 0; j < nj; ++j) {
            Mat dW{}, dE{}, fW{}, fE{}, rj{};
            double g[NV] = {0, 0, 0, 0}, cW[NV] = {0, 0, 0, 0}, cE[NV] = {0, 0, 0, 0}, gW[NV] = {0, 0, 0, 0},
                   gE[NV] = {0, 0, 0, 0}, i, di[NV];
            for (int k = 0; k < NV; ++k) {
                if (j > 0) { cW[k] = aW[j] * C(j, k) + (1.0 - aW[j]) * C(j - 1, k); gW[k] = bW[j] * (C(j, k) - C(j - 1, k)); }
                if (j < nj - 1) { cE[k] = aE[j] * C(j + 1, k) + (1.0 - aE[j]) * C(j, k); gE[k] = bE[j] * (C(j + 1, k) - C(j, k)); }
            }
            rate_derivs(C(j, IC), C(j, ICS), C(j, IP1), C(j, IP2), i, di);
            const double dxj = dx[j];

            if (j == 0) {  // Li-foil face
                dE[IC][IC] = -(eps_sep_face * dcat_s);
                fE[IC][IC] = -(eps_sep_face * p.z_plus * ucat_s * F * gE[IP2]);
                dE[IC][IP2] = -(eps_sep_face * p.z_plus * ucat_s * F * cE[IC]);
                g[IC] = -i_app / F + (dE[IC][IC] * gE[IC] + fE[IC][IC] * cE[IC]);
                rj[ICS][ICS] = 0.0 - 1.0 * vf_AM / dt;
                dE[IP1][IP1] = -(1.0 - eps_sep_face) * sig;
                g[IP1] = phi1_sign * dE[IP1][IP1] * gE[IP1];
                rj[IP2][IP2] = 1.0;
                g[IP2] = 0.0 - C(j, IP2);
            } else if (j < s) {  // separator interior
                face(eps_sep, dcat_s, ucat_s, dan_s, uan_s, cW[IC], gW[IP2], dW, fW);
                face(eps_sep, dcat_s, ucat_s, dan_s, uan_s, cE[IC], gE[IP2], dE, fE);
                rj[IC][IC] = -(eps_sep / dt * dxj);
                g[IC] = 0.0 - (fW[IC][IC] * cW[IC] + dW[IC][IC] * gW[IC]) + (fE[IC][IC] * cE[IC] + dE[IC][IC] * gE[IC]);
                rj[ICS][ICS] = -(1.0 * (1.0 - eps_sep) / dt);
                dW[IP1][IP1] = -(1.0 - eps_sep) * sig;
                dE[IP1][IP1] = -(1.0 - eps_sep) * sig;
                g[IP1] = 0.0 - (fW[IP1][IP1] * cW[IP1] + dW[IP1][IP1] * gW[IP1]) + (fE[IP1][IP1] * cE[IP1] + dE[IP1][IP1] * gE[IP1]);
                g[IP2] = 0.0 - current(fW, dW, cW, gW) + current(fE, dE, cE, gE);
            } else if (j == s) {  // separator/cathode interface
                face(eps_sep_face, dcat_s, ucat_s, dan_s, uan_s, cW[IC], gW[IP2], dW, fW);
                face(eps, dcat_c, ucat_c, dan_c, uan_c, cE[IC], gE[IP2], dE, fE);
                g[IC] = 0.0 - (fW[IC][IC] * cW[IC] + dW[IC][IC] * gW[IC]) + (fE[IC][IC] * cE[IC] + dE[IC][IC] * gE[IC]);
                dW[IP1][IP1] = -(1.0 - eps_sep_face) * sig;
                dE[IP1][IP1] = -(1.0 - eps) * sig;
                g[IP1] = 0.0 - (fW[IP1][IP1] * cW[IP1] + dW[IP1][IP1] * gW[IP1]) + (fE[IP1][IP1] * cE[IP1] + dE[IP1][IP1] * gE[IP1]);
                g[IP2] = 0.0 - current(fW, dW, cW, gW) + current(fE, dE, cE, gE);
                solid_row(i, di, dt, rj, g);
            } else if (j < nj - 1) {  // cathode interior
                face(eps, dcat_c, ucat_c, dan_c, uan_c, cW[IC], gW[IP2], dW, fW);
                face(eps, dcat_c, ucat_c, dan_c, uan_c, cE[IC], gE[IP2], dE, fE);
                for (int k = 0; k < NV; ++k) {
                    rj[IC][k] = (a * di[k] / F) * dxj;
                    rj[IP1][k] = -((a * di[k]) * dxj);
                    rj[IP2][k] = (a * di[k]) * dxj;
                }
                rj[IC][IC] = (a * di[IC] / F) * dxj - (eps / dt) * dxj;
                g[IC] = -((a * i / F) * dxj) - (fW[IC][IC] * cW[IC] + dW[IC][IC] * gW[IC]) + (fE[IC][IC] * cE[IC] + dE[IC][IC] * gE[IC]);
                dW[IP1][IP1] = -(1.0 - eps) * sig;
                dE[IP1][IP1] = -(1.0 - eps) * sig;
                g[IP1] = (a * i) * dxj - (fW[IP1][IP1] * cW[IP1] + dW[IP1][IP1] * gW[IP1]) + (fE[IP1][IP1] * cE[IP1] + dE[IP1][IP1] * gE[IP1]);
                g[IP2] = -((a * i) * dxj) - current(fW, dW, cW, gW) + current(fE, dE, cE, gE);
                solid_row(i, di, dt, rj, g);
            } else {  // current collector
                face(eps, dcat_c, ucat_c, dan_c, uan_c, cW[IC], gW[IP2], dW, fW);
                g[IC] = 0.0 - dW[IC][IC] * gW[IC] - fW[IC][IC] * cW[IC];
                dW[IP1][IP1] = -(1.0 - eps) * sig;
                g[IP1] = i_app - dW[IP1][IP1] * gW[IP1];
                g[IP2] = 0.0 - dW[IP2][IC] * gW[IC] - fW[IP2][IC] * cW[IC];
                solid_row(i, di, dt, rj, g);
            }

            // control-volume blocks: A dc[j-1] + B dc[j] + D dc[j+1] = G
            const std::size_t o = static_cast<std::size_t>(j) * NV * NV;
            for (int r = 0; r < NV; ++r) {
                for (int k = 0; k < NV; ++k) {
                    const std::size_t q = o + r * NV + k;
                    if (j == 0) {
                        B[q] = rj[r][k] - (1.0 - aE[j]) * fE[r][k] + bE[j] * dE[r][k];
                        Dm[q] = -(aE[j] * fE[r][k]) - bE[j] * dE[r][k];
                    } else if (j == nj - 1) {
                        A[q] = (1.0 - aW[j]) * fW[r][k] - bW[j] * dW[r][k];
                        B[q] = rj[r][k] + bW[j] * dW[r][k] + aW[j] * fW[r][k];
                    } else {
                        A[q] = (1.0 - aW[j]) * fW[r][k] - bW[j] * dW[r][k];
                        B[q] = rj[r][k] + bW[j] * dW[r][k] + aW[j] * fW[r][k] - (1.0 - aE[j]) * fE[r][k] + bE[j] * dE[r][k];
                        Dm[q] = -(aE[j] * fE[r][k]) - bE[j] * dE[r][k];
                    }
                }
                G[static_cast<std::size_t>(j) * NV + r] = g[r];
            }
        }
    }

    // ---- corrected mode (log variables; docs/model.md section 10) ----
    // Butler-Volmer: i_n and d i_n / d(u, phi1, phi2, s). U is the arctangent fit plus the
    // thermodynamic tails (D-16) and the electrolyte Nernst term (D-14).
    void log_rate(const double* x, double& i, double di[NV]) const {
        const double rtf = p.R * p.T / p.F, ff = 1.0 / rtf, th = sigm(x[ICS]), om = sigm(-x[ICS]);
        const double z1 = -(ocp_c[2] * th) + ocp_c[3], z2 = -(ocp_c[5] * th) + ocp_c[6];
        const double fit = ocp_c[0] + ocp_c[1] * std::atan(z1) - ocp_c[4] * std::atan(z2);
        const double dfit = ocp_c[1] * (-ocp_c[2]) / (1.0 + z1 * z1) - ocp_c[4] * (-ocp_c[5]) / (1.0 + z2 * z2);
        const double theta_tail = 1.0e-4, se = std::log(theta_tail / (1.0 - theta_tail));
        const double tail = rtf * (softplus(se - x[ICS]) - softplus(x[ICS] + se));
        const double dtail = -rtf * (sigm(se - x[ICS]) + sigm(x[ICS] + se));
        const double dU_ds = dfit * th * om + dtail;
        const double ln_i0 = std::log(p.F * p.k_rxn * std::pow(p.c_bulk, p.alpha_a) * std::pow(cs_max(), p.alpha_a + p.alpha_c))
                             + p.alpha_a * x[IC] - p.alpha_a * softplus(x[ICS]) - p.alpha_c * softplus(-x[ICS]);
        const double eta = x[IP1] - x[IP2] - (fit + tail + rtf * x[IC]);
        const double ea = std::exp(ln_i0 + p.alpha_a * ff * eta), ec = std::exp(ln_i0 - p.alpha_c * ff * eta);
        i = ea - ec;
        const double di_deta = p.alpha_a * ff * ea + p.alpha_c * ff * ec;
        di[IC] = i * p.alpha_a - di_deta * rtf;
        di[IP1] = di_deta;
        di[IP2] = -di_deta;
        di[ICS] = i * (-p.alpha_a * th + p.alpha_c * om) - di_deta * dU_ds;
    }
    // face fluxes (N+, i1, i2) from xa to xb (Scharfetter-Gummel) and derivatives dFa, dFb [3][NV]
    void sg_face(const double* xa, const double* xb, double g, double gs, double* Fv, double* dFa, double* dFb) const {
        const double F = p.F, ff = F / (p.R * p.T);
        const double ca = p.c_bulk * std::exp(xa[IC]), cb = p.c_bulk * std::exp(xb[IC]), d = ff * (xb[IP2] - xa[IP2]);
        const double Bp = bern(d), Bm = bern(-d), dBp = bern_p(d), dBm = bern_p(-d);
        const double Np = g * dplus * (Bp * ca - Bm * cb), Nm = g * dminus * (Bm * ca - Bp * cb);
        const double dNp_dd = g * dplus * (dBp * ca + dBm * cb), dNm_dd = g * dminus * (-dBm * ca - dBp * cb);
        const double kb = g * p.kappa_bg;
        Fv[0] = Np;
        Fv[1] = -gs * (xb[IP1] - xa[IP1]);
        Fv[2] = F * (Np - Nm) - kb * (xb[IP2] - xa[IP2]);
        for (int z = 0; z < 3 * NV; ++z) { dFa[z] = 0.0; dFb[z] = 0.0; }
        dFa[0 * NV + IC] = g * dplus * Bp * ca;
        dFb[0 * NV + IC] = -g * dplus * Bm * cb;
        dFa[0 * NV + IP2] = -ff * dNp_dd;
        dFb[0 * NV + IP2] = ff * dNp_dd;
        dFa[1 * NV + IP1] = gs;
        dFb[1 * NV + IP1] = -gs;
        dFa[2 * NV + IC] = F * (dFa[0 * NV + IC] - g * dminus * Bm * ca);
        dFb[2 * NV + IC] = F * (dFb[0 * NV + IC] + g * dminus * Bp * cb);
        dFa[2 * NV + IP2] = -F * ff * (dNp_dd - dNm_dd) + kb;
        dFb[2 * NV + IP2] = F * ff * (dNp_dd - dNm_dd) - kb;
    }
    // residual and blocks (G = -R) in the log variables, with the particles at nodes s..nj-1
    void log_assemble(const std::vector<double>& x, const std::vector<double>& xold, double dt, double Ia,
                      std::vector<double>& A, std::vector<double>& B, std::vector<double>& Dm, std::vector<double>& G) const {
        const int nj = p.nj;
        const double F = p.F;
        std::fill(A.begin(), A.end(), 0.0); std::fill(B.begin(), B.end(), 0.0); std::fill(Dm.begin(), Dm.end(), 0.0);
        std::vector<double> R(static_cast<std::size_t>(nj) * NV, 0.0), Fv(3 * static_cast<std::size_t>(nj - 1)),
            dFa(3 * NV * static_cast<std::size_t>(nj - 1)), dFb(dFa.size());
        for (int k = 0; k < nj - 1; ++k) {
            const double h = (dx[k] + dx[k + 1]) / 2.0;
            const bool sep = k < s;
            const double g = (sep ? p.eps_sep / p.tau_sep : p.eps / tortuosity) / h;
            const double gs = (sep ? 1.0 - p.eps_sep : 1.0 - p.eps) * p.sigma / h;
            sg_face(&x[static_cast<std::size_t>(k) * NV], &x[static_cast<std::size_t>(k + 1) * NV], g, gs, &Fv[3 * k],
                    &dFa[3 * NV * k], &dFb[3 * NV * k]);
        }
        auto X = [&](int j, int v) { return x[static_cast<std::size_t>(j) * NV + v]; };
        auto blk = [&](std::vector<double>& M, int j, int r, int v) -> double& { return M[(static_cast<std::size_t>(j) * NV + r) * NV + v]; };
        auto FA = [&](int k, int f, int v) { return dFa[3 * NV * k + f * NV + v]; };
        auto FB = [&](int k, int f, int v) { return dFb[3 * NV * k + f * NV + v]; };
        auto Rr = [&](int j, int r) -> double& { return R[static_cast<std::size_t>(j) * NV + r]; };
        // foil face: N+ = I/F, zero electronic current, phi2 = 0 (the gauge)
        Rr(0, IC) = Fv[0] - Ia / F;
        for (int v = 0; v < NV; ++v) { blk(B, 0, IC, v) = FA(0, 0, v); blk(Dm, 0, IC, v) = FB(0, 0, v); }
        Rr(0, IP1) = X(1, IP1) - X(0, IP1);
        blk(B, 0, IP1, IP1) = -1.0; blk(Dm, 0, IP1, IP1) = 1.0;
        Rr(0, IP2) = X(0, IP2);
        blk(B, 0, IP2, IP2) = 1.0;
        // separator, interface and cathode
        for (int j = 1; j < nj - 1; ++j) {
            for (int f = 0; f < 3; ++f) {
                Rr(j, f) = Fv[3 * j + f] - Fv[3 * (j - 1) + f];
                for (int v = 0; v < NV; ++v) {
                    blk(B, j, f, v) += FA(j, f, v) - FB(j - 1, f, v);
                    blk(Dm, j, f, v) += FB(j, f, v);
                    blk(A, j, f, v) += -FA(j - 1, f, v);
                }
            }
            const double e = j < s ? p.eps_sep : p.eps;
            const double cj = p.c_bulk * std::exp(X(j, IC)), cjo = p.c_bulk * std::exp(xold[static_cast<std::size_t>(j) * NV + IC]);
            Rr(j, IC) += e * dx[j] * (cj - cjo) / dt;
            blk(B, j, IC, IC) += e * dx[j] * cj / dt;
        }
        // collector: no salt flux, no ionic current, electronic current I
        const int n = nj - 1, kl = nj - 2;
        Rr(n, IC) = -Fv[3 * kl];
        Rr(n, IP1) = Fv[3 * kl + 1] - Ia;
        Rr(n, IP2) = Fv[3 * kl + 2];
        for (int v = 0; v < NV; ++v) {
            blk(A, n, IC, v) = -FA(kl, 0, v); blk(B, n, IC, v) = -FB(kl, 0, v);
            blk(A, n, IP1, v) = FA(kl, 1, v); blk(B, n, IP1, v) = FB(kl, 1, v);
            blk(A, n, IP2, v) = FA(kl, 2, v); blk(B, n, IP2, v) = FB(kl, 2, v);
        }
        // the S column: fixed outside the cathode, the particles at s..nj-1
        for (int j = 0; j < nj; ++j) {
            Rr(j, ICS) = X(j, ICS) - xold[static_cast<std::size_t>(j) * NV + ICS];
            for (int v = 0; v < NV; ++v) blk(B, j, ICS, v) = 0.0;
            blk(B, j, ICS, ICS) = 1.0;
        }
        if (crystal) {  // the S column stays fixed; xtal_newton adds the crystals' reaction
            for (std::size_t q = 0; q < G.size(); ++q) G[q] = -R[q];
            return;
        }
        for (int j = s; j < nj; ++j) {
            double i, di[NV];
            log_rate(&x[static_cast<std::size_t>(j) * NV], i, di);
            const double a = spec_a;
            Rr(j, IC) -= a * i * dx[j] / F;
            Rr(j, IP1) += a * i * dx[j];
            Rr(j, IP2) -= a * i * dx[j];
            for (int v = 0; v < NV; ++v) {
                blk(B, j, IC, v) -= a * di[v] * dx[j] / F;
                blk(B, j, IP1, v) += a * di[v] * dx[j];
                blk(B, j, IP2, v) -= a * di[v] * dx[j];
            }
            const double th = sigm(X(j, ICS)), tho = sigm(xold[static_cast<std::size_t>(j) * NV + ICS]);
            Rr(j, ICS) = vf_AM * cs_max() * (th - tho) / dt + a * i / F;
            for (int v = 0; v < NV; ++v) blk(B, j, ICS, v) = a * di[v] / F;
            blk(B, j, ICS, ICS) += vf_AM * cs_max() * th * sigm(-X(j, ICS)) / dt;
        }
        for (std::size_t q = 0; q < G.size(); ++q) G[q] = -R[q];
    }
};

// Newton convergence: the update is below tol, or it has stagnated at the round-off floor
// (within 1e3*tol and down by less than half since the last iteration)
bool converged(double upd, double prev, double tol) {
    return upd <= tol || (upd <= 1.0e3 * tol && upd >= 0.5 * prev);
}

// Newton step length <= 1 limiting |du| <= 1, |dphi| <= 0.1 V and |ds| <= 2 per iteration
double lbounded(const std::vector<double>& d) {
    const double caps[NV] = {1.0, 0.1, 0.1, 2.0};
    double lam = 1.0;
    for (int k = 0; k < NV; ++k) {
        double mx = 0.0;
        for (std::size_t q = k; q < d.size(); q += NV) mx = std::max(mx, std::abs(d[q]));
        if (mx > caps[k]) lam = std::min(lam, caps[k] / mx);
    }
    return lam;
}

// the Newton update in the physical variables: max of e^u |du|, |dphi| and theta(1-theta) |ds|
double phys_update(const std::vector<double>& x, const std::vector<double>& d) {
    double u = 0.0;
    for (std::size_t q = 0; q < x.size(); q += NV)
        u = std::max({u, std::exp(x[q + IC]) * std::abs(d[q + IC]), std::abs(d[q + IP1]), std::abs(d[q + IP2]),
                      sigm(x[q + ICS]) * sigm(-x[q + ICS]) * std::abs(d[q + ICS])});
    return u;
}

// ============================== output ==============================
std::string fixed12(double v) {
    char b[64];
    if (std::isnan(v)) std::snprintf(b, sizeof b, "%12s", "NaN");
    else std::snprintf(b, sizeof b, "%12.5f", v);
    return b;
}
std::string sci15(double v) {
    char b[64];
    if (std::isnan(v)) std::snprintf(b, sizeof b, "%15s", "NaN");
    else std::snprintf(b, sizeof b, "%15.5E", v);
    return b;
}
std::string fit(const std::string& s, int w) {  // Fortran A<w>: right-justify, keep leftmost w characters
    const std::string t = s.substr(0, w);
    return std::string(w - t.size(), ' ') + t;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        Model m(read_input(argc > 1 ? argv[1] : "lfp.nml"));
        const Params& p = m.p;
        const int nj = p.nj;
        std::vector<Step> steps;
        if (!m.faithful) steps = parse_protocol(p.steps, p.C_rate, p.V_min, p.V_max, p.cycles);
        std::FILE* out = std::fopen(p.file.c_str(), "w");
        if (!out) throw std::runtime_error("cannot write " + p.file);

        std::vector<double> c(static_cast<std::size_t>(nj) * NV), dc(c.size(), 0.0), A(c.size() * NV), B(A.size()),
            Dm(A.size()), G(c.size());
        for (int j = 0; j < nj; ++j) {
            c[j * NV + IC] = p.c_bulk; c[j * NV + IP1] = p.phi1_init; c[j * NV + IP2] = p.phi2_init; c[j * NV + ICS] = p.cs_init;
        }
        if (!m.faithful) {  // corrected mode: u = ln(c/c_bulk) and the particles' log-odds
            const double th0 = p.cs_init / m.cs_max();
            for (int j = 0; j < nj; ++j) { c[j * NV + IC] = 0.0; c[j * NV + ICS] = std::log(th0 / (1.0 - th0)); }
        }
        const int nc = p.nj_crystal;
        std::vector<double> xc(static_cast<std::size_t>(nc) * m.nl, std::log((p.cs_init / m.cs_max()) / (1.0 - p.cs_init / m.cs_max())));
        double t = 0.0, mAhg = 0.0, dt = p.t_max / static_cast<double>(p.n_steps);
        std::string exit_reason = "max_steps";
        int nsolve = 0;

        auto header_lines = [&](bool extended) {
            std::vector<std::string> h1 = {"State", "Time", "Voltage", "Equivalence", "Anode_Eta", "anode_exchange_c", "Edge_c0"};
            std::vector<std::string> h2 = {"CDR", "hours", "Volts", "electron_equivs", "mV", "mA/cm2", "mol/cm3"};
            if (extended) { h1.push_back("Current"); h1.push_back("Step"); h1.push_back("Li_Nernst");
                          h2.push_back("mA/cm2"); h2.push_back("#"); h2.push_back("mV"); }
            for (auto* h : {&h1, &h2}) {
                std::string l = fit((*h)[0], 5) + " " + fit((*h)[1], 12) + " " + fit((*h)[2], 12);
                for (std::size_t k = 3; k < h->size(); ++k) l += " " + fit((*h)[k], 15);
                std::fprintf(out, "%s\n", l.c_str());
            }
        };
        // electrolyte concentration at the foil (corrected mode stores u = ln(c/c_bulk))
        auto c_foil = [&]() { return m.faithful ? c[IC] : p.c_bulk * std::exp(c[IC]); };
        auto i0_li = [&]() { return p.F * p.k_Li * std::pow(c_foil(), 0.5) * std::pow(p.c_Li_ref, 0.5); };
        // Li counter-electrode overpotential; corrected mode: symmetric Butler-Volmer (D-12)
        auto li_eta = [&](char state) {
            const double alpha = 0.5;
            if (!m.faithful) return -(p.R * p.T / (alpha * p.F)) * std::asinh(m.i_app / (2.0 * i0_li()));
            if (state == 'C') return 0.5 * std::log(m.i_app / i0_li()) / (alpha * p.F / (p.R * p.T));
            if (state == 'D') return -(0.5 * std::log(m.i_app / i0_li())) / (alpha * p.F / (p.R * p.T));
            return 0.0;
        };
        // corrected mode: Nernst potential of the lithium foil, (RT/F) ln(c/c_Li_ref)
        auto li_nernst = [&]() { return p.R * p.T / p.F * std::log(c_foil() / p.c_Li_ref); };
        // Corrected mode: cell voltage against the lithium foil (0 V). The solver fixes the gauge with
        // phi2 = 0 at the foil face; the equations depend only on potential differences, so the
        // foil-referenced potentials are the solved ones minus U_Li + eta_Li. (Imposing the foil
        // reference as the boundary condition makes the first block singular at rest.)
        auto cell_voltage = [&]() { return c[(nj - 1) * NV + IP1] + li_eta('D') - li_nernst(); };

        if (m.faithful) {
            // ---------------- the original program's constant-current discharge ----------------
            const char state = p.C_rate < 0 ? 'C' : 'D';
            const double write_every = p.t_max / p.n_steps / 200;
            long last_write = 0;  // an integer in the original (D-4)
            auto write_row = [&](bool header) {
                if (header) header_lines(false);
                const double eta = li_eta(state);
                std::fprintf(out, "%5c %s %s %s %s %s %s\n", state, fixed12(t / static_cast<double>(3600)).c_str(),
                             fixed12(c[(nj - 1) * NV + IP1] + eta).c_str(), sci15(mAhg * p.M * m.lit36 / p.F).c_str(),
                             sci15(eta * 1.0e3).c_str(), sci15(i0_li() * 1.0e3).c_str(), sci15(c[IC]).c_str());
            };
            for (int it = 1; it <= p.n_steps; ++it) {
                if (it == 1) {
                    write_row(true);
                } else if ((t - last_write) / 3600 >= write_every) {
                    write_row(false);
                    last_write = static_cast<long>(t - dt);
                } else if (it >= p.n_steps) {
                    write_row(false);
                } else if (c[(nj - 1) * NV + IP1] >= 99.0 && state == 'C') {
                    write_row(false); exit_reason = "end_of_charge"; break;
                } else if (std::isnan(dc[IC])) {
                    write_row(false); exit_reason = "nan"; break;
                } else if (t >= 99.0 * 3600.0) {
                    write_row(false); exit_reason = "max_time"; break;
                }
                if (state == 'D') mAhg = mAhg + 1000.0 * m.i_spec * dt / 3600.0;
                else mAhg = mAhg - 1000.0 * m.i_spec * dt / 3600.0;
                m.assemble(c, dt, A, B, Dm, G);
                if (!band::solve(NV, nj, A, B, Dm, G, dc, band::Pivot::legacy))
                    std::fill(dc.begin(), dc.end(), std::numeric_limits<double>::quiet_NaN());
                for (std::size_t k = 0; k < c.size(); ++k) c[k] = c[k] + dc[k];
                nsolve = it;
                dt = p.t_max / static_cast<double>(p.n_steps);
                t = t + dt;
            }
        } else {
            // ---------------- corrected mode: protocol of cc / cv / rest steps ----------------
            auto write_row = [&](bool header, int step) {
                if (header) header_lines(true);
                const char st = m.i_app > 0 ? 'D' : (m.i_app < 0 ? 'C' : 'R');
                const double eta = li_eta(st);
                std::fprintf(out, "%5c %s %s %s %s %s %s %s %15d %s\n", st, fixed12(t / 3600.0).c_str(),
                             fixed12(cell_voltage()).c_str(), sci15(mAhg * p.M * 3.6 / p.F).c_str(),
                             sci15(eta * 1.0e3).c_str(), sci15(i0_li() * 1.0e3).c_str(), sci15(c_foil()).c_str(),
                             sci15(m.i_app * 1.0e3).c_str(), step, sci15(li_nernst() * 1.0e3).c_str());
            };
            // scale every equation by the largest entry of its row in B (the solution is unchanged)
            auto equilibrate = [&]() {
                for (int j = 0; j < nj; ++j)
                    for (int r = 0; r < NV; ++r) {
                        const std::size_t o = (static_cast<std::size_t>(j) * NV + r) * NV;
                        double sc = 0.0;
                        for (int k = 0; k < NV; ++k) sc = std::max(sc, std::abs(B[o + k]));
                        if (sc == 0.0) sc = 1.0;
                        for (int k = 0; k < NV; ++k) { A[o + k] = A[o + k] / sc; B[o + k] = B[o + k] / sc; Dm[o + k] = Dm[o + k] / sc; }
                        G[static_cast<std::size_t>(j) * NV + r] = G[static_cast<std::size_t>(j) * NV + r] / sc;
                    }
            };
            // crystal l's rows without the reaction (diffusion_rows in Python): storage plus the net outward
            // diffusive flux of each control volume, and the tridiagonal dR/ds
            auto xtal_rows = [&](int l, double dtt, const std::vector<double>& xcold, std::vector<double>& Rc,
                                 std::vector<double>& Bc, std::vector<double>& Dq, std::vector<double>& Aq) {
                std::vector<double> th(nc), dth(nc), gx(nc - 1), Jf(nc - 1);
                const double csm = m.cs_max();
                const std::size_t o = static_cast<std::size_t>(l) * nc;
                for (int q = 0; q < nc; ++q) {
                    th[q] = sigm(xc[o + q]);
                    dth[q] = th[q] * sigm(-xc[o + q]);
                    Rc[q] = csm * m.xtal_V[q] * (th[q] - sigm(xcold[o + q])) / dtt;
                    Bc[q] = csm * m.xtal_V[q] * dth[q] / dtt;
                }
                for (int q = 0; q < nc - 1; ++q) {
                    gx[q] = p.D_c * csm * m.xtal_A[q] / m.xtal_h;
                    Jf[q] = -gx[q] * (th[q + 1] - th[q]);
                }
                std::fill(Dq.begin(), Dq.end(), 0.0);
                std::fill(Aq.begin(), Aq.end(), 0.0);
                for (int q = 0; q < nc - 1; ++q) {
                    Rc[q] = Rc[q] + Jf[q];
                    Bc[q] = Bc[q] + gx[q] * dth[q];
                    Dq[q] = -gx[q] * dth[q + 1];
                }
                for (int q = 1; q < nc; ++q) {
                    Rc[q] = Rc[q] - Jf[q - 1];
                    Bc[q] = Bc[q] + gx[q - 1] * dth[q];
                    Aq[q] = -gx[q - 1] * dth[q - 1];
                }
            };
            // one backward-Euler step with crystals, solved by the condensed Newton iteration (see
            // CrystalModel.newton_step in Python and xtal_newton in the Fortran program)
            auto xtal_newton = [&](double h) {
                const std::size_t ntot = static_cast<std::size_t>(nc) * m.nl;
                const std::vector<double> c_old = c, xc_old = xc;
                const double sgn[3] = {-1.0 / p.F, 1.0, -1.0};  // reaction sign in rows IC, IP1, IP2
                std::vector<double> gcpl(3 * static_cast<std::size_t>(m.nl)), Jce(gcpl.size());
                std::vector<double> Rc(nc), Bc(nc), Dq(nc), Aq(nc), A1(ntot), B1(ntot), D1(ntot), dxc(ntot);
                std::vector<std::vector<double>> rhs(4, std::vector<double>(ntot)), sol(4);
                double prev = std::numeric_limits<double>::infinity();
                for (int k = 0; k < p.newton_max_iter; ++k) {
                    m.log_assemble(c, c_old, h, m.i_app, A, B, Dm, G);
                    for (auto& r : rhs) std::fill(r.begin(), r.end(), 0.0);
                    for (int l = 1; l <= m.nl; ++l) {
                        const int j = m.s + l;
                        double y[NV], i, di[NV];
                        for (int v = 0; v < NV; ++v) y[v] = c[static_cast<std::size_t>(j) * NV + v];
                        y[ICS] = xc[static_cast<std::size_t>(l - 1) * nc + nc - 1];
                        m.log_rate(y, i, di);
                        const double w = m.a_x * m.dx[j];
                        for (int fl = 0; fl < 3; ++fl) {
                            G[static_cast<std::size_t>(j) * NV + fl] = G[static_cast<std::size_t>(j) * NV + fl] - (w * i) * sgn[fl];
                            for (int v = 0; v < 3; ++v) {
                                const std::size_t o = (static_cast<std::size_t>(j) * NV + fl) * NV + v;
                                B[o] = B[o] + (w * sgn[fl]) * di[v];
                            }
                            gcpl[static_cast<std::size_t>(l - 1) * 3 + fl] = (w * di[ICS]) * sgn[fl];
                            Jce[static_cast<std::size_t>(l - 1) * 3 + fl] = m.xtal_AR * di[fl] / p.F;
                        }
                        xtal_rows(l - 1, h, xc_old, Rc, Bc, Dq, Aq);
                        Rc[nc - 1] = Rc[nc - 1] + m.xtal_AR * i / p.F;
                        Bc[nc - 1] = Bc[nc - 1] + m.xtal_AR * di[ICS] / p.F;
                        for (int q = 0; q < nc; ++q) {
                            double sc = std::abs(Bc[q]);
                            if (sc == 0.0) sc = 1.0;
                            const std::size_t o = static_cast<std::size_t>(l - 1) * nc + q;
                            A1[o] = Aq[q] / sc;
                            B1[o] = Bc[q] / sc;
                            D1[o] = Dq[q] / sc;
                            rhs[0][o] = -Rc[q] / sc;
                            if (q == nc - 1)
                                for (int kk = 1; kk <= 3; ++kk) rhs[kk][o] = -Jce[static_cast<std::size_t>(l - 1) * 3 + kk - 1] / sc;
                        }
                    }
                    bool good = true;
                    for (int kk = 0; kk <= 3 && good; ++kk)
                        good = band::solve(1, static_cast<int>(ntot), A1, B1, D1, rhs[kk], sol[kk], band::Pivot::partial);
                    if (!good) break;
                    // condense the crystals into the electrode's diagonal blocks
                    for (int l = 1; l <= m.nl; ++l) {
                        const int j = m.s + l;
                        const std::size_t last = static_cast<std::size_t>(l) * nc - 1;
                        for (int fl = 0; fl < 3; ++fl) {
                            const double gc = gcpl[static_cast<std::size_t>(l - 1) * 3 + fl];
                            for (int kk = 0; kk < 3; ++kk) {
                                const std::size_t o = (static_cast<std::size_t>(j) * NV + fl) * NV + kk;
                                B[o] = B[o] + gc * sol[kk + 1][last];
                            }
                            G[static_cast<std::size_t>(j) * NV + fl] = G[static_cast<std::size_t>(j) * NV + fl] - gc * sol[0][last];
                        }
                    }
                    equilibrate();
                    if (!band::solve(NV, nj, A, B, Dm, G, dc, band::Pivot::partial)) break;
                    for (int l = 1; l <= m.nl; ++l) {
                        const std::size_t j = static_cast<std::size_t>(m.s + l) * NV;
                        for (int q = 0; q < nc; ++q) {
                            const std::size_t o = static_cast<std::size_t>(l - 1) * nc + q;
                            dxc[o] = sol[0][o] + sol[1][o] * dc[j + 0] + sol[2][o] * dc[j + 1] + sol[3][o] * dc[j + 2];
                        }
                    }
                    double lam = lbounded(dc), mx = 0.0;
                    for (double v : dxc) mx = std::max(mx, std::abs(v));
                    if (mx > 2.0) lam = std::min(lam, 2.0 / mx);
                    // divergence is judged on the damped electrode update: near theta = 0 or 1 a linearized
                    // log-odds update of the crystals is legitimately huge, and through the condensation it
                    // also inflates the undamped electrode update; the step limit scales both down
                    double raw = 0.0;
                    for (double v : dc) raw = std::max(raw, std::abs(v));
                    raw = lam * raw;
                    for (std::size_t q = 0; q < c.size(); ++q) c[q] = c[q] + lam * dc[q];
                    for (std::size_t q = 0; q < ntot; ++q) xc[q] = xc[q] + lam * dxc[q];
                    double upd = 0.0;
                    for (std::size_t q = 0; q < c.size(); q += NV)
                        upd = std::max({upd, std::exp(c[q + IC]) * std::abs(dc[q + IC]), std::abs(dc[q + IP1]), std::abs(dc[q + IP2])});
                    for (std::size_t q = 0; q < ntot; ++q) {
                        const double th = sigm(xc[q]);
                        upd = std::max(upd, th * sigm(-xc[q]) * std::abs(dxc[q]));
                    }
                    if (!std::isfinite(raw) || raw > 1.0e3) break;
                    if (lam == 1.0 && converged(upd, prev, p.newton_tol)) return true;
                    prev = upd;
                }
                c = c_old;
                xc = xc_old;
                return false;
            };
            // one backward-Euler step of length h at the current m.i_app, solved with Newton (log variables)
            auto newton_step = [&](double h) {
                if (m.crystal) return xtal_newton(h);
                const std::vector<double> c_old = c;
                double prev = std::numeric_limits<double>::infinity();
                for (int k = 0; k < p.newton_max_iter; ++k) {
                    m.log_assemble(c, c_old, h, m.i_app, A, B, Dm, G);
                    equilibrate();
                    if (!band::solve(NV, nj, A, B, Dm, G, dc, band::Pivot::partial)) break;
                    const double lam = lbounded(dc);
                    double raw = 0.0;
                    for (double v : dc) raw = std::max(raw, std::abs(v));
                    for (std::size_t q = 0; q < c.size(); ++q) c[q] = c[q] + lam * dc[q];
                    const double upd = phys_update(c, dc);
                    if (!std::isfinite(raw) || raw > 1.0e3) break;
                    if (lam == 1.0 && converged(upd, prev, p.newton_tol)) return true;
                    prev = upd;
                }
                c = c_old;
                return false;
            };
            // advance by dt; with check, locate a voltage-cutoff crossing to within 0.1 mV
            auto advance = [&](double dtt, double& t_done, bool& stopped, bool check, double vlo, double vhi) {
                const double min_dt = 1.0e-10, event_dv = 1.0e-4, event_min_dt = 1.0e-12;
                const int max_failures = 200;  // Newton failures allowed within one step
                int failures = 0;
                double hh = dtt;
                t_done = 0.0;
                stopped = false;
                while (t_done < dtt) {
                    hh = std::min(hh, dtt - t_done);
                    const std::vector<double> c_save = c, xc_save = xc;
                    if (!newton_step(hh)) {
                        if (hh / 2 < min_dt || ++failures >= max_failures) {
                            c = c_save;  // give up: keep the last converged sub-step (t_done advanced to it)
                            xc = xc_save;
                            return false;
                        }
                        hh = hh / 2;
                        continue;
                    }
                    if (check) {
                        // a discharge ends at vlo, a charge at vhi (the other bound is not its cutoff)
                        const double vv = cell_voltage();
                        const double mg = m.i_app > 0 ? vv - vlo : m.i_app < 0 ? vhi - vv : std::min(vv - vlo, vhi - vv);
                        if (mg < 0.0) {
                            if (mg < -event_dv && hh / 2 >= event_min_dt) { c = c_save; xc = xc_save; hh = hh / 2; continue; }
                            t_done = t_done + hh;
                            stopped = true;
                            return true;
                        }
                    }
                    t_done = t_done + hh;
                    hh = 2.0 * hh;  // grow back after a success (up to dt, by the min above)
                }
                return true;
            };
            // the physical limit the state has reached, reported as the exit reason when a step cannot be
            // solved: electrolyte below 1e-3*c_bulk anywhere, or particles within 1e-3 of full or empty
            auto limit_reason = [&]() -> std::string {
                double cmin = std::numeric_limits<double>::infinity(), thmin = cmin, thmax = -cmin;
                for (int j = 0; j < nj; ++j) cmin = std::min(cmin, p.c_bulk * std::exp(c[j * NV + IC]));
                if (m.crystal) {  // the crystal surfaces, where the reaction is
                    for (int l = 0; l < m.nl; ++l) {
                        thmin = std::min(thmin, sigm(xc[static_cast<std::size_t>(l) * nc + nc - 1]));
                        thmax = std::max(thmax, sigm(xc[static_cast<std::size_t>(l) * nc + nc - 1]));
                    }
                } else {
                    for (int j = m.s; j < nj; ++j) {
                        thmin = std::min(thmin, sigm(c[j * NV + ICS]));
                        thmax = std::max(thmax, sigm(c[j * NV + ICS]));
                    }
                }
                if (cmin < 1.0e-3 * p.c_bulk) return "electrolyte_depleted";
                if (thmax > 1.0 - 1.0e-3) return "particles_full";
                if (thmin < 1.0e-3) return "particles_empty";
                return "solver_fail";
            };
            // one constant-voltage time step: find I with V(I) = V_set (see simulate.cv_step in Python)
            auto cv_step = [&](double h, double V_set, double& I) {
                const double tol = 1.0e-9, inf = std::numeric_limits<double>::infinity();
                const std::vector<double> c_start = c, xc_start = xc;
                auto f = [&](double Itry, bool& good) {
                    c = c_start;
                    xc = xc_start;
                    m.i_app = Itry;
                    good = newton_step(h);
                    if (!good) return Itry < 0 ? inf : -inf;
                    return cell_voltage() - V_set;
                };
                bool good;
                double fI = f(I, good);
                if (good && std::abs(fI) <= tol) return true;
                double grow = std::max(std::abs(I), 1.0e-2 * m.i_1C), a = 0, b = 0, fa = 0, fb = 0;
                bool have_a = false, have_b = false;
                for (int it = 0; it < 60; ++it) {
                    if (fI > 0) {
                        a = I; fa = fI; have_a = true;
                        if (have_b) break;
                        I = I < 0 ? 0.0 : I + grow;
                    } else {
                        b = I; fb = fI; have_b = true;
                        if (have_a) break;
                        I = I > 0 ? 0.0 : I - grow;
                    }
                    grow *= 2.0;
                    fI = f(I, good);
                    if (good && std::abs(fI) <= tol) return true;
                }
                if (!(have_a && have_b)) { c = c_start; xc = xc_start; return false; }
                int side = 0;
                for (int it = 0; it < 200; ++it) {
                    if (std::isfinite(fa) && std::isfinite(fb)) {
                        I = (a * fb - b * fa) / (fb - fa);
                        if (!(a < I && I < b)) I = 0.5 * (a + b);
                    } else {
                        I = 0.5 * (a + b);
                    }
                    fI = f(I, good);
                    if (std::abs(fI) <= tol || (b - a) <= 1.0e-14 * m.i_1C) {
                        // a collapsed bracket can sit on the edge of the currents for which a step converges,
                        // next to a state far from V_set: accept only a state on the set voltage
                        const bool accept = good && std::abs(fI) <= 1.0e-6;
                        if (!accept) { c = c_start; xc = xc_start; }
                        return accept;
                    }
                    if (fI > 0) {
                        a = I; fa = fI;
                        if (side == 1 && std::isfinite(fb)) fb *= 0.5;
                        side = 1;
                    } else {
                        b = I; fb = fI;
                        if (side == -1 && std::isfinite(fa)) fa *= 0.5;
                        side = -1;
                    }
                }
                c = c_start;
                xc = xc_start;
                return false;
            };
            // a constant-voltage sub-step: cv_step over h, halved (down to 1e-10 s) while no current holds V_set
            // for that long; on failure c, xc and I are as on entry and h_done = 0
            auto cv_advance = [&](double h, double V_set, double& I, double& h_done) {
                const double I_guess = I;
                const std::vector<double> c_begin = c, xc_begin = xc;
                h_done = h;
                while (true) {
                    c = c_begin;
                    xc = xc_begin;
                    I = I_guess;
                    if (cv_step(h_done, V_set, I)) return true;
                    if (h_done / 2 < 1.0e-10) break;
                    h_done = h_done / 2;
                }
                c = c_begin;
                xc = xc_begin;
                I = I_guess;
                h_done = 0.0;
                return false;
            };

            double I = steps[0].kind == Kind::cc ? steps[0].C * m.i_1C : 0.0;
            m.i_app = I;
            write_row(true, 1);
            double last_write = t;
            std::string reason;
            int n_done = 0;
            bool finished = true;
            for (std::size_t k = 0; k < steps.size() && finished; ++k) {
                const Step& st = steps[k];
                const int step_no = static_cast<int>(k) + 1;
                double t_step = 0.0;
                if (st.kind == Kind::cc) I = st.C * m.i_1C;
                else if (st.kind == Kind::rest) I = 0.0;
                while (true) {
                    const double h = st.t < 0 ? dt : std::min(dt, st.t - t_step);
                    double h_done = 0.0;
                    bool stopped = false, ok;
                    std::string why;
                    if (st.kind == Kind::cv) {
                        ok = cv_advance(h, st.V, I, h_done);
                        stopped = st.Imin >= 0 && std::abs(I) <= st.Imin * m.i_1C;
                        why = "current_limit";
                    } else {
                        m.i_app = I;
                        ok = advance(h, h_done, stopped, st.kind == Kind::cc, st.Vmin, st.Vmax);
                        why = stopped && cell_voltage() <= st.Vmin ? "cutoff_low" : "cutoff_high";
                    }
                    m.i_app = I;
                    if (!ok) {
                        // keep the sub-steps completed before the failure: the limit is judged where it was reached
                        if (st.kind != Kind::cv && h_done > 0.0) {
                            mAhg = mAhg + 1000.0 * (I / m.mass_area) * h_done / 3600.0;
                            t = t + h_done;
                            ++n_done;
                        }
                        write_row(false, step_no); exit_reason = limit_reason(); finished = false; break;
                    }
                    mAhg = mAhg + 1000.0 * (I / m.mass_area) * h_done / 3600.0;
                    t = t + h_done;
                    t_step = t_step + h_done;
                    ++n_done;
                    bool bad = false;
                    for (double x : c) bad = bad || std::isnan(x);
                    for (double x : xc) bad = bad || std::isnan(x);
                    if (bad) { write_row(false, step_no); exit_reason = "nan"; finished = false; break; }
                    if (stopped || (st.t >= 0 && t_step >= st.t * (1.0 - 1.0e-12))) {
                        write_row(false, step_no);
                        last_write = t;
                        reason = stopped ? why : "duration";
                        break;
                    }
                    if (t - last_write >= p.write_interval) { write_row(false, step_no); last_write = t; }
                    if (t >= 99.0 * 3600.0) { write_row(false, step_no); exit_reason = "max_time"; finished = false; break; }
                }
            }
            if (finished) exit_reason = steps.size() == 1 ? reason : "end_of_protocol";
            nsolve = n_done;
        }
        std::fclose(out);
        std::printf("%s run, C-rate %.17g: exit %s after %d steps; wrote %s\n", p.mode.c_str(), p.C_rate,
                    exit_reason.c_str(), nsolve, p.file.c_str());
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 2;
    }
}
