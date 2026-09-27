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

struct Params {
    double L_cath = 24.0e-4, L_sep = 25.0e-4;
    int nj = 101, sep_node = 22;
    double eps = 0.5, eps_AM = 0.8, eps_sep = 0.39, tau_sep = 4.0, bruggeman = -0.5;
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
        {"l_cath", &p.L_cath}, {"l_sep", &p.L_sep}, {"eps", &p.eps}, {"eps_am", &p.eps_AM},
        {"eps_sep", &p.eps_sep}, {"tau_sep", &p.tau_sep}, {"bruggeman", &p.bruggeman}, {"d", &p.D},
        {"t_plus", &p.t_plus}, {"c_bulk", &p.c_bulk}, {"z_plus", &p.z_plus}, {"z_minus", &p.z_minus},
        {"sigma", &p.sigma}, {"m", &p.M}, {"rho", &p.rho}, {"q_th", &p.Q_th}, {"r_p", &p.R_p},
        {"k_rxn", &p.k_rxn}, {"alpha_a", &p.alpha_a}, {"alpha_c", &p.alpha_c}, {"k_li", &p.k_Li},
        {"c_li_ref", &p.c_Li_ref}, {"r", &p.R}, {"t", &p.T}, {"f", &p.F}, {"c_rate", &p.C_rate},
        {"phi1_init", &p.phi1_init}, {"phi2_init", &p.phi2_init}, {"cs_init", &p.cs_init},
        {"t_max", &p.t_max}, {"fd_step", &p.fd_step}, {"v_min", &p.V_min}, {"v_max", &p.V_max},
        {"newton_tol", &p.newton_tol}, {"write_interval", &p.write_interval}};
    std::map<std::string, int*> ints = {{"nj", &p.nj}, {"sep_node", &p.sep_node}, {"n_steps", &p.n_steps},
                                        {"newton_max_iter", &p.newton_max_iter}, {"cycles", &p.cycles}};
    std::map<std::string, std::string*> strs = {{"mode", &p.mode}, {"file", &p.file}, {"steps", &p.steps}};
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
constexpr double THETA_REG = 1.0e-6;  // D-13 regularization threshold

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
    double lit36, ocp_c[7], spec_a, tortuosity, i_spec, eps_sep_face, phi1_sign, mass_area, i_1C;
    double i_app;  // applied current density [A/cm2]; changes step by step in corrected mode
    bool full_current;
    double dcat_s, dan_s, ucat_s, uan_s, dcat_c, dan_c, ucat_c, uan_c;
    int s;  // 0-based interface node
    std::vector<double> dx, aW, aE, bW, bE;

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
        spec_a = 3 * p.eps_AM / p.R_p;
        tortuosity = std::pow(p.eps, p.bruggeman);
        i_spec = p.Q_th * p.C_rate;
        i_app = i_spec * p.L_cath * p.eps_AM * p.rho;
        mass_area = p.L_cath * p.eps_AM * p.rho;
        i_1C = p.Q_th * mass_area;
        const double t_an = 1.0 - p.t_plus;
        const double d_cat = p.D * (1.0 + (t_an / p.t_plus)) / (2.0 * t_an / p.t_plus);
        const double d_an = d_cat * t_an / p.t_plus;
        const double u_cat = d_cat / (p.R * p.T), u_an = d_an / (p.R * p.T);
        dcat_s = d_cat / p.tau_sep; dan_s = d_an / p.tau_sep; ucat_s = u_cat / p.tau_sep; uan_s = u_an / p.tau_sep;
        dcat_c = d_cat / tortuosity; dan_c = d_an / tortuosity; ucat_c = u_cat / tortuosity; uan_c = u_an / tortuosity;

        const int nj = p.nj;
        s = p.sep_node - 1;
        dx.assign(nj, 0.0); aW.assign(nj, 0.0); aE.assign(nj, 0.0); bW.assign(nj, 0.0); bE.assign(nj, 0.0);
        const double h_sep = p.L_sep / static_cast<double>(p.sep_node - 2);
        const double h_cat = p.L_cath / static_cast<double>(nj - p.sep_node - 1);
        for (int j = 1; j < s; ++j) dx[j] = h_sep;
        for (int j = s + 1; j < nj - 1; ++j) dx[j] = h_cat;
        for (int j = 1; j < nj; ++j) { aW[j] = dx[j - 1] / (dx[j - 1] + dx[j]); bW[j] = 2.0 / (dx[j - 1] + dx[j]); }
        for (int j = 0; j < nj - 1; ++j) { aE[j] = dx[j] / (dx[j + 1] + dx[j]); bE[j] = 2.0 / (dx[j] + dx[j + 1]); }
    }

    // ---- kinetics ----
    double cs_max() const { return (p.rho / p.M) * p.M * p.Q_th * 1000.0 * lit36 / p.F; }
    double ocp(double cs) const {
        const double th = (cs / (p.rho / p.M)) / (p.M * p.Q_th * 1000.0 * lit36 / p.F);
        return ocp_c[0] + ocp_c[1] * std::atan(-(ocp_c[2] * th) + ocp_c[3]) - ocp_c[4] * std::atan(-(ocp_c[5] * th) + ocp_c[6]);
    }
    double rate(double c, double cs, double p1, double p2) const {
        const double eta = p1 - p2 - ocp(cs);
        double i0 = p.F * p.k_rxn * std::pow(c, p.alpha_a) * std::pow(cs_max() - cs, p.alpha_a) * std::pow(cs, p.alpha_c);
        if (faithful) i0 = r32(i0);  // D-4
        return i0 * (std::exp(p.alpha_a * p.F * eta / (p.R * p.T)) - std::exp(-(p.alpha_c * p.F * eta / (p.R * p.T))));
    }
    // rate and finite-difference derivatives w.r.t. (c, phi1, phi2, cs) (D-6)
    double ocp_slope(double cs) const {  // dU/dcs [V cm3/mol]
        const double th = (cs / (p.rho / p.M)) / (p.M * p.Q_th * 1000.0 * lit36 / p.F);
        const double dth = 1.0 / ((p.rho / p.M) * (p.M * p.Q_th * 1000.0 * lit36 / p.F));
        const double x1 = -(ocp_c[2] * th) + ocp_c[3], x2 = -(ocp_c[5] * th) + ocp_c[6];
        const double du = ocp_c[1] * (-ocp_c[2]) / (1.0 + x1 * x1) - ocp_c[4] * (-ocp_c[5]) / (1.0 + x2 * x2);
        return du * dth;
    }
    // rate and exact derivatives; corrected mode (fixes D-6)
    // x^alpha, replaced below delta by a C1 quadratic with g(0) = 0 and a finite slope (D-13)
    static void power_reg(double x, double alpha, double delta, double& g, double& dg) {
        if (x < delta) {
            const double u = x / delta;
            g = std::pow(delta, alpha) * ((2.0 - alpha) * u + (alpha - 1.0) * u * u);
            dg = std::pow(delta, alpha - 1.0) * ((2.0 - alpha) + 2.0 * (alpha - 1.0) * u);
        } else {
            g = std::pow(x, alpha);
            dg = alpha * std::pow(x, alpha - 1.0);
        }
    }
    void rate_derivs_exact(double c, double cs, double p1, double p2, double& i, double di[NV]) const {
        const double rt = p.R * p.T, aa = p.alpha_a * p.F / rt, bb = p.alpha_c * p.F / rt;
        const double eta = p1 - p2 - ocp(cs);
        double gv, dgv, gs, dgs;
        power_reg(cs_max() - cs, p.alpha_a, THETA_REG * cs_max(), gv, dgv);
        power_reg(cs, p.alpha_c, THETA_REG * cs_max(), gs, dgs);
        const double pre = p.F * p.k_rxn * std::pow(c, p.alpha_a);
        const double i0 = pre * gv * gs, di0 = pre * (gs * -dgv + gv * dgs);
        const double ea = std::exp(aa * eta), ec = std::exp(-bb * eta);
        i = i0 * (ea - ec);
        const double di_deta = i0 * (aa * ea + bb * ec);
        di[IC] = p.alpha_a * i / c;
        di[ICS] = di0 * (ea - ec) - di_deta * ocp_slope(cs);
        di[IP1] = di_deta;
        di[IP2] = -di_deta;
    }
    void rate_derivs(double c, double cs, double p1, double p2, double& i, double di[NV]) const {
        if (!faithful) { rate_derivs_exact(c, cs, p1, p2, i, di); return; }
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
    std::vector<double> time_terms(double dt) const {
        const int nj = p.nj;
        std::vector<double> T(static_cast<std::size_t>(nj) * NV, 0.0);
        T[ICS] = -(p.eps_AM / dt);
        for (int j = 1; j < s; ++j) { T[j * NV + IC] = -(p.eps_sep / dt * dx[j]); T[j * NV + ICS] = -((1.0 - p.eps_sep) / dt); }
        for (int j = s + 1; j < nj - 1; ++j) T[j * NV + IC] = -((p.eps / dt) * dx[j]);
        for (int j = s; j < nj; ++j) T[j * NV + ICS] = -(p.eps_AM / dt);
        return T;
    }
    // largest step <= 1 keeping 0 < c and 0 < cs < cs_max (at most 90 % of the way to a bound)
    // and changing no potential by more than 0.1 V (Butler-Volmer exponentials make Newton overshoot)
    double bounded_step(const std::vector<double>& c, const std::vector<double>& d) const {
        const double keep = 0.9, csm = cs_max(), max_dphi = 0.1;  // potential change cap per iteration [V]
        double lam = 1.0, dphi = 0.0;
        for (int j = 0; j < p.nj; ++j)
            dphi = std::max({dphi, std::abs(d[j * NV + IP1]), std::abs(d[j * NV + IP2])});
        if (dphi > max_dphi) lam = max_dphi / dphi;
        for (int j = 0; j < p.nj; ++j) {
            const double dcc = d[j * NV + IC], dcs = d[j * NV + ICS];
            if (dcc < 0) lam = std::min(lam, keep * (c[j * NV + IC] - 0.0) / -dcc);
            if (dcs < 0) lam = std::min(lam, keep * (c[j * NV + ICS] - 0.0) / -dcs);
            if (dcs > 0) lam = std::min(lam, keep * (csm - c[j * NV + ICS]) / dcs);
        }
        return std::max(lam, 0.0);
    }

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
        rj[ICS][ICS] = -(spec_a * di[ICS] / p.F) - 1.0 * p.eps_AM / dt;
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
                rj[ICS][ICS] = 0.0 - 1.0 * p.eps_AM / dt;
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
};

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
        double t = 0.0, mAhg = 0.0, dt = p.t_max / static_cast<double>(p.n_steps);
        std::string exit_reason = "max_steps";
        int nsolve = 0;

        auto header_lines = [&](bool extended) {
            std::vector<std::string> h1 = {"State", "Time", "Voltage", "Equivalence", "Anode_Eta", "anode_exchange_c", "Edge_c0"};
            std::vector<std::string> h2 = {"CDR", "hours", "Volts", "electron_equivs", "mV", "mA/cm2", "mol/cm3"};
            if (extended) { h1.push_back("Current"); h1.push_back("Step"); h2.push_back("mA/cm2"); h2.push_back("#"); }
            for (auto* h : {&h1, &h2}) {
                std::string l = fit((*h)[0], 5) + " " + fit((*h)[1], 12) + " " + fit((*h)[2], 12);
                for (std::size_t k = 3; k < h->size(); ++k) l += " " + fit((*h)[k], 15);
                std::fprintf(out, "%s\n", l.c_str());
            }
        };
        auto i0_li = [&]() { return p.F * p.k_Li * std::pow(c[IC], 0.5) * std::pow(p.c_Li_ref, 0.5); };
        // Li counter-electrode overpotential; corrected mode: symmetric Butler-Volmer (D-12)
        auto li_eta = [&](char state) {
            const double alpha = 0.5;
            if (!m.faithful) return -(p.R * p.T / (alpha * p.F)) * std::asinh(m.i_app / (2.0 * i0_li()));
            if (state == 'C') return 0.5 * std::log(m.i_app / i0_li()) / (alpha * p.F / (p.R * p.T));
            if (state == 'D') return -(0.5 * std::log(m.i_app / i0_li())) / (alpha * p.F / (p.R * p.T));
            return 0.0;
        };
        auto cell_voltage = [&]() { return c[(nj - 1) * NV + IP1] + li_eta('D'); };

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
                std::fprintf(out, "%5c %s %s %s %s %s %s %s %15d\n", st, fixed12(t / 3600.0).c_str(),
                             fixed12(c[(nj - 1) * NV + IP1] + eta).c_str(), sci15(mAhg * p.M * 3.6 / p.F).c_str(),
                             sci15(eta * 1.0e3).c_str(), sci15(i0_li() * 1.0e3).c_str(), sci15(c[IC]).c_str(),
                             sci15(m.i_app * 1.0e3).c_str(), step);
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
            // one backward-Euler step of length h at the current m.i_app, solved with Newton
            auto newton_step = [&](double h) {
                const std::vector<double> c_old = c, T = m.time_terms(h);
                const double scale[NV] = {p.c_bulk, 1.0, 1.0, m.cs_max()};
                for (int k = 0; k < p.newton_max_iter; ++k) {
                    m.assemble(c, h, A, B, Dm, G);
                    for (std::size_t q = 0; q < G.size(); ++q) G[q] = G[q] - T[q] * (c[q] - c_old[q]);
                    equilibrate();
                    if (!band::solve(NV, nj, A, B, Dm, G, dc, band::Pivot::partial)) break;
                    const double lam = m.bounded_step(c, dc);
                    double worst = 0.0;
                    for (std::size_t q = 0; q < c.size(); ++q) {
                        c[q] = c[q] + lam * dc[q];
                        worst = std::max(worst, std::abs(dc[q]) / scale[q % NV]);
                    }
                    if (lam == 1.0 && worst <= p.newton_tol) return true;
                }
                c = c_old;
                return false;
            };
            // advance by dt; with check, locate a voltage-cutoff crossing to within 0.1 mV
            auto advance = [&](double dtt, double& t_done, bool& stopped, bool check, double vlo, double vhi) {
                const double min_dt = 1.0e-6, event_dv = 1.0e-4, event_min_dt = 1.0e-9;
                double hh = dtt;
                t_done = 0.0;
                stopped = false;
                while (t_done < dtt) {
                    hh = std::min(hh, dtt - t_done);
                    const std::vector<double> c_save = c;
                    if (!newton_step(hh)) {
                        if (hh / 2 < min_dt) return false;
                        hh = hh / 2;
                        continue;
                    }
                    if (check) {
                        const double vv = cell_voltage(), mg = std::min(vv - vlo, vhi - vv);
                        if (mg < 0.0) {
                            if (mg < -event_dv && hh / 2 >= event_min_dt) { c = c_save; hh = hh / 2; continue; }
                            t_done = t_done + hh;
                            stopped = true;
                            return true;
                        }
                    }
                    t_done = t_done + hh;
                }
                return true;
            };
            // one constant-voltage time step: find I with V(I) = V_set (see simulate.cv_step in Python)
            auto cv_step = [&](double h, double V_set, double& I) {
                const double tol = 1.0e-9, inf = std::numeric_limits<double>::infinity();
                const std::vector<double> c_start = c;
                auto f = [&](double Itry, bool& good) {
                    c = c_start;
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
                if (!(have_a && have_b)) { c = c_start; return false; }
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
                        if (!good) c = c_start;
                        return good;
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
                        ok = cv_step(h, st.V, I);
                        h_done = h;
                        stopped = st.Imin >= 0 && std::abs(I) <= st.Imin * m.i_1C;
                        why = "current_limit";
                    } else {
                        m.i_app = I;
                        ok = advance(h, h_done, stopped, st.kind == Kind::cc, st.Vmin, st.Vmax);
                        why = stopped && cell_voltage() <= st.Vmin ? "cutoff_low" : "cutoff_high";
                    }
                    m.i_app = I;
                    if (!ok) { write_row(false, step_no); exit_reason = "solver_fail"; finished = false; break; }
                    mAhg = mAhg + 1000.0 * (I / m.mass_area) * h_done / 3600.0;
                    t = t + h_done;
                    t_step = t_step + h_done;
                    ++n_done;
                    bool bad = false;
                    for (double x : c) bad = bad || std::isnan(x);
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
