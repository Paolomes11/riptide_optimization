/*
 * Copyright 2026 Giulio Mesini
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 */

#include "homography.hpp"

#include <nlohmann/json.hpp>

#include <array>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>

namespace riptide::exp3 {

// ---------------------------------------------------------------------------
// Jacobi eigendecomposition per una matrice simmetrica 9×9.
// Restituisce gli eigenvalue in 'evals' e gli eigenvector come colonne di 'evecs'.
// ---------------------------------------------------------------------------
static void jacobi_eigen9(std::array<double, 81>& A,  // input/output (distrutta)
                           std::array<double, 9>&  evals,
                           std::array<double, 81>& evecs) {
    // Inizializza evecs come identità
    evecs.fill(0.0);
    for (int i = 0; i < 9; ++i)
        evecs[static_cast<size_t>(i) * 9 + static_cast<size_t>(i)] = 1.0;

    const int MAX_SWEEP = 200;
    const double TOL    = 1e-10;

    for (int sweep = 0; sweep < MAX_SWEEP; ++sweep) {
        // Controlla convergenza: max |off-diagonal|
        double off = 0.0;
        for (int i = 0; i < 9; ++i)
            for (int j = i + 1; j < 9; ++j)
                off = std::max(off, std::abs(A[static_cast<size_t>(i) * 9 + static_cast<size_t>(j)]));
        if (off < TOL)
            break;

        // Spazza tutti i piani (i,j) con i<j
        for (int p = 0; p < 9; ++p) {
            for (int q = p + 1; q < 9; ++q) {
                double apq = A[static_cast<size_t>(p) * 9 + static_cast<size_t>(q)];
                if (std::abs(apq) < 1e-15)
                    continue;
                double app = A[static_cast<size_t>(p) * 9 + static_cast<size_t>(p)];
                double aqq = A[static_cast<size_t>(q) * 9 + static_cast<size_t>(q)];
                double theta = 0.5 * (aqq - app) / apq;
                double t     = (theta >= 0.0)
                               ? 1.0 / (theta + std::sqrt(1.0 + theta * theta))
                               : 1.0 / (theta - std::sqrt(1.0 + theta * theta));
                double c = 1.0 / std::sqrt(1.0 + t * t);
                double s = t * c;

                // Aggiorna A con la rotazione di Givens G(p,q,theta)
                A[static_cast<size_t>(p) * 9 + static_cast<size_t>(p)] = app - t * apq;
                A[static_cast<size_t>(q) * 9 + static_cast<size_t>(q)] = aqq + t * apq;
                A[static_cast<size_t>(p) * 9 + static_cast<size_t>(q)] = 0.0;
                A[static_cast<size_t>(q) * 9 + static_cast<size_t>(p)] = 0.0;

                for (int r = 0; r < 9; ++r) {
                    if (r == p || r == q)
                        continue;
                    double arp = A[static_cast<size_t>(r) * 9 + static_cast<size_t>(p)];
                    double arq = A[static_cast<size_t>(r) * 9 + static_cast<size_t>(q)];
                    A[static_cast<size_t>(r) * 9 + static_cast<size_t>(p)] =  c * arp - s * arq;
                    A[static_cast<size_t>(p) * 9 + static_cast<size_t>(r)] =  c * arp - s * arq;
                    A[static_cast<size_t>(r) * 9 + static_cast<size_t>(q)] =  s * arp + c * arq;
                    A[static_cast<size_t>(q) * 9 + static_cast<size_t>(r)] =  s * arp + c * arq;
                }

                // Aggiorna gli eigenvector
                for (int r = 0; r < 9; ++r) {
                    double erp = evecs[static_cast<size_t>(r) * 9 + static_cast<size_t>(p)];
                    double erq = evecs[static_cast<size_t>(r) * 9 + static_cast<size_t>(q)];
                    evecs[static_cast<size_t>(r) * 9 + static_cast<size_t>(p)] =  c * erp - s * erq;
                    evecs[static_cast<size_t>(r) * 9 + static_cast<size_t>(q)] =  s * erp + c * erq;
                }
            }
        }
    }

    for (int i = 0; i < 9; ++i)
        evals[static_cast<size_t>(i)] = A[static_cast<size_t>(i) * 9 + static_cast<size_t>(i)];
}

// ---------------------------------------------------------------------------
// Normalizzazione isotropica dei punti (Hartley 1997).
// Restituisce la matrice di trasformazione T (3×3, row-major, array<double,9>).
// ---------------------------------------------------------------------------
static std::array<double, 9> normalize_points(std::vector<CalibPoint>& pts, bool use_disp) {
    double cx = 0.0, cy = 0.0;
    for (const auto& p : pts) {
        cx += use_disp ? p.disp_x : p.sens_x;
        cy += use_disp ? p.disp_y : p.sens_y;
    }
    double n = static_cast<double>(pts.size());
    cx /= n;
    cy /= n;

    double scale = 0.0;
    for (const auto& p : pts) {
        double dx = (use_disp ? p.disp_x : p.sens_x) - cx;
        double dy = (use_disp ? p.disp_y : p.sens_y) - cy;
        scale += std::sqrt(dx * dx + dy * dy);
    }
    scale /= n;
    if (scale < 1e-12)
        scale = 1.0;
    scale = std::sqrt(2.0) / scale;

    // Applica trasformazione
    for (auto& p : pts) {
        if (use_disp) {
            p.disp_x = (p.disp_x - cx) * scale;
            p.disp_y = (p.disp_y - cy) * scale;
        } else {
            p.sens_x = (p.sens_x - cx) * scale;
            p.sens_y = (p.sens_y - cy) * scale;
        }
    }

    // T = [[scale, 0, -scale*cx], [0, scale, -scale*cy], [0, 0, 1]]
    std::array<double, 9> T{};
    T[0] = scale;  T[1] = 0.0;   T[2] = -scale * cx;
    T[3] = 0.0;    T[4] = scale; T[5] = -scale * cy;
    T[6] = 0.0;    T[7] = 0.0;   T[8] = 1.0;
    return T;
}

// ---------------------------------------------------------------------------
// Moltiplicazione matrice 3×3 (row-major): C = A * B
// ---------------------------------------------------------------------------
static std::array<double, 9> mat3_mul(const std::array<double, 9>& A,
                                      const std::array<double, 9>& B) {
    std::array<double, 9> C{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                C[static_cast<size_t>(i) * 3 + static_cast<size_t>(j)] +=
                    A[static_cast<size_t>(i) * 3 + static_cast<size_t>(k)] *
                    B[static_cast<size_t>(k) * 3 + static_cast<size_t>(j)];
    return C;
}

// ---------------------------------------------------------------------------
// Inversa di matrice di normalizzazione T (struttura diagonale con traslazione).
// T^{-1} = [[1/s, 0, cx], [0, 1/s, cy], [0, 0, 1]]
// ---------------------------------------------------------------------------
static std::array<double, 9> norm_mat_inv(const std::array<double, 9>& T) {
    double s = T[0];
    std::array<double, 9> Ti{};
    Ti[0] = 1.0 / s; Ti[1] = 0.0;     Ti[2] = -T[2] / s;
    Ti[3] = 0.0;     Ti[4] = 1.0 / s; Ti[5] = -T[5] / s;
    Ti[6] = 0.0;     Ti[7] = 0.0;     Ti[8] = 1.0;
    return Ti;
}

// ---------------------------------------------------------------------------
// Proietta un punto display con H e restituisce le coordinate nel sensore.
// ---------------------------------------------------------------------------
static void project_point(const std::array<double, 9>& H,
                           double xd, double yd,
                           double& xs_out, double& ys_out) {
    double w  = H[6] * xd + H[7] * yd + H[8];
    xs_out = (H[0] * xd + H[1] * yd + H[2]) / w;
    ys_out = (H[3] * xd + H[4] * yd + H[5]) / w;
}

// ---------------------------------------------------------------------------
Homography compute_homography(const std::vector<CalibPoint>& points) {
    if (points.size() < 4)
        throw std::invalid_argument("compute_homography: servono almeno 4 punti (n=" +
                                    std::to_string(points.size()) + ")");

    // Copia lavorativa per la normalizzazione
    std::vector<CalibPoint> pts = points;

    auto T_d = normalize_points(pts, /*use_disp=*/true);
    auto T_s = normalize_points(pts, /*use_disp=*/false);

    // Costruisce A^T·A (9×9)
    std::array<double, 81> ATA{};
    ATA.fill(0.0);

    for (const auto& p : pts) {
        double xd = p.disp_x, yd = p.disp_y;
        double xs = p.sens_x, ys = p.sens_y;

        // Riga 1: [-xd, -yd, -1, 0, 0, 0, xs*xd, xs*yd, xs]
        // Riga 2: [0, 0, 0, -xd, -yd, -1, ys*xd, ys*yd, ys]
        double r1[9] = {-xd, -yd, -1.0,  0.0,  0.0,  0.0,  xs * xd, xs * yd, xs};
        double r2[9] = { 0.0,  0.0,  0.0, -xd, -yd, -1.0,  ys * xd, ys * yd, ys};

        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 9; ++j) {
                ATA[static_cast<size_t>(i) * 9 + static_cast<size_t>(j)] +=
                    r1[i] * r1[j] + r2[i] * r2[j];
            }
    }

    // Jacobi eigendecomposition
    std::array<double, 9>  evals{};
    std::array<double, 81> evecs{};
    jacobi_eigen9(ATA, evals, evecs);

    // Trova l'indice dell'eigenvalue minimo
    int min_idx = 0;
    for (int i = 1; i < 9; ++i)
        if (evals[static_cast<size_t>(i)] < evals[static_cast<size_t>(min_idx)])
            min_idx = i;

    // Estrai eigenvector corrispondente (colonna min_idx di evecs)
    std::array<double, 9> h{};
    for (int i = 0; i < 9; ++i)
        h[static_cast<size_t>(i)] = evecs[static_cast<size_t>(i) * 9 + static_cast<size_t>(min_idx)];

    // Reshape in H 3×3 e denormalizza: H_final = T_s^{-1} * H_norm * T_d
    auto T_s_inv = norm_mat_inv(T_s);
    auto H_final = mat3_mul(mat3_mul(T_s_inv, h), T_d);

    // Normalizza per H[8] = 1
    if (std::abs(H_final[8]) > 1e-12) {
        double w = H_final[8];
        for (auto& v : H_final)
            v /= w;
    }

    // Calcola RMS residui sui punti originali
    double rms = 0.0;
    for (const auto& p : points) {
        double xs_proj, ys_proj;
        project_point(H_final, p.disp_x, p.disp_y, xs_proj, ys_proj);
        double dx = xs_proj - p.sens_x;
        double dy = ys_proj - p.sens_y;
        rms += dx * dx + dy * dy;
    }
    rms = std::sqrt(rms / static_cast<double>(points.size()));

    Homography result;
    result.H            = H_final;
    result.rms_residual = rms;
    result.n_points     = static_cast<int>(points.size());
    return result;
}

PhysicalPoint apply_homography(const Homography& H, double disp_px_x, double disp_px_y,
                               double sensor_pixel_size_mm) {
    double xs, ys;
    project_point(H.H, disp_px_x, disp_px_y, xs, ys);
    return {xs * sensor_pixel_size_mm, ys * sensor_pixel_size_mm};
}

void save_homography(const Homography& H, const std::filesystem::path& path) {
    nlohmann::json j;
    j["H"]   = std::vector<double>(H.H.begin(), H.H.end());
    j["rms"] = H.rms_residual;
    j["n"]   = H.n_points;
    std::filesystem::create_directories(path.parent_path());
    std::ofstream ofs(path);
    if (!ofs)
        throw std::runtime_error("save_homography: impossibile aprire " + path.string());
    ofs << j.dump(4) << '\n';
}

Homography load_homography(const std::filesystem::path& path) {
    std::ifstream ifs(path);
    if (!ifs)
        throw std::runtime_error("load_homography: file non trovato: " + path.string());
    nlohmann::json j;
    ifs >> j;
    Homography H;
    auto hv = j["H"].get<std::vector<double>>();
    if (hv.size() != 9)
        throw std::runtime_error("load_homography: vettore H non ha 9 elementi");
    std::copy(hv.begin(), hv.end(), H.H.begin());
    H.rms_residual = j.value("rms", 0.0);
    H.n_points     = j.value("n", 0);
    return H;
}

} // namespace riptide::exp3
