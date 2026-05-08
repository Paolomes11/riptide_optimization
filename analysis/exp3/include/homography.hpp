/* Copyright 2026 Giulio Mesini — Apache License 2.0 */

#ifndef RIPTIDE_EXP3_HOMOGRAPHY_HPP
#define RIPTIDE_EXP3_HOMOGRAPHY_HPP

#include <array>
#include <filesystem>
#include <stdexcept>
#include <vector>

namespace riptide::exp3 {

/// Omografia proiettiva 3×3 (row-major): H[3*i+j] = elemento (i,j).
struct Homography {
  std::array<double, 9> H{};
  double rms_residual = 0.0;
  int    n_points     = 0;
};

/// Corrispondenza pixel display → pixel sensore per la calibrazione.
struct CalibPoint {
  double disp_x = 0.0;
  double disp_y = 0.0;
  double sens_x = 0.0;
  double sens_y = 0.0;
};

/// Punto fisico nel piano sensore [mm].
struct PhysicalPoint {
  double x_mm = 0.0;
  double y_mm = 0.0;
};

/**
 * Calcola l'omografia H con DLT (Direct Linear Transform) da almeno 4 punti.
 *
 * Implementazione self-contained: costruisce la matrice A (2N×9), calcola
 * A^T·A (9×9) e trova l'eigenvector con eigenvalue minimo via Jacobi.
 * I punti sono normalizzati isotropicamente prima del calcolo (stabilità numerica).
 *
 * @param points  Lista di corrispondenze display→sensore (almeno 4)
 * @return        Omografia calibrata con RMS residui
 * @throws        std::invalid_argument se n_points < 4
 */
Homography compute_homography(const std::vector<CalibPoint>& points);

/**
 * Applica H per trasformare coordinate pixel display in coordinate fisiche sensore.
 *
 * La trasformazione è:
 *   [wx, wy, w]^T = H · [disp_px_x, disp_px_y, 1]^T
 *   pixel_sens_x  = wx/w,  pixel_sens_y = wy/w
 *   x_mm = pixel_sens_x * sensor_pixel_size_mm
 *   y_mm = pixel_sens_y * sensor_pixel_size_mm
 *
 * @param H                    Omografia calibrata
 * @param disp_px_x            Coordinata x pixel nel display
 * @param disp_px_y            Coordinata y pixel nel display
 * @param sensor_pixel_size_mm Dimensione fisica pixel sensore [mm/pixel]
 * @return                     Coordinata fisica nel piano sensore [mm]
 */
PhysicalPoint apply_homography(const Homography& H, double disp_px_x, double disp_px_y,
                               double sensor_pixel_size_mm);

/**
 * Serializza l'omografia in JSON: {"H":[h00,...,h22],"rms":...,"n":...}.
 *
 * @param H    Omografia da serializzare
 * @param path Path del file JSON di output
 */
void save_homography(const Homography& H, const std::filesystem::path& path);

/**
 * Deserializza l'omografia da un file JSON prodotto da save_homography.
 *
 * @param path Path del file JSON
 * @return     Omografia deserializzata
 * @throws     std::runtime_error se il file non esiste o il formato è errato
 */
Homography load_homography(const std::filesystem::path& path);

} // namespace riptide::exp3

#endif // RIPTIDE_EXP3_HOMOGRAPHY_HPP
