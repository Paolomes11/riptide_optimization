/* Copyright 2026 Giulio Mesini — Apache License 2.0 */

#ifndef RIPTIDE_EXP_COMMON_FITS_IO_HPP
#define RIPTIDE_EXP_COMMON_FITS_IO_HPP

/**
 * fits_io — Lettura di file .fit / .fits a 16 bit (grayscale)
 *
 * Implementazione minimale FITS senza dipendenze esterne (cfitsio).
 * Supporta:
 *   - SIMPLE FITS con NAXIS=2 (immagini 2D grayscale)
 *   - BITPIX = 16 (int16 big-endian con BSCALE/BZERO)
 *   - Header FITS standard a blocchi di 2880 byte
 */

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace riptide::fits {

/// Struttura header FITS (subset keyword standard utilizzate in RIPTIDE).
struct FitsHeader {
  int naxis1    = 0;
  int naxis2    = 0;
  int bitpix    = 0;
  double bscale = 1.0;
  double bzero  = 0.0;
  std::string instrument;
  std::string object;
  std::string date_obs;
  double exptime = 0.0;
};

/// Frame FITS: header + dati pixel in layout row-major.
struct FitsFrame {
  FitsHeader header;
  std::vector<float> data;

  /// Accesso sicuro al pixel con controllo bounds.
  double pixel(int x, int y) const {
    if (x < 0 || x >= header.naxis1 || y < 0 || y >= header.naxis2)
      return 0.0;
    return static_cast<double>(
        data[static_cast<size_t>(y) * static_cast<size_t>(header.naxis1) + static_cast<size_t>(x)]);
  }

  /// Larghezza immagine [pixel].
  int width() const {
    return header.naxis1;
  }

  /// Altezza immagine [pixel].
  int height() const {
    return header.naxis2;
  }

  /// Numero totale di pixel.
  size_t npixels() const {
    return static_cast<size_t>(header.naxis1) * static_cast<size_t>(header.naxis2);
  }
};

/**
 * Legge un singolo file .fit / .fits e restituisce un FitsFrame.
 *
 * @param path  Percorso al file FITS
 * @return      FitsFrame con header e dati pixel
 * @throws      std::runtime_error se il file non è leggibile o il formato non è supportato
 */
FitsFrame read_fits(const std::filesystem::path& path);

/**
 * Restituisce tutti i path .fit/.fits in una directory, ordinati.
 *
 * @param dir  Cartella contenente i file FITS
 * @return     Lista di file FITS ordinata lessicograficamente
 * @throws     std::runtime_error se la cartella non esiste o non contiene file FITS
 */
std::vector<std::filesystem::path> list_fits_files(const std::filesystem::path& dir);

/**
 * Legge tutti i file .fit / .fits in una directory (ordine lessicografico).
 *
 * @param dir         Cartella contenente i file FITS
 * @param max_frames  Numero massimo di frame da leggere (0 = tutti)
 * @return            Vettore di FitsFrame
 * @throws            std::runtime_error se la cartella non esiste o non contiene file FITS
 */
std::vector<FitsFrame> read_fits_stack(const std::filesystem::path& dir, size_t max_frames = 0);

} // namespace riptide::fits

#endif // RIPTIDE_EXP_COMMON_FITS_IO_HPP
