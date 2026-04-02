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

#ifndef EXP1_FIT_READER_HPP
#define EXP1_FIT_READER_HPP

/**
 * fit_reader — Lettura di file .fit / .fits a 16 bit (grayscale)
 *
 * Implementazione minimale FITS senza dipendenze esterne (cfitsio).
 * Supporta:
 *   - SIMPLE FITS con NAXIS=2 (immagini 2D grayscale)
 *   - BITPIX = 16 (uint16 big-endian, valore massimo 65535 = 2^16 - 1)
 *   - Header FITS standard a blocchi di 2880 byte
 *
 * I dati vengono restituiti come matrice flat row-major (y * width + x)
 * con valori in virgola mobile [0.0, 65535.0].
 *
 * Limitazioni intenzionali (non necessarie per questo progetto):
 *   - Non supporta BSCALE/BZERO (assumiamo dati raw dal sensore)
 *   - Non supporta immagini 3D o cubi
 *   - Non supporta BITPIX = 8 o 32
 */

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace exp1 {

// Struttura header FITS

struct FitsHeader {
  int naxis1    = 0;      // width  (numero di colonne, asse NAXIS1)
  int naxis2    = 0;      // height (numero di righe,   asse NAXIS2)
  int bitpix    = 0;      // bit per pixel (16 per uint16)
  double bscale = 1.0;    // scala applicata ai dati raw
  double bzero  = 0.0;    // offset applicato ai dati raw
  std::string instrument; // INSTRUME keyword (opzionale)
  std::string object;     // OBJECT keyword (opzionale)
  std::string date_obs;   // DATE-OBS keyword (opzionale)
  double exptime = 0.0;   // EXPTIME keyword [s] (opzionale)
};

// Frame: header + dati pixel

struct FitsFrame {
  FitsHeader header;

  // Dati pixel in virgola mobile (float per risparmiare memoria), layout row-major:
  //   pixel(x, y) = data[y * header.naxis1 + x]
  // Range fisico: [0.0, 65535.0] per sensori a 16 bit senza BSCALE/BZERO.
  std::vector<float> data;

  // Accesso sicuro con controllo bounds
  double pixel(int x, int y) const {
    if (x < 0 || x >= header.naxis1 || y < 0 || y >= header.naxis2)
      return 0.0;
    return static_cast<double>(
        data[static_cast<size_t>(y) * static_cast<size_t>(header.naxis1) + static_cast<size_t>(x)]);
  }

  int width() const {
    return header.naxis1;
  }
  int height() const {
    return header.naxis2;
  }
  size_t npixels() const {
    return static_cast<size_t>(header.naxis1) * static_cast<size_t>(header.naxis2);
  }
};

// API pubblica

/**
 * Legge un singolo file .fit / .fits e restituisce un FitsFrame.
 *
 * @param path  Percorso al file FITS
 * @return      FitsFrame con header e dati pixel
 * @throws      std::runtime_error se il file non è leggibile o il formato non è supportato
 */
FitsFrame read_fits(const std::filesystem::path& path);

/**
 * Legge tutti i file .fit / .fits in una directory (ordine lessicografico).
 *
 * @param dir      Cartella contenente i file FITS
 * @param max_frames  Numero massimo di frame da leggere (0 = tutti)
 * @return         Vettore di FitsFrame
 * @throws         std::runtime_error se la cartella non esiste o non contiene file FITS
 */
std::vector<FitsFrame> read_fits_stack(const std::filesystem::path& dir, size_t max_frames = 0);

/**
 * Restituisce tutti i path .fit/.fits in una directory, ordinati.
 */
std::vector<std::filesystem::path> list_fits_files(const std::filesystem::path& dir);

} // namespace exp1

#endif // EXP1_FIT_READER_HPP