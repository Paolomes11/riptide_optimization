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

#include <filesystem>
#include <string>
#include <vector>

#include "fits_io.hpp"

namespace exp1 {

using FitsHeader = riptide::fits::FitsHeader;
using FitsFrame  = riptide::fits::FitsFrame;

// API pubblica

/**
 * Legge un singolo file .fit / .fits e restituisce un FitsFrame.
 *
 * @param path  Percorso al file FITS
 * @return      FitsFrame con header e dati pixel
 * @throws      std::runtime_error se il file non è leggibile o il formato non è supportato
 */
using riptide::fits::read_fits;

/**
 * Legge tutti i file .fit / .fits in una directory (ordine lessicografico).
 *
 * @param dir      Cartella contenente i file FITS
 * @param max_frames  Numero massimo di frame da leggere (0 = tutti)
 * @return         Vettore di FitsFrame
 * @throws         std::runtime_error se la cartella non esiste o non contiene file FITS
 */
using riptide::fits::read_fits_stack;

/**
 * Restituisce tutti i path .fit/.fits in una directory, ordinati.
 */
using riptide::fits::list_fits_files;

} // namespace exp1

#endif // EXP1_FIT_READER_HPP
