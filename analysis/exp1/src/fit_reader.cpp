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

#include "fit_reader.hpp"

#include <algorithm>
#include <array>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace exp1 {

// Costanti FITS

static constexpr size_t FITS_BLOCK_SIZE  = 2880; // byte per blocco header/data
static constexpr size_t FITS_CARD_SIZE   = 80;   // byte per card header
static constexpr size_t FITS_CARDS_BLOCK = FITS_BLOCK_SIZE / FITS_CARD_SIZE; // 36

// Helpers big-endian

// Converte 2 byte big-endian in int16 (poi cast a uint16)
static int16_t be_to_int16(const uint8_t* p) {
  return static_cast<int16_t>((static_cast<uint16_t>(p[0]) << 8) | static_cast<uint16_t>(p[1]));
}

// Parsing card FITS

// Estrae il valore intero da una card "KEYWORD =      VALORE / commento"
static bool card_int(const std::string& card, const std::string& key, int& out) {
  if (card.substr(0, key.size()) != key)
    return false;
  std::string val = card.substr(10, 20);
  try {
    out = std::stoi(val);
    return true;
  } catch (...) {
    return false;
  }
}

static bool card_double(const std::string& card, const std::string& key, double& out) {
  if (card.substr(0, key.size()) != key)
    return false;
  std::string val = card.substr(10, 20);
  // rimuovi spazi
  val.erase(std::remove(val.begin(), val.end(), ' '), val.end());
  try {
    out = std::stod(val);
    return true;
  } catch (...) {
    return false;
  }
}

static bool card_string(const std::string& card, const std::string& key, std::string& out) {
  if (card.substr(0, key.size()) != key)
    return false;
  // Il valore stringa è tra apici singoli: '...'
  size_t q1 = card.find('\'');
  size_t q2 = card.find('\'', q1 + 1);
  if (q1 == std::string::npos || q2 == std::string::npos)
    return false;
  out = card.substr(q1 + 1, q2 - q1 - 1);
  // trim trailing spaces
  while (!out.empty() && out.back() == ' ')
    out.pop_back();
  return true;
}

// Lettura header FITS

static FitsHeader parse_header(std::ifstream& f) {
  FitsHeader hdr;
  bool end_found = false;
  int naxis      = 0;

  while (!end_found && f.good()) {
    // Leggi un blocco di 2880 byte (36 card da 80 byte)
    std::array<char, FITS_BLOCK_SIZE> block{};
    f.read(block.data(), static_cast<std::streamsize>(FITS_BLOCK_SIZE));
    if (f.gcount() < static_cast<std::streamsize>(FITS_BLOCK_SIZE))
      throw std::runtime_error("fit_reader: file FITS troncato nell'header");

    for (size_t c = 0; c < FITS_CARDS_BLOCK; ++c) {
      // Ogni card è esattamente 80 caratteri ASCII (senza null terminator nello standard)
      std::string card(block.data() + c * FITS_CARD_SIZE, FITS_CARD_SIZE);

      if (card.substr(0, 3) == "END") {
        end_found = true;
        break;
      }

      // Parsing keyword per keyword
      int ival    = 0;
      double dval = 0.0;
      std::string sval;

      if (card_int(card, "BITPIX  ", ival))
        hdr.bitpix = ival;
      else if (card_int(card, "NAXIS   ", ival))
        naxis = ival;
      else if (card_int(card, "NAXIS1  ", ival))
        hdr.naxis1 = ival;
      else if (card_int(card, "NAXIS2  ", ival))
        hdr.naxis2 = ival;
      else if (card_double(card, "BSCALE  ", dval))
        hdr.bscale = dval;
      else if (card_double(card, "BZERO   ", dval))
        hdr.bzero = dval;
      else if (card_double(card, "EXPTIME ", dval))
        hdr.exptime = dval;
      else if (card_string(card, "INSTRUME", sval))
        hdr.instrument = sval;
      else if (card_string(card, "OBJECT  ", sval))
        hdr.object = sval;
      else if (card_string(card, "DATE-OBS", sval))
        hdr.date_obs = sval;
    }
  }

  if (!end_found)
    throw std::runtime_error("fit_reader: keyword END non trovato nell'header FITS");
  if (naxis != 2)
    throw std::runtime_error("fit_reader: NAXIS=" + std::to_string(naxis)
                             + " non supportato (richiesto NAXIS=2)");
  if (hdr.bitpix != 16)
    throw std::runtime_error("fit_reader: BITPIX=" + std::to_string(hdr.bitpix)
                             + " non supportato (richiesto BITPIX=16)");
  if (hdr.naxis1 <= 0 || hdr.naxis2 <= 0)
    throw std::runtime_error("fit_reader: dimensioni immagine non valide (NAXIS1="
                             + std::to_string(hdr.naxis1) + ", NAXIS2=" + std::to_string(hdr.naxis2)
                             + ")");
  return hdr;
}

// Lettura dati pixel

static std::vector<float> read_data(std::ifstream& f, const FitsHeader& hdr) {
  const size_t n = static_cast<size_t>(hdr.naxis1) * static_cast<size_t>(hdr.naxis2);

  // I dati FITS a 16 bit sono int16 big-endian con offset BZERO = 32768
  // per rappresentare uint16. La convenzione standard è:
  //   valore_fisico = BSCALE * raw + BZERO
  // Per sensori a 16 bit unsigned raw senza BSCALE/BZERO personalizzati:
  //   raw = int16 (firmato), valore_fisico = raw + 32768 → [0, 65535]

  // Leggi raw bytes
  std::vector<uint8_t> raw(n * 2);
  f.read(reinterpret_cast<char*>(raw.data()), static_cast<std::streamsize>(n * 2));
  if (static_cast<size_t>(f.gcount()) < n * 2)
    throw std::runtime_error("fit_reader: file FITS troncato nei dati pixel");

  std::vector<float> data(n);
  // FITS è big-endian. Il valore uint16 si ottiene da int16 + 32768 (con BZERO=32768)
  // oppure direttamente interpretando come uint16 big-endian.
  // Per sicurezza usiamo BSCALE e BZERO dall'header.
  for (size_t i = 0; i < n; ++i) {
    int16_t raw16 = be_to_int16(raw.data() + i * 2);
    // Applica BSCALE e BZERO standard
    double val = hdr.bscale * static_cast<double>(raw16) + hdr.bzero;
    // Clamp a [0, 65535] per robustezza
    if (val < 0.0)
      val = 0.0;
    if (val > 65535.0)
      val = 65535.0;
    data[i] = static_cast<float>(val);
  }

  return data;
}

// Implementazioni API pubbliche

FitsFrame read_fits(const std::filesystem::path& path) {
  std::ifstream f(path, std::ios::binary);
  if (!f.is_open())
    throw std::runtime_error("fit_reader: impossibile aprire " + path.string());

  FitsFrame frame;
  frame.header = parse_header(f);
  frame.data   = read_data(f, frame.header);

  std::cout << "[FIT] Letto: " << path.filename().string() << "  (" << frame.header.naxis1 << "×"
            << frame.header.naxis2 << ", BITPIX=" << frame.header.bitpix << ")\n";

  return frame;
}

std::vector<std::filesystem::path> list_fits_files(const std::filesystem::path& dir) {
  if (!std::filesystem::exists(dir) || !std::filesystem::is_directory(dir))
    throw std::runtime_error("fit_reader: cartella non trovata: " + dir.string());

  std::vector<std::filesystem::path> files;
  for (const auto& entry : std::filesystem::directory_iterator(dir)) {
    if (!entry.is_regular_file())
      continue;
    auto ext = entry.path().extension().string();
    // normalizza in minuscolo per confronto
    std::transform(ext.begin(), ext.end(), ext.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (ext == ".fit" || ext == ".fits" || ext == ".fts")
      files.push_back(entry.path());
  }

  std::sort(files.begin(), files.end());

  if (files.empty())
    throw std::runtime_error("fit_reader: nessun file .fit/.fits trovato in " + dir.string());

  return files;
}

std::vector<FitsFrame> read_fits_stack(const std::filesystem::path& dir, size_t max_frames) {
  auto files = list_fits_files(dir);

  if (max_frames > 0 && files.size() > max_frames)
    files.resize(max_frames);

  std::cout << "[FIT] Lettura stack da: " << dir.string() << " (" << files.size() << " frame)\n";

  std::vector<FitsFrame> stack;
  stack.reserve(files.size());

  for (const auto& p : files)
    stack.push_back(read_fits(p));

  // Verifica consistenza dimensioni
  if (stack.size() > 1) {
    int w0 = stack[0].width(), h0 = stack[0].height();
    for (size_t i = 1; i < stack.size(); ++i) {
      if (stack[i].width() != w0 || stack[i].height() != h0)
        throw std::runtime_error("fit_reader: dimensioni incoerenti nel frame " + std::to_string(i)
                                 + " (" + std::to_string(stack[i].width()) + "×"
                                 + std::to_string(stack[i].height()) + " vs " + std::to_string(w0)
                                 + "×" + std::to_string(h0) + ")");
    }
  }

  return stack;
}

} // namespace exp1