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

#ifndef RIPTIDE_FOCUS_MAP_HPP
#define RIPTIDE_FOCUS_MAP_HPP

#include <cmath>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace riptide {

// Chiave per lookup (x1,x2): arrotonda a 1 decimale per robustezza floating-point
inline std::string make_focus_key(double x1, double x2) {
  auto r = [](double v) -> long { return static_cast<long>(std::round(v * 10.0)); };
  return std::to_string(r(x1)) + "_" + std::to_string(r(x2));
}

// Mappa (x1,x2)->x_focus caricata dal TSV prodotto da dof_map
using FocusMap = std::unordered_map<std::string, double>;

// Legge il file TSV e restituisce la mappa.
// Cerca le colonne "x1", "x2", "x_focus" dall'header; lancia runtime_error se assenti.
inline FocusMap load_focus_map(const std::string& tsv_path) {
  std::ifstream in(tsv_path);
  if (!in.is_open())
    throw std::runtime_error("focus_map: impossibile aprire " + tsv_path);

  std::string line;
  if (!std::getline(in, line))
    throw std::runtime_error("focus_map: file vuoto: " + tsv_path);

  // Trova indici colonne dall'header
  std::vector<std::string> headers;
  {
    std::istringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, '\t'))
      headers.push_back(tok);
  }
  auto find_col = [&](const std::string& name) -> int {
    for (int i = 0; i < static_cast<int>(headers.size()); ++i)
      if (headers[i] == name)
        return i;
    throw std::runtime_error("focus_map: colonna '" + name + "' non trovata in " + tsv_path);
  };
  int col_x1     = find_col("x1");
  int col_x2     = find_col("x2");
  int col_xfocus = find_col("x_focus");

  FocusMap m;
  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    std::vector<std::string> fields;
    {
      std::istringstream ss(line);
      std::string tok;
      while (std::getline(ss, tok, '\t'))
        fields.push_back(tok);
    }
    int ncols = static_cast<int>(fields.size());
    if (col_x1 >= ncols || col_x2 >= ncols || col_xfocus >= ncols)
      continue;
    double x1     = std::stod(fields[col_x1]);
    double x2     = std::stod(fields[col_x2]);
    double xfocus = std::stod(fields[col_xfocus]);
    m[make_focus_key(x1, x2)] = xfocus;
  }
  return m;
}

// Lookup con tolleranza: cerca la chiave nella mappa; restituisce nullopt se assente.
// Il confronto è già incorporato nella chiave (arrotondamento a 0.1 mm).
inline std::optional<double> lookup_focus(const FocusMap& m, double x1, double x2) {
  auto it = m.find(make_focus_key(x1, x2));
  if (it == m.end())
    return std::nullopt;
  return it->second;
}

// Estrae tutte le coppie (x1,x2) dalla mappa, ordinate.
// Usato in modalità focus-tsv per generare le coppie da simulare direttamente dal TSV.
inline std::vector<std::pair<double, double>> get_pairs_from_focus_map(const FocusMap& m) {
  std::vector<std::pair<double, double>> pairs;
  pairs.reserve(m.size());
  for (const auto& [key, _] : m) {
    auto sep = key.find('_');
    if (sep == std::string::npos)
      continue;
    double x1 = static_cast<double>(std::stol(key.substr(0, sep))) / 10.0;
    double x2 = static_cast<double>(std::stol(key.substr(sep + 1))) / 10.0;
    pairs.push_back({x1, x2});
  }
  std::sort(pairs.begin(), pairs.end());
  return pairs;
}

// Mappa lens-pair -> FocusMap, per TSV con colonne extra lens1_id/lens2_id
// (fuoco thin-lens accurato per singola coppia, invece di un'unica curva
// condivisa/mediata su tutte le coppie).
using PerPairFocusMap = std::unordered_map<std::string, FocusMap>;

inline std::string make_pair_key(const std::string& l1_id, const std::string& l2_id) {
  return l1_id + "|" + l2_id;
}

// True se l'header del TSV contiene le colonne lens1_id/lens2_id, cioè se il
// file è nel formato per-coppia invece del formato (x1,x2,x_focus) legacy.
inline bool focus_tsv_has_pair_columns(const std::string& tsv_path) {
  std::ifstream in(tsv_path);
  if (!in.is_open())
    return false;
  std::string line;
  if (!std::getline(in, line))
    return false;
  return line.find("lens1_id") != std::string::npos && line.find("lens2_id") != std::string::npos;
}

// Legge un TSV con colonne x1, x2, x_focus, lens1_id, lens2_id e raggruppa
// per coppia di lenti. Lancia runtime_error se una colonna richiesta manca.
inline PerPairFocusMap load_per_pair_focus_map(const std::string& tsv_path) {
  std::ifstream in(tsv_path);
  if (!in.is_open())
    throw std::runtime_error("focus_map: impossibile aprire " + tsv_path);

  std::string line;
  if (!std::getline(in, line))
    throw std::runtime_error("focus_map: file vuoto: " + tsv_path);

  std::vector<std::string> headers;
  {
    std::istringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, '\t'))
      headers.push_back(tok);
  }
  auto find_col = [&](const std::string& name) -> int {
    for (int i = 0; i < static_cast<int>(headers.size()); ++i)
      if (headers[i] == name)
        return i;
    throw std::runtime_error("focus_map: colonna '" + name + "' non trovata in " + tsv_path);
  };
  int col_x1     = find_col("x1");
  int col_x2     = find_col("x2");
  int col_xfocus = find_col("x_focus");
  int col_l1     = find_col("lens1_id");
  int col_l2     = find_col("lens2_id");

  PerPairFocusMap m;
  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    std::vector<std::string> fields;
    {
      std::istringstream ss(line);
      std::string tok;
      while (std::getline(ss, tok, '\t'))
        fields.push_back(tok);
    }
    int ncols = static_cast<int>(fields.size());
    if (col_x1 >= ncols || col_x2 >= ncols || col_xfocus >= ncols || col_l1 >= ncols ||
        col_l2 >= ncols)
      continue;
    double x1     = std::stod(fields[col_x1]);
    double x2     = std::stod(fields[col_x2]);
    double xfocus = std::stod(fields[col_xfocus]);
    m[make_pair_key(fields[col_l1], fields[col_l2])][make_focus_key(x1, x2)] = xfocus;
  }
  return m;
}

} // namespace riptide

#endif // RIPTIDE_FOCUS_MAP_HPP
