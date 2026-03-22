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

#ifndef RIPTIDE_PSF_INTERPOLATOR_HPP
#define RIPTIDE_PSF_INTERPOLATOR_HPP

#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace riptide {

//  Strutture dati

// Un punto PSF: media e matrice di covarianza 2x2 per un dato raggio r
struct PSFPoint {
  double y_source;  // raggio (= y_source dalla simulazione) [mm]
  double mu_y;      // media y sul detector [mm]
  double mu_z;      // media z sul detector [mm]
  double cov_yy;    // varianza y [mm^2]
  double cov_yz;    // covarianza y-z [mm^2]
  double cov_zz;    // varianza z [mm^2]
  bool on_detector; // se true, la PSF è interamente sul fotocatodo
};

// Matrice di covarianza 2x2
struct Cov2 {
  double yy, yz, zz;
};

// Risultato dell'interpolazione per un singolo punto
struct PSFValue {
  double mu_y;
  double mu_z;
  Cov2 cov;
  bool on_detector;
};

// Un punto della traccia sul detector
struct TracePoint {
  double t;    // parametro lungo la traccia [mm]
  double r;    // distanza radiale dall'asse ottico [mm]
  double mu_y; // posizione media y sul detector [mm]
  double mu_z; // posizione media z sul detector [mm]
  Cov2 cov;    // matrice di covarianza associata [mm^2]
  bool valid;  // true se on_detector per questo punto
};

// Chiave di configurazione lenti (x1, x2)
struct LensConfig {
  double x1, x2;
  bool operator<(const LensConfig& o) const {
    const double eps = 1e-4;
    if (std::abs(x1 - o.x1) > eps)
      return x1 < o.x1;
    if (std::abs(x2 - o.x2) > eps)
      return x2 < o.x2;
    return false;
  }
};

// Database PSF: mappa (x1,x2) -> lista ordinata di PSFPoint per y_source crescente
using PSFDatabase = std::map<LensConfig, std::vector<PSFPoint>>;

/**
 * Risultato del fit ODR pesato della traccia media sul detector.
 * Modello: z = a*y + b
 */
struct LineFitResult {
  double a;
  double b;
  double sigma_a;
  double sigma_b;
  double cov_ab;
  double chi2;
  int ndof;
  double chi2_ndof;
  int n_iter;
  bool converged;

  // Residui perpendicolari
  std::vector<double> residuals;
  std::vector<double> residual_sig;
  std::vector<double> pull;

  // Quanti punti validi sono stati usati nel fit
  int n_points_used;
};

//  Funzioni pubbliche

/**
 * Carica psf_data.root in un PSFDatabase in memoria.
 * @param root_path  Path al file psf_data.root
 * @return           PSFDatabase popolato
 * @throws std::runtime_error se il file non esiste o i TTree mancano
 */
PSFDatabase load_psf_database(const std::string& root_path);

/**
 * Trova la configurazione (x1, x2) più vicina a quella richiesta nel database.
 * La distanza è euclidea nello spazio (x1, x2).
 * @param cfg  Configurazione richiesta
 * @param db   Database PSF
 * @return     Configurazione più vicina presente nel database
 */
LensConfig find_nearest_config(const LensConfig& cfg, const PSFDatabase& db);

/**
 * Interpola la PSF per un raggio r arbitrario e una configurazione (x1, x2).
 * Usa interpolazione lineare tra i due punti della griglia adiacenti.
 * @param r    Distanza radiale dall'asse ottico [mm]
 * @param cfg  Configurazione delle lenti
 * @param db   Database PSF caricato con load_psf_database()
 * @return     PSFValue con mu_y, mu_z e matrice di covarianza interpolati
 * @throws std::out_of_range se r è fuori dal range simulato o cfg non esiste
 */
PSFValue interpolate(double r, const LensConfig& cfg, const PSFDatabase& db);

/**
 * Costruisce la traccia media sul detector per una traccia ideale a distanza y0.
 * La traccia è parametrizzata da t in [-L/2, +L/2] con step dt.
 * Per ogni punto calcola r(t) = sqrt(y0^2 + t^2) e interpola la PSF.
 * @param y0   Distanza della traccia dall'asse ottico [mm]
 * @param cfg  Configurazione delle lenti
 * @param db   Database PSF
 * @param L    Lunghezza della traccia [mm], default 10.0
 * @param dt   Step di campionamento [mm], default 0.1
 * @return     Vettore di TracePoint ordinati per t crescente
 */
std::vector<TracePoint> build_trace(double y0, const LensConfig& cfg, const PSFDatabase& db,
                                    double L = 10.0, double dt = 0.1);

/**
 * Verifica se una traccia è valida in base alla frazione di punti on_detector.
 * @param trace              Traccia da verificare
 * @param point_valid_fraction  Frazione minima di punti on_detector
 * @return                   true se la traccia è valida, false altrimenti
 */
bool is_trace_valid(const std::vector<TracePoint>& trace, double point_valid_fraction = 0.75);

/**
 * Fit lineare pesato ODR della traccia media sul detector.
 *
 * Esegue il fit della retta z = a*y + b sui punti {(mu_y_i, mu_z_i)}
 * usando le matrici di covarianza Sigma_i come peso statistico, tramite
 * Orthogonal Distance Regression iterativa.
 *
 * Algoritmo:
 *   1. Stima iniziale (a, b) con fit non pesato (OLS su z vs y).
 *   2. Per ogni iterazione:
 *      a. Calcola n_hat = (-a, 1) / sqrt(1 + a^2)  (normale alla retta z=ay+b)
 *      b. Calcola sigma_{d,i}^2 = n_hat^T Sigma_i n_hat  per ogni punto
 *      c. Risolve il sistema lineare pesato 2x2 per (a, b):
 *         minimizza sum_i [ d_i^2 / sigma_{d,i}^2 ]
 *         con d_i = distanza perpendicolare (linearizzata attorno a n_hat corrente)
 *      d. Aggiorna (a, b); controlla convergenza su delta_a
 *   3. Calcola chi^2, ndof, residui e pull finali.
 *
 * Nota sulla linearizzazione: la distanza perpendicolare esatta e' non lineare
 * in (a, b). Per renderla lineare si usa la proiezione del residuo vettoriale
 * sul vettore normale corrente:
 *   d_i_lin = n_y * (mu_y_i - y_fit_i) + n_z * (mu_z_i - z_fit_i)
 *           = n_y * mu_y_i + n_z * mu_z_i - (n_y * y_hat_i + n_z * z_hat_i)
 * dove il punto sulla retta piu' vicino a (mu_y_i, mu_z_i) viene approssimato
 * con la proiezione ortogonale calcolata con la stima corrente di (a, b).
 * Questo e' l'approccio standard del IRLS per ODR (Boggs & Rogers 1990).
 *
 * @param trace     Traccia media prodotta da build_trace()
 * @param max_iter  Numero massimo di iterazioni (default: 20)
 * @param tol       Soglia di convergenza su |delta_a| (default: 1e-8)
 * @return          LineFitResult con tutti i parametri del fit
 * @throws std::invalid_argument se trace ha meno di 3 punti
 */
LineFitResult fit_trace(const std::vector<TracePoint>& trace, int max_iter = 20, double tol = 1e-8);

// Funzione di qualità Q(x1, x2)
// Parametri di configurazione per il calcolo di Q.
struct QConfig {
  std::vector<double> y0_values{};

  double y0_min = 0.0;
  double y0_max = 10.0;
  double dy0    = 0.1;

  double trace_L  = 10.0;
  double trace_dt = 0.1;

  int fit_max_iter = 20;
  double fit_tol   = 1e-8;

  // Soglie di validità 75%/75%
  double point_valid_fraction = 0.75; // frazione minima di punti on_detector per traccia valida
  double trace_valid_fraction = 0.75; // frazione minima di tracce valide per config valida

  bool verbose = false;
};

// Descrizione di una singola esclusione all'interno di compute_Q
struct QWarning {
  enum class Kind { FitNotConverged, BuildTraceFailed, FitFailed, TraceInvalid, ConfigInvalid };
  Kind kind;
  double y0;
  double a_final;
  std::string message;
};

// Risultato del calcolo di Q per una singola configurazione lenti.
struct QResult {
  double Q;
  int n_traces;      // tracce valide usate nella somma
  int n_failed;      // tracce scartate per errore nel fit
  int n_invalid;     // tracce scartate per validità PSF insufficiente
  bool config_valid; // false se la config non supera la soglia trace_valid_fraction

  std::vector<double> y0_used;
  std::vector<double> chi2_per_y0;
  std::vector<double> chi2_ndof_per_y0;
  std::vector<bool> trace_valid_flags; // una entry per ogni y0 campionato

  std::vector<QWarning> warnings;
};

/**
 * Calcola la funzione di qualità Q(x1, x2) per una configurazione di lenti.
 *
 * Per ogni valore di y0 in QConfig::y0_values (o generato dall'intervallo):
 *   1. Costruisce la traccia media build_trace(y0, cfg, db, L, dt)
 *   2. Esegue il fit ODR pesato fit_trace(trace)
 *   3. Accumula chi^2 → Q = sum_i chi^2(y0_i)
 *
 * Le tracce con fit non convergente vengono conteggiate in QResult::n_failed
 * ma NON contribuiscono a Q (comportamento configurabile tramite
 * include_non_converged — se true, vengono sommate comunque).
 *
 * @param cfg                  Configurazione lenti (deve essere presente nel db)
 * @param db                   Database PSF caricato con load_psf_database()
 * @param qcfg                 Parametri di campionamento e fit
 * @param include_non_converged Se true, include chi^2 anche di fit non convergenti
 * @return                     QResult con il valore di Q e i contributi per-traccia
 * @throws std::invalid_argument se db non contiene cfg
 */
QResult compute_Q(const LensConfig& cfg, const PSFDatabase& db, const QConfig& qcfg = QConfig{},
                  bool include_non_converged = false);

} // namespace riptide

#endif // RIPTIDE_PSF_INTERPOLATOR_HPP
