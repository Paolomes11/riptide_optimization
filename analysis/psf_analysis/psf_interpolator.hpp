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
  double y_source; // raggio (= y_source dalla simulazione) [mm]
  double mu_y;     // media y sul detector [mm]
  double mu_z;     // media z sul detector [mm]
  double cov_yy;   // varianza y [mm^2]
  double cov_yz;   // covarianza y-z [mm^2]
  double cov_zz;   // varianza z [mm^2]
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
};

// Un punto della traccia sul detector
struct TracePoint {
  double t;    // parametro lungo la traccia [mm]
  double r;    // distanza radiale dall'asse ottico [mm]
  double mu_y; // posizione media y sul detector [mm]
  double mu_z; // posizione media z sul detector [mm]
  Cov2 cov;    // matrice di covarianza associata [mm^2]
};

// Chiave di configurazione lenti (x1, x2) con confronto per std::map
// Usa una tolleranza di 1e-4 mm per evitare problemi di floating point
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
 *
 * Modello: z = a*y + b
 *
 * Il fit minimizza il chi^2 con distanza perpendicolare pesata dalla matrice
 * di covarianza di ciascun punto (Orthogonal Distance Regression con Sigma_i
 * completa). Il problema e' non lineare perche' il vettore normale n_hat
 * dipende da 'a'; si risolve per iterazione.
 */
struct LineFitResult {
  // Parametri della retta
  double a; // pendenza [adimensionale]
  double b; // intercetta [mm]

  // Incertezza sui parametri
  double sigma_a; // incertezza su a
  double sigma_b; // incertezza su b
  double cov_ab;  // covarianza tra a e b

  // Chi^2 del fit
  double chi2;
  int ndof;         // gradi di libertà
  double chi2_ndof; // chi^2 ridotto

  // Valori per iterazione del fit
  int n_iter; // numero di iterazioni fatte
  bool converged;   // se è convergente

  // Residui perpendicolari pesati per ciascun punto
  // d_i = distanza perpendicolare del punto i dalla retta stimata [mm]
  // sigma_d_i = sqrt(n_hat^T Sigma_i n_hat) [mm]
  std::vector<double> residuals;    // d_i
  std::vector<double> residual_sig; // sigma_{d,i}
  std::vector<double> pull;         // d_i / sigma_{d,i}
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

} // namespace riptide

#endif // RIPTIDE_PSF_INTERPOLATOR_HPP
