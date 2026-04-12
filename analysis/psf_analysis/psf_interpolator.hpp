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

// Un punto PSF: media e matrice di covarianza 2x2 per un dato (x_source, y_source)
struct PSFPoint {
  double x_source;  // x_source from simulation [mm]
  double y_source;  // y_source from simulation [mm] (equivale a raggio r nella vecchia versione)
  double mu_y;      // media y sul detector [mm]
  double mu_z;      // media z sul detector [mm]
  double cov_yy;    // varianza y [mm^2]
  double cov_yz;    // covarianza y-z [mm^2]
  double cov_zz;    // varianza z [mm^2]
  bool on_detector; // mantenuto per compatibilità (true se n_hits >= soglia)
  double n_hits;    // numero di fotoni arrivati (pesato, usato per efficienza locale)
  double n_hits_count;
};

// Matrice di covarianza 2x2
struct Cov2 {
  double yy, yz, zz;
};

// Risultato dell'interpolazione per un singolo punto (x, y)
struct PSFValue {
  double mu_y;
  double mu_z;
  Cov2 cov;
  bool on_detector;
  double n_hits_interp;
  double n_hits_count_interp;
};

// Un punto della traccia sul detector
struct TracePoint {
  double t;     // parametro lungo la traccia [mm]
  double r;     // mantenuto per compatibilità (distanza radiale x=0)
  double x_src; // posizione x sorgente associata [mm]
  double y_src; // posizione y sorgente associata [mm]
  double z_src; // posizione z sorgente associata [mm]
  double mu_y;  // posizione media y sul detector [mm]
  double mu_z;  // posizione media z sul detector [mm]
  // Errore standard sulla posizione media (mu_y, mu_z) [mm^2].
  // = Sigma_distribuzione / n_hits_count.
  // NON e' la covarianza della distribuzione dei fotoni sulla macchia PSF.
  Cov2 cov;
  bool valid; // mantenuto per compatibilità
  double n_hits;
  double n_hits_count;
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

/// Asse scelto automaticamente da fit_trace() in base allo spread dei dati.
/// ZvsY → modello z = a·y + b  (default, tracce oblique o quasi-orizzontali)
/// YvsZ → modello y = a·z + b  (attivato quando spread_z > 3·spread_y)
enum class FitAxis { ZvsY, YvsZ };

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
  FitAxis axis = FitAxis::ZvsY; // asse selezionato automaticamente

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
// Interpolazione 2D (x_src, y_src)
PSFValue interpolate(double x, double y, const LensConfig& cfg, const PSFDatabase& db);

struct Point3D {
  double x, y, z;
};

/**
 * Costruisce la traccia media sul detector campionando i punti lungo il segmento P1-P2.
 * Per ogni punto, sfrutta la simmetria circolare ruotando i risultati della PSF(x, r).
 */
std::vector<TracePoint> build_trace_3d(const Point3D& p1, const Point3D& p2, const LensConfig& cfg,
                                       const PSFDatabase& db, double dt);

/**
 * Wrapper di compatibilità per costruire una traccia rettilinea lungo l'asse X.
 */
std::vector<TracePoint> build_trace(double y0, const LensConfig& cfg, const PSFDatabase& db,
                                    double L = 10.0, double dt = 0.1);

/**
 * Verifica se una traccia è valida in base alla frazione di punti validi.
 * @param trace                 Traccia da verificare
 * @param point_valid_fraction  Frazione minima di punti validi (pt.valid)
 * @return                      true se la traccia è valida, false altrimenti
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
LineFitResult fit_trace(const std::vector<TracePoint>& trace, double min_hits_per_point = 10.0,
                        int max_iter = 20, double tol = 1e-8);

// Funzione di qualità Q(x1, x2)
// Parametri di configurazione per il calcolo di Q tramite tracce casuali 3D.
struct QConfig {
  // Dimensioni scintillatore [mm] (x=lunghezza, y,z=sezione)
  // Default: 60x20x20 (x in [-30, 30], y,z in [-10, 10])
  double scint_x = 60.0;
  double scint_y = 20.0;
  double scint_z = 20.0;

  int n_tracks = 100; // Numero di tracce casuali per config

  double trace_dt = 0.1; // Passo campionamento lungo la traccia [mm]

  int fit_max_iter = 20;
  double fit_tol   = 1e-8;

  // Soglie di validità
  double min_hits_per_point   = 10.0; // hit minime per considerare un punto PSF valido
  double trace_valid_fraction = 0.75; // frazione minima di tracce valide per config valida

  // Temporal unfolding
  // Se > 0: offset fisso in mm applicato per passo della traccia.
  // Se == 0 (default): calcolato automaticamente come L_traccia / (N-1),
  //   in modo che l'offset totale sia pari alla lunghezza della traccia.
  double z_unfold_step          = 0.0;  // [mm/passo], 0 = automatico
  bool apply_temporal_unfolding = true; // default attivo

  bool verbose = false;
};

// Descrizione di una singola esclusione all'interno di compute_Q
struct QWarning {
  enum class Kind { FitNotConverged, BuildTraceFailed, FitFailed, TraceInvalid, ConfigInvalid };
  Kind kind;
  int track_id;
  double a_final;
  std::string message;
};

// Risultato del calcolo di Q per una singola configurazione lenti.
struct QResult {
  double Q;          // media del Chi-squared ridotto sulle tracce valide
  int n_traces;      // tracce valide usate nella media
  int n_failed;      // tracce scartate per errore nel fit o Chi2 non valido
  int n_invalid;     // tracce scartate per validità PSF insufficiente
  bool config_valid; // false se la config non supera la soglia trace_valid_fraction

  std::vector<double> chi2_per_trace;
  std::vector<double> chi2_ndof_per_trace;
  std::vector<bool> trace_valid_flags; // una entry per ogni traccia campionata

  std::vector<QWarning> warnings;
};

/**
 * Calcola la funzione di qualità Q(x1, x2) per una configurazione di lenti.
 *
 * Genera n_tracks casuali nello scintillatore, ne calcola la proiezione 3D
 * sul detector tramite simmetria circolare e ne esegue il fit lineare.
 */
QResult compute_Q(const LensConfig& cfg, const PSFDatabase& db, const QConfig& qcfg = QConfig{},
                  bool include_non_converged = false);

/**
 * Risultato del calcolo di copertura per una singola configurazione lenti.
 *
 * "Copertura" = frazione media di punti on_detector su tutte le tracce campionate
 * e su tutti i punti di ogni traccia.  Misura quanta luce arriva effettivamente
 * sul fotocatodo, indipendentemente dalla qualità del fit ODR.
 *
 * Definizione formale:
 *
 *   coverage(x1,x2) = (1/N_y0) * Σ_{y0} [ n_valid(y0) / N_pts_per_trace ]
 *
 * dove:
 *   n_valid(y0) = numero di TracePoint con valid==true per la traccia a y0
 *   N_pts_per_trace = numero totale di punti per traccia (= round(L/dt)+1)
 *   N_y0 = numero di tracce per cui build_trace ha avuto successo
 *
 * Un valore di 1.0 significa che tutti i fotoni simulati raggiungono il
 * fotocatodo per ogni posizione della sorgente nell'intervallo campionato.
 * Un valore di 0.0 significa copertura nulla (lenti totalmente fuori fuoco
 * o configurazione geometricamente impossibile).
 */
struct CoverageResult {
  double coverage;    // frazione media on_detector ∈ [0, 1]
  int n_y0_evaluated; // y0 per cui build_trace ha avuto successo
  int n_y0_requested; // y0 totali richiesti dalla configurazione
  bool config_valid;  // false se n_y0_evaluated == 0
};

/**
 * Calcola la mappa di copertura geometrica per una configurazione di lenti.
 *
 * Per ogni valore di y0 nel campionamento (stessa griglia di compute_Q):
 *   1. Costruisce la traccia con build_trace(y0, cfg, db, L, dt)
 *   2. Conta i TracePoint con valid == true
 *   3. Accumula coverage += n_valid / N_total
 * Restituisce la media su tutti i y0 campionati con successo.
 *
 * Non esegue nessun fit ODR: il costo computazionale è trascurabile.
 * Ideale come primo passo diagnostico prima di eseguire q_map.
 *
 * @param cfg   Configurazione lenti (deve essere presente nel db)
 * @param db    Database PSF
 * @param qcfg  Parametri di campionamento (usa y0_min/max/dy0, trace_L, trace_dt)
 * @return      CoverageResult
 * @throws std::invalid_argument se db non contiene cfg
 */
CoverageResult compute_coverage(const LensConfig& cfg, const PSFDatabase& db,
                                const QConfig& qcfg = QConfig{});

} // namespace riptide

#endif // RIPTIDE_PSF_INTERPOLATOR_HPP
