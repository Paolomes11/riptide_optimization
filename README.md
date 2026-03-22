# riptide_optimization

Simulazione Monte Carlo Geant4 per l'ottimizzazione del posizionamento di un sistema ottico a doppia lente (UVFS) davanti a un fotocatodo GaAsP. Il progetto fa parte dell'esperimento **RIPTIDE** e comprende due programmi principali — `optimization` e `lens_simulation` — cinque strumenti di analisi ROOT, e una libreria di analisi PSF (`psf_analysis`) con test unitari integrati.

---

## Indice

1. [Requisiti](#requisiti)
2. [Struttura del progetto](#struttura-del-progetto)
3. [Setup fisico](#setup-fisico)
4. [Build](#build)
5. [Programma: optimization](#programma-optimization)
6. [Programma: lens\_simulation](#programma-lens_simulation)
7. [Output su SSD esterna](#output-su-ssd-esterna)
8. [Parallelizzazione a processi](#parallelizzazione-a-processi)
9. [Strumenti di analisi](#strumenti-di-analisi)
10. [Libreria psf\_analysis](#libreria-psf_analysis)
11. [Test unitari](#test-unitari)
12. [Configurazione raccomandata per 'q_map'](#Configurazione-raccomandata-per-`q_map`)
13. [File di configurazione](#file-di-configurazione)
14. [Formato dei file ROOT](#formato-dei-file-root)
15. [Workflow completo](#workflow-completo)

---

## Requisiti

| Dipendenza | Versione minima | Note |
|---|---|---|
| CMake | 3.26 | |
| C++ | 17 | |
| Geant4 | 11.x | Con supporto `ui_all`, `vis_all`, GDML, ottica |
| ROOT | 6.x | Componenti: `Core`, `Hist`, `RIO`, `Graf` |
| spdlog | qualsiasi | Logging strutturato |
| nlohmann/json | (header-only, incluso) | Lettura configurazione |
| lyra | (header-only, incluso in `external/`) | Parsing CLI |

---

## Struttura del progetto

```
riptide_optimization/
├── scripts/
│   ├── run.sh                           # Lancia lens_simulation o optimization (locale o SSD)
│   └── monitor.sh                       # Monitoraggio progresso chunk paralleli
│
├── programs/
│   ├── optimization_main.cpp            # Main per la scansione di efficienza geometrica
│   └── lens_simulation_main.cpp         # Main per la simulazione beam scan PSF
│
├── src/
│   ├── common/                          # Codice condiviso
│   │   ├── physics_list.cpp
│   │   └── primary_generator_action.cpp
│   ├── optimization/                    # Sorgenti programma optimization
│   │   ├── detector_construction.cpp
│   │   ├── action_initialization.cpp
│   │   ├── event_action.cpp
│   │   ├── run_action.cpp
│   │   ├── sensitive_detector.cpp
│   │   └── optimizer.cpp
│   └── lens_simulation/                 # Sorgenti programma lens_simulation
│       ├── detector_construction.cpp
│       ├── action_initialization.cpp
│       ├── event_action.cpp
│       ├── run_action.cpp
│       ├── sensitive_detector.cpp
│       └── lens_scan.cpp
│
├── include/                             # Header files (specchi di src/)
│   ├── common/
│   ├── optimization/
│   └── lens_simulation/
│
├── analysis/
│   ├── plot2D.cpp                       # Mappa 2D efficienza geometrica
│   ├── beam_scan_plot.cpp               # Posizione fotoni vs posizione sorgente
│   ├── m_c_creator.cpp                  # Istogramma 2D hit su detector per un run
│   ├── psf_extractor.cpp                # Estrazione media e covarianza PSF per tutti i run
│   ├── CMakeLists.txt
│   └── psf_analysis/                    # Libreria PSF + analisi traccia + ottimizzazione
│       ├── psf_interpolator.hpp         # API pubblica: strutture dati, dichiarazioni
│       ├── psf_interpolator.cpp         # load, interpolate, build_trace, fit_trace, compute_Q
│       ├── trace_viewer.cpp             # Visualizza traccia media con ellissi di covarianza
│       ├── q_map.cpp                    # ★ Mappa 2D di Q(x1,x2) — ottimizzazione PSF-based
│       ├── test_fit_trace.cpp           # Test unitari per fit_trace() e compute_Q()
│       └── CMakeLists.txt
│
├── geometry/
│   ├── main.gdml
│   ├── define.xml
│   ├── materials.xml
│   ├── solids.xml
│   └── structure.xml
│
├── macros/
│   ├── optimization.mac
│   ├── lens_simulation.mac
│   ├── run.mac
│   └── vis.mac
│
├── config/
│   └── config.json
│
├── external/
│   └── lyra/
│
└── output/
    ├── events.root
    ├── efficiency2D.png
    ├── lens_simulation/
    │   └── lens.root
    ├── mean_covariance_maps/
    ├── psf/
    │   └── psf_data.root
    └── psf_analysis/
        ├── q_map.png                    # ★ output di q_map
        └── trace_*.png
```

---

## Setup fisico

Il sistema simulato è composto da tre elementi posizionati lungo l'asse X, all'interno di un volume mondo cubico di 1000×1000×1000 mm³ riempito d'aria:

```
Sorgente fotoni  →  [Lente 75mm]  →  [Lente 60mm]  →  [Fotocatodo GaAsP 16×16mm]
      (GPS)           (UVFS)           (UVFS)              (sensore)
         x=0        x ≈ 14–170mm    x ≈ 45–186mm           x ≈ 180mm
```

**Lente 75 mm** (`lens75`): ellissoide in UV Fused Silica (UVFS), raggio 38.6 mm, spessore 12.5 mm.

**Lente 60 mm** (`lens60`): ellissoide UVFS, raggio 30.9 mm, spessore 16.3 mm.

**Fotocatodo GaAsP**: lastra quadrata 16×16×0.01 mm, indice di rifrazione 3.5–3.8 nel range 2–4 eV, lunghezza di assorbimento ~1 µm.

**Sorgente**: fotoni ottici a 2.5 eV (≈ 496 nm), generati da GPS Geant4 su un disco con distribuzione angolare isotropa diretta lungo +X.

---

## Build

```bash
# Configurazione (Debug o Release)
cmake -S . -B build/ -G "Ninja Multi-Config"

# Compilazione
cmake --build build/ --config Release

# Eseguibili prodotti:
#   build/Release/optimization_main
#   build/Release/lens_simulation_main
#   build/analysis/Release/plot2D
#   build/analysis/Release/beam_scan_plot
#   build/analysis/Release/m_c_creator
#   build/analysis/Release/psf_extractor
#   build/analysis/psf_analysis/Release/trace_viewer
#   build/analysis/psf_analysis/Release/q_map           ★ nuovo
#   build/analysis/psf_analysis/Release/test_fit_trace
```

> **Nota**: CMake crea automaticamente le cartelle `output/`, `output/lens_simulation/`, `output/mean_covariance_maps/`, `output/psf/` e `output/psf_analysis/` durante la configurazione.

---

## Programma: optimization

**Eseguibile**: `optimization_main`

**Scopo**: scansione grezza dell'efficienza geometrica del sistema a doppia lente. Per ogni coppia di posizioni (x1, x2) delle lenti, esegue una simulazione con 10000 fotoni emessi da una sorgente circolare a r=10 mm e conta quanti raggiungono il fotocatodo.

### Utilizzo

```bash
./build/Release/optimization_main -g geometry/main.gdml -o
./build/Release/optimization_main -g geometry/main.gdml -v
./build/Release/optimization_main --help
```

### Opzioni CLI

| Flag | Tipo | Descrizione |
|---|---|---|
| `-g`, `--geometry` | path | **Obbligatorio.** Percorso al file GDML |
| `-m`, `--macro` | path | Macro Geant4 (default: `macros/optimization.mac`) |
| `-v`, `--visualize` | flag | Visualizzazione interattiva OpenGL/Qt |
| `-b`, `--batch` | flag | Modalità batch senza UI |
| `-o`, `--optimize` | flag | Avvia la scansione di ottimizzazione |
| `--output` | path | File ROOT di output (default: `output/events.root`) |
| `--config` | path | File `config.json` (default: `config/config.json`) |
| `--ssd` | flag | Output sull'SSD esterna con timestamp automatico |
| `--ssd-mount` | path | Mount point dell'SSD (default: `/mnt/external_ssd`) |

### Output

File ROOT: `output/events.root` — TTree `events` con una riga per ogni fotone rilevato:

| Branch | Tipo | Descrizione |
|---|---|---|
| `x1` | `Double_t` | Posizione della lente 75mm [mm] |
| `x2` | `Double_t` | Posizione della lente 60mm [mm] |
| `config_id` | `Int_t` | Indice della configurazione |

---

## Programma: lens\_simulation

**Eseguibile**: `lens_simulation_main`

**Scopo**: simulazione dettagliata del beam scan per la caratterizzazione della PSF. Per ogni configurazione di lenti e per ogni posizione della sorgente, registra la posizione (y, z) di ogni fotone sul fotocatodo.

### Utilizzo

```bash
./build/Release/lens_simulation_main -g geometry/main.gdml -l
./build/Release/lens_simulation_main -g geometry/main.gdml -v
./build/Release/lens_simulation_main --help
```

### Opzioni CLI

| Flag | Tipo | Descrizione |
|---|---|---|
| `-g`, `--geometry` | path | **Obbligatorio.** Percorso al file GDML |
| `-m`, `--macro` | path | Macro Geant4 (default: `macros/lens_simulation.mac`) |
| `-v`, `--visualize` | flag | Visualizzazione interattiva |
| `-b`, `--batch` | flag | Modalità batch |
| `-l`, `--lens-sim` | flag | Avvia il beam scan completo |
| `--output` | path | File ROOT di output (default: `output/lens_simulation/lens.root`) |
| `--config` | path | File `config.json` (default: `config/config.json`) |
| `--ssd` | flag | Output sull'SSD esterna |
| `--ssd-mount` | path | Mount point dell'SSD |

### Output

File ROOT: `output/lens_simulation/lens.root` — tre TTree:

**`Configurations`**: una riga per coppia di posizioni lenti.

| Branch | Tipo | Descrizione |
|---|---|---|
| `config_id` | `Int_t` | Indice configurazione |
| `x1` | `Double_t` | Posizione lente 75mm [mm] |
| `x2` | `Double_t` | Posizione lente 60mm [mm] |

**`Runs`**: una riga per ogni esecuzione (configurazione × posizione sorgente).

| Branch | Tipo | Descrizione |
|---|---|---|
| `run_id` | `Int_t` | Indice run globale |
| `config_id` | `Int_t` | Configurazione di appartenenza |
| `x_source` | `Float_t` | Posizione Y della sorgente [mm] |
| `n_hits` | `Int_t` | Fotoni rilevati in questo run |

**`Hits`**: una riga per ogni fotone rilevato.

| Branch | Tipo | Descrizione |
|---|---|---|
| `y_hit` | `Float_t` | Coordinata Y sul fotocatodo [mm] |
| `z_hit` | `Float_t` | Coordinata Z sul fotocatodo [mm] |

> **Mapping hits → run**: le hit sono scritte sequenzialmente. La riga i-esima di `Runs` corrisponde alle `n_hits[i]` righe di `Hits` a partire dall'offset cumulativo `sum(n_hits[0..i-1])`.

---

## Output su SSD esterna

### Trovare e montare il device

```bash
lsblk
sudo mkdir -p /mnt/external_ssd
sudo mount -o noatime,nodiratime,discard /dev/nvme1n1p1 /mnt/external_ssd
mountpoint -q /mnt/external_ssd && echo "OK"
```

### Utilizzo con `--ssd`

```bash
./build/Release/lens_simulation_main -g geometry/main.gdml -b -l --ssd
./build/Release/optimization_main    -g geometry/main.gdml -b -o --ssd
./build/Release/lens_simulation_main -g geometry/main.gdml -b -l \
    --ssd --ssd-mount /mnt/myusb
```

Il path di output viene generato automaticamente con timestamp:
```
/mnt/external_ssd/riptide/runs/run_20260315_094512/lens.root
```

---

## Parallelizzazione a processi

`lens_simulation` supporta la parallelizzazione tramite processi indipendenti. Lo spazio delle configurazioni viene diviso in N chunk, ognuno con il proprio file ROOT. Al termine i chunk vengono uniti con `hadd`.

```bash
./scripts/run.sh lens local --jobs 4
./scripts/run.sh lens ssd   --jobs 8
./scripts/run.sh lens ssd   --jobs $(nproc --all)
```

I `config_id` e `run_id` nel file finale sono globalmente unici e contigui grazie agli offset calcolati dallo script.

---

## Strumenti di analisi

### plot2D

Genera la mappa 2D di efficienza geometrica da `output/events.root`.

```bash
./build/analysis/Release/plot2D
# → output/efficiency2D.png
```

### beam\_scan\_plot

Per una coppia (x1, x2) fissata, calcola posizione media e deviazione standard dei fotoni al variare della posizione sorgente. Applica filtro outlier 2σ iterativo.

```bash
./build/analysis/Release/beam_scan_plot <x1> <x2> [file.root]
./build/analysis/Release/beam_scan_plot 94.9 186.4
# → output/lens_simulation/beam_x1_94.90_x2_186.40.png
```

### m\_c\_creator

Istogramma 2D delle hit per una terna (x1, x2, y0) fissata.

```bash
./build/analysis/Release/m_c_creator <x1> <x2> <y0>
./build/analysis/Release/m_c_creator 94.9 186.4 5.0
# → output/mean_covariance_maps/detector_hits_config_<id>_y0_5.0.png
```

### psf\_extractor

Estrae media bidimensionale e matrice di covarianza della PSF per tutti i run. Applica filtro outlier ellittico a 2σ (4 iterazioni).

```bash
./build/analysis/Release/psf_extractor [input.root] [output.root]
./build/analysis/Release/psf_extractor \
    output/lens_simulation/lens.root \
    output/psf/psf_data.root
```

Output `psf_data.root` — TTree `PSF`:

| Branch | Tipo | Descrizione |
|---|---|---|
| `config_id` | `Int_t` | Indice configurazione |
| `x1`, `x2` | `Double_t` | Posizioni lenti [mm] |
| `y_source` | `Float_t` | Posizione sorgente [mm] |
| `mean_y`, `mean_z` | `Double_t` | Media PSF [mm] |
| `cov_yy`, `cov_yz`, `cov_zz` | `Double_t` | Matrice di covarianza [mm²] |
| `n_hits_filtered` | `Int_t` | Hit dopo filtro 2σ |
| `n_hits_raw` | `Int_t` | Hit prima del filtro |

---

## Libreria psf\_analysis

La libreria `psf_analysis` (in `analysis/psf_analysis/`) implementa la catena analitica PSF-based. È compilata come libreria statica e usata da `trace_viewer`, `q_map` e `test_fit_trace`.

### API pubblica (`psf_interpolator.hpp`)

#### Strutture dati principali

| Struttura | Descrizione |
|---|---|
| `PSFPoint` | Un punto del database PSF: `y_source`, `mu_y`, `mu_z`, covarianza |
| `Cov2` | Matrice di covarianza 2×2: `yy`, `yz`, `zz` |
| `PSFValue` | Risultato interpolato: `mu_y`, `mu_z`, `cov` |
| `TracePoint` | Punto della traccia sul detector: `t`, `r`, `mu_y`, `mu_z`, `cov` |
| `LensConfig` | Chiave di configurazione `(x1, x2)` con confronto a tolleranza 1e-4 mm |
| `LineFitResult` | Risultato completo del fit ODR |
| `QConfig` | Parametri di campionamento per `compute_Q` |
| `QResult` | Risultato di `compute_Q`: valore di Q, contributi per-traccia |

#### Funzioni

**`load_psf_database(path)`** — carica `psf_data.root` in memoria come `PSFDatabase`.

**`find_nearest_config(cfg, db)`** — trova la configurazione più vicina nel database (distanza euclidea in `(x1, x2)`).

**`interpolate(r, cfg, db)`** — interpola `mu_y`, `mu_z` e `Σ` per un raggio `r` arbitrario. Usa interpolazione lineare; clamp agli estremi senza estrapolazione.

**`build_trace(y0, cfg, db, L=10, dt=0.1)`** — costruisce la traccia media sul detector per una traccia ideale a distanza `y0` dall'asse ottico, parametrizzata da `t ∈ [-L/2, +L/2]` con step `dt`.

**`fit_trace(trace, max_iter=20, tol=1e-8)`** — fit lineare pesato ODR della traccia.

**`compute_Q(cfg, db, qcfg, include_non_converged=false)`** — calcola la funzione di qualità Q(x1,x2).

#### `fit_trace` — Fit ODR iterativo

Esegue il fit della retta `z = a·y + b` sui punti `{(mu_y_i, mu_z_i)}` usando le matrici di covarianza `Σ_i` come peso statistico, tramite Orthogonal Distance Regression (ODR) con schema IRLS.

**Algoritmo:**

1. Stima iniziale `(a, b)` con OLS non pesato.
2. Per ogni iterazione: calcola `n̂`, poi i pesi `w_i = 1/(n̂ᵀ Σ_i n̂)`, risolve il WLS 2×2 in forma chiusa, controlla convergenza su `|Δa| < tol`.
3. Calcola χ², residui e pull finali con distanza perpendicolare esatta.

#### `compute_Q` — Funzione di qualità

Implementa la funzione di qualità definita nella tesi (Sez. 8):

```
Q(x1, x2) = Σ_{y0} chi²(y0, x1, x2)
```

Per ogni valore di `y0` nel campionamento: costruisce la traccia con `build_trace`, esegue `fit_trace`, accumula χ².

**`QConfig` — campi:**

| Campo | Default | Descrizione |
|---|---|---|
| `y0_values` | `{}` | Lista esplicita di y0; se vuota, usa l'intervallo sotto |
| `y0_min` | `0.0` | Inizio intervallo y0 [mm] |
| `y0_max` | `10.0` | Fine intervallo y0 [mm] |
| `dy0` | `0.1` | Passo di campionamento [mm] |
| `trace_L` | `10.0` | Lunghezza traccia [mm] |
| `trace_dt` | `0.1` | Step traccia [mm] |
| `fit_max_iter` | `20` | Iterazioni massime IRLS |
| `fit_tol` | `1e-8` | Soglia convergenza |

**`QResult` — campi:**

| Campo | Descrizione |
|---|---|
| `Q` | Valore totale della funzione di qualità |
| `n_traces` | Numero di tracce usate nella somma |
| `n_failed` | Tracce scartate (fit non convergente o dati mancanti) |
| `y0_used` | Valori di y0 effettivamente usati |
| `chi2_per_y0` | χ² del fit ODR per ogni y0 |
| `chi2_ndof_per_y0` | χ²/ndof per ogni y0 |

**Nota di design:** `compute_Q` richiede che la configurazione sia *esatta* nel database. Usa `find_nearest_config` prima se operi su configurazioni non campionate.

**Utilizzo tipico** in un loop di ottimizzazione:

```cpp
auto db = riptide::load_psf_database("output/psf/psf_data.root");

riptide::QConfig qcfg;
qcfg.y0_min = 0.0;  qcfg.y0_max = 10.0;  qcfg.dy0 = 0.1;

double Q_best = std::numeric_limits<double>::max();
riptide::LensConfig best_cfg{};

for (const auto& [cfg, _] : db) {
    auto res = riptide::compute_Q(cfg, db, qcfg);
    if (res.n_traces > 0 && res.Q < Q_best) {
        Q_best   = res.Q;
        best_cfg = cfg;
    }
}
std::cout << "Ottimale: x1=" << best_cfg.x1
          << " x2=" << best_cfg.x2
          << "  Q=" << Q_best << "\n";
```

### trace\_viewer

Visualizza la traccia media sul detector con ellissi di covarianza, colorate per parametro `t` lungo la traccia (palette Rainbow: blu = inizio, rosso = fine).

```bash
./build/analysis/psf_analysis/Release/trace_viewer \
    --x1 94.9 --x2 186.4 --y0 5.0

./build/analysis/psf_analysis/Release/trace_viewer \
    --x1 94.9 --x2 186.4 --y0 5.0    \
    --psf    output/psf/psf_data.root  \
    --output output/psf_analysis/mia_traccia.png \
    --dt 0.1   \
    --L  10.0  \
    --sigma 1.0
# → output/psf_analysis/trace_x1_94.9_x2_186.4_y0_5.0.png
```

| Opzione | Default | Descrizione |
|---|---|---|
| `--x1`, `--x2` | — | **Obbligatori.** Posizioni lenti [mm] |
| `--y0` | — | **Obbligatorio.** Distanza radiale sorgente [mm] |
| `--psf` | `output/psf/psf_data.root` | Database PSF |
| `--output` | auto | Path immagine PNG |
| `--dt` | `0.1` | Step traccia [mm] |
| `--L` | `10.0` | Lunghezza traccia [mm] |
| `--sigma` | `1.0` | Scala delle ellissi di covarianza (in unità σ) |

### q\_map ★

Genera la mappa 2D della funzione di qualità `Q(x1, x2)` su tutte le configurazioni presenti nel database PSF. Individua e marca automaticamente il minimo globale di Q.

```bash
# Uso minimo (tutti i parametri dai default)
./build/analysis/psf_analysis/Release/q_map

# Uso completo
./build/analysis/psf_analysis/Release/q_map \
    --psf    output/psf/psf_data.root        \
    --config config/config.json              \
    --output output/psf_analysis/q_map.png   \
    --y0-min 0.0  --y0-max 10.0  --dy0 0.1  \
    --L 10.0  --dt 0.1                       \
    --log                                    \
    --norm                                   \
    --tsv output/psf_analysis/q_map.tsv
# → output/psf_analysis/q_map.png
# → output/psf_analysis/q_map.tsv  (se --tsv specificato)
```

| Opzione | Default | Descrizione |
|---|---|---|
| `--psf` | `output/psf/psf_data.root` | Database PSF |
| `--config` | `config/config.json` | Parametri griglia lenti |
| `--output` | `output/psf_analysis/q_map.png` | Immagine PNG di output |
| `--tsv` | (disabilitato) | Esporta valori Q in formato TSV |
| `--y0-min` | `0.0` | Minimo y0 per il campionamento tracce [mm] |
| `--y0-max` | `10.0` | Massimo y0 [mm] |
| `--dy0` | `0.1` | Passo di campionamento y0 [mm] |
| `--L` | `10.0` | Lunghezza traccia ideale [mm] |
| `--dt` | `0.1` | Step di campionamento della traccia [mm] |
| `--log` | off | Scala logaritmica sull'asse Z (colori) |
| `--norm` | off | Normalizza Q per numero di tracce (mostra χ² medio) |

**Output grafico:** mappa TH2D con palette `kBird`, scala colori troncata ai percentili 5–95 per massimizzare il contrasto nell'area di interesse. Marker rosso a stella sulla configurazione ottimale. Info panel inferiore con parametri di campionamento e coordinate del minimo.

**Output TSV** (con `--tsv`): tabella `x1 \t x2 \t Q \t n_traces \t n_failed` — compatibile con qualsiasi tool di analisi esterno (Python, gnuplot, Excel).

---

## Test unitari

### Esecuzione

```bash
# Esecuzione diretta
./build/analysis/psf_analysis/Release/test_fit_trace

# Tramite CTest
cd build && ctest -R fit_trace_unit -V
```

L'eseguibile ritorna `0` se tutti i test passano, `1` altrimenti.

### Note sull'inizializzazione dei dati sintetici

I test costruiscono `TracePoint` e `PSFPoint` sintetici in memoria, bypassando i file ROOT. Due invarianti sono **essenziali**:

- **`TracePoint::valid = true`** — `fit_trace()` filtra i punti con `valid == false` (li considera fuori dal fotocatodo). Se non impostato esplicitamente, rimane `false` per zero-inizializzazione e `fit_trace()` lancia `std::invalid_argument` per assenza di punti utilizzabili.
- **`PSFPoint::on_detector = true`** — `build_trace()` copia questo flag in `TracePoint::valid`. Se `false`, `is_trace_valid()` ritorna `false` e `compute_Q()` non accumula nessun χ².

### Casi di test

I test costruiscono dati sintetici con soluzione analitica nota — nessun file ROOT necessario.

#### Suite T — `fit_trace`

| Test | Scenario | Verifica |
|---|---|---|
| **T1** | Retta `z = y`, cov isotropa uniforme | `a=1`, `b=0`, `χ²≈0`, pull tutti ≈ 0, `n_points_used=21` |
| **T2** | Retta `z = 0.5y + 3`, cov isotropa | `a=0.5`, `b=3`, `χ²≈0`, `σ_a > 0`, `n_points_used=31` |
| **T3** | Retta `z = 2y − 1`, cov anisotropa (`σ_y >> σ_z`) | Parametri invariati: con residui nulli la soluzione non dipende dai pesi |
| **T4** | Retta `z = 0` con outlier `Δz = 0.5 mm` al centro | Pull outlier ≈ `Δz/σ`, pull degli altri < 2, `χ²/ndof` dominato dall'outlier |
| **T5** | Traccia con 0, 2 punti validi, o N punti con tutti `valid=false` | Eccezione `std::invalid_argument`; nessuna eccezione per N≥3 con `valid=true` |
| **T6** | Covarianza degenere (`var = 0`) | Nessun crash grazie al floor `sd²_floor = 1e-6 mm²`; convergenza e parametri corretti |
| **T7** | Retta `z = 0`, cov non diagonale; outlier `Δz` al centro | `σ_{d,outlier} ≈ √(cov_zz)` (verifica formula `n̂ᵀ Σ n̂` con `a ≈ 0`) |

#### Suite TV — `is_trace_valid`

| Test | Scenario | Verifica |
|---|---|---|
| **TV1** | Tutti i punti `valid=true` | Valida per soglie 0.75 e 1.0 |
| **TV2** | 4 punti su 10 `valid=false` (60% validi) | Valida per soglia 50%, invalida per 75% e 100% |
| **TV3** | Traccia vuota | Ritorna `false` |

#### Suite TQ — `compute_Q`

| Test | Scenario | Verifica |
|---|---|---|
| **TQ1** | PSF ideale (`mu_y = r`, `mu_z = 0`, `on_detector=true`), 11 tracce, unfolding OFF | `Q ≈ 0`, `n_traces = 11`, `n_failed = 0`, `Q == Σ chi2_per_y0` |
| **TQ2** | Lista `y0_values` esplicita con 3 valori, unfolding OFF | `n_traces = 3`, `Q ≈ 0` |
| **TQ3** | Configurazione non presente nel database | Eccezione `std::invalid_argument` |
| **TQ4** | PSF con curvatura quadratica in z; due DB con σ_z diversa, unfolding OFF | `Q(σ_piccola) > Q(σ_grande)` — PSF stretta penalizza di più la non-linearità |
| **TQ5** | DB ripiegato (`mu_y` non monotona) vs DB lineare; confronto con/senza unfolding | `Q(ripiegata) > Q(lineare)` con unfolding ON; differenza ridotta con unfolding OFF |

---

## Configurazione raccomandata per `q_map`

Dopo ottimizzazione empirica sul dataset reale (`psf_data.root`), la configurazione seguente produce mappe stabili e sensibili alle non-linearità ottiche:

```bash
./build/analysis/psf_analysis/Release/q_map \
    --psf    output/psf/psf_data_1.root \
    --dy0    0.1                        \
    --unfold-dz 0.000002                \
    --point-frac 0.5                    \
    --trace-frac 0.50                   \
    --log
```

| Parametro | Valore | Motivazione |
|---|---|---|
| `--dy0 0.1` | Passo 0.1 mm su y₀ | Campionamento denso: ≈142 tracce per configurazione, sufficiente per mediare le fluttuazioni statistiche della PSF |
| `--unfold-dz 0.000002` | δz = 2×10⁻⁶ mm/passo | Valore molto piccolo rispetto al passo di campionamento della traccia: introduce un segnale di srotolamento appena sufficiente a rompere la degenerazione delle tracce ripiegate, senza distorcere il χ² delle tracce lineari. Con `trace_dt = 0.1` mm e ~101 punti, l'offset totale è ≈ 0.0002 mm, molto inferiore alla σ tipica della PSF (~0.1–0.5 mm) |
| `--point-frac 0.5` | Soglia 50% punti on-detector | Permette di valutare configurazioni in cui la PSF è parzialmente fuori dal fotocatodo — utile per mappare l'intera griglia senza perdere troppe celle |
| `--trace-frac 0.50` | Soglia 50% tracce valide | Coerente con `--point-frac`; filtra le configurazioni chiaramente inutilizzabili mantenendo quelle borderline |
| `--log` | Scala logaritmica | Q copre tipicamente 2–3 ordini di grandezza tra la configurazione ottimale e quelle scadenti; la scala log rende leggibile l'intera mappa |

> **Nota sul valore di `--unfold-dz`:** il default automatico (`δz = trace_L / (N-1) ≈ 0.1 mm`) è ottimale per rilevare ripiegamenti netti ma può sovrastimare Q per configurazioni con leggere curvature. Il valore `0.000002` è conservativo: preserva il ranking relativo delle configurazioni senza modificare significativamente i valori assoluti di Q.

---

## File di configurazione

`config/config.json` — letto da tutti i programmi di simulazione e analisi:

```json
{
  "x_min": 33.0,
  "x_max": 171.0,
  "dx": 3.0,
  "r1": 38.6,
  "h1": 12.5,
  "r2": 30.9,
  "h2": 16.3,
  "lower_percentile": 0.45,
  "upper_percentile": 0.0
}
```

| Parametro | Descrizione |
|---|---|
| `x_min`, `x_max` | Range di scansione [mm] |
| `dx` | Passo della scansione [mm] |
| `r1`, `h1` | Raggio e spessore della lente 75mm [mm] |
| `r2`, `h2` | Raggio e spessore della lente 60mm [mm] |
| `lower_percentile` | Soglia bassa per `plot2D` (esclude il 45% inferiore dell'efficienza) |
| `upper_percentile` | Soglia alta per `plot2D` (0 = nessun taglio superiore) |

I vincoli geometrici applicati durante la scansione:

```
x1_min = x_min - r1 + h1
x1_max = x_max - h2 - 3 - r1
x2_min = x1 + r1 + r2 + 3       (gap minimo 3 mm tra le lenti)
x2_max = x_max + r2 - h2
```

---

## Formato dei file ROOT

### events.root (optimization)

```
events.root
└── TTree "events"
    ├── x1        Double_t
    ├── x2        Double_t
    └── config_id Int_t
```

### lens.root (lens_simulation)

```
lens.root
├── TTree "Configurations"    (≈703 righe)
│   ├── config_id  Int_t
│   ├── x1         Double_t  [mm]
│   └── x2         Double_t  [mm]
│
├── TTree "Runs"              (≈71 000 righe)
│   ├── run_id     Int_t
│   ├── config_id  Int_t
│   ├── x_source   Float_t   [mm]
│   └── n_hits     Int_t
│
└── TTree "Hits"              (≈14M righe)
    ├── y_hit      Float_t   [mm]
    └── z_hit      Float_t   [mm]
```

### psf\_data.root (psf\_extractor)

```
psf_data.root
└── TTree "PSF"    (≈71 000 righe)
    ├── config_id         Int_t
    ├── x1, x2            Double_t   [mm]
    ├── y_source          Float_t    [mm]
    ├── mean_y, mean_z    Double_t   [mm]
    ├── cov_yy, cov_yz, cov_zz  Double_t  [mm²]
    ├── n_hits_filtered   Int_t
    └── n_hits_raw        Int_t
```

---

## Workflow completo

```bash
# ── 1. Build ────────────────────────────────────────────────────────────────
cmake -S . -B build/ -G "Ninja Multi-Config"
cmake --build build/ --config Release

# ── 2. Test unitari ─────────────────────────────────────────────────────────
./build/analysis/psf_analysis/Release/test_fit_trace
# oppure: cd build && ctest -R fit_trace_unit -V

# ── 3. Scansione efficienza geometrica (opzionale, complementare a Q) ───────
./build/Release/optimization_main -g geometry/main.gdml -o
./build/analysis/Release/plot2D
# → output/efficiency2D.png

# ── 4. Beam scan PSF (parallelizzato) ───────────────────────────────────────
./scripts/run.sh lens local --jobs $(nproc --all)
# → output/lens_simulation/lens_<timestamp>.root

# ── 5. Estrazione PSF ────────────────────────────────────────────────────────
./build/analysis/Release/psf_extractor \
    output/lens_simulation/lens.root \
    output/psf/psf_data.root
# → output/psf/psf_data.root

# ── 6. Mappa Q — ottimizzazione PSF-based ★ ─────────────────────────────────
./build/analysis/psf_analysis/Release/q_map \
    --y0-min 0.0 --y0-max 10.0 --dy0 0.1
# → output/psf_analysis/q_map.png
# Stampa su stdout: x1*, x2*, Q_min

# Con scala log e normalizzazione per una visione alternativa:
./build/analysis/psf_analysis/Release/q_map --log --norm \
    --output output/psf_analysis/q_map_norm_log.png

# ── 7. Analisi traccia per la configurazione ottimale ───────────────────────
# (sostituisci x1* e x2* con i valori stampati da q_map)
./build/analysis/psf_analysis/Release/trace_viewer \
    --x1 <x1*> --x2 <x2*> --y0 5.0
# → output/psf_analysis/trace_x1_<x1*>_x2_<x2*>_y0_5.0.png

# ── 8. Strumenti di diagnostica ─────────────────────────────────────────────
./build/analysis/Release/beam_scan_plot <x1*> <x2*>
./build/analysis/Release/m_c_creator   <x1*> <x2*> 5.0

# ── 9. Visualizzazione interattiva geometria ─────────────────────────────────
./build/Release/lens_simulation_main -g geometry/main.gdml -v
```
