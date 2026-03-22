# riptide_optimization

Simulazione Monte Carlo Geant4 per l'ottimizzazione del posizionamento di un sistema ottico a doppia lente (UVFS) davanti a un fotocatodo GaAsP. Il progetto fa parte dell'esperimento **RIPTIDE** e comprende due programmi principali — `optimization` e `lens_simulation` — quattro strumenti di analisi ROOT, e una libreria di analisi PSF (`psf_analysis`) con test unitari integrati.

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
12. [File di configurazione](#file-di-configurazione)
13. [Formato dei file ROOT](#formato-dei-file-root)
14. [Workflow completo](#workflow-completo)

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
├── scripts/                             # Script di lancio e monitoraggio
│   ├── run.sh                           # Lancia lens_simulation o optimization (locale o SSD)
│   └── monitor.sh                       # Monitoraggio progresso chunk paralleli
│
├── programs/                        # Entry point dei due eseguibili
│   ├── optimization_main.cpp        # Main per la scansione di efficienza
│   └── lens_simulation_main.cpp     # Main per la simulazione beam scan
│
├── src/
│   ├── common/                      # Codice condiviso tra i due programmi
│   │   ├── physics_list.cpp
│   │   └── primary_generator_action.cpp
│   ├── optimization/                # Sorgenti del programma optimization
│   │   ├── detector_construction.cpp
│   │   ├── action_initialization.cpp
│   │   ├── event_action.cpp
│   │   ├── run_action.cpp
│   │   ├── sensitive_detector.cpp
│   │   └── optimizer.cpp
│   └── lens_simulation/             # Sorgenti del programma lens_simulation
│       ├── detector_construction.cpp
│       ├── action_initialization.cpp
│       ├── event_action.cpp
│       ├── run_action.cpp
│       ├── sensitive_detector.cpp
│       └── lens_scan.cpp
│
├── include/                         # Header files (specchi di src/)
│   ├── common/
│   ├── optimization/
│   └── lens_simulation/
│
├── analysis/                        # Strumenti di analisi ROOT (standalone)
│   ├── plot2D.cpp                   # Mappa 2D efficienza geometrica
│   ├── beam_scan_plot.cpp           # Grafici posizione fotoni vs posizione sorgente
│   ├── m_c_creator.cpp              # Istogramma 2D hit su detector per un run
│   ├── psf_extractor.cpp            # Media e matrice di covarianza PSF per tutti i run
│   ├── CMakeLists.txt
│   └── psf_analysis/                # Libreria PSF + analisi traccia
│       ├── psf_interpolator.hpp     # API pubblica: strutture dati, dichiarazioni
│       ├── psf_interpolator.cpp     # Implementazione: load, interpolate, build_trace, fit_trace
│       ├── trace_viewer.cpp         # Eseguibile: visualizza traccia media con ellissi di cov
│       ├── test_fit_trace.cpp       # Test unitari analitici per fit_trace()
│       └── CMakeLists.txt
│
├── geometry/                        # Geometria GDML
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
    ├── lens_simulation/
    │   └── lens.root
    ├── mean_covariance_maps/
    └── psf/
        ├── psf_data.root
        └── psf_analysis/            # Output di trace_viewer
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

# Gli eseguibili vengono creati in:
#   build/Release/optimization_main
#   build/Release/lens_simulation_main
#   build/analysis/Release/plot2D
#   build/analysis/Release/beam_scan_plot
#   build/analysis/Release/m_c_creator
#   build/analysis/Release/psf_extractor
#   build/analysis/psf_analysis/Release/trace_viewer
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
lsblk                                          # identifica il device
sudo mkdir -p /mnt/external_ssd
sudo mount -o noatime,nodiratime,discard /dev/nvme1n1p1 /mnt/external_ssd
mountpoint -q /mnt/external_ssd && echo "OK"
```

Opzioni di mount consigliate per I/O intensivo: `noatime,nodiratime,discard`.

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

La libreria `psf_analysis` (in `analysis/psf_analysis/`) implementa la catena analitica PSF-based descritta nella documentazione di tesi. È compilata come libreria statica e usata da `trace_viewer` e `test_fit_trace`.

### API pubblica (`psf_interpolator.hpp`)

#### Strutture dati principali

| Struttura | Descrizione |
|---|---|
| `PSFPoint` | Un punto del database PSF: `y_source`, `mu_y`, `mu_z`, covarianza |
| `Cov2` | Matrice di covarianza 2×2: `yy`, `yz`, `zz` |
| `PSFValue` | Risultato interpolato: `mu_y`, `mu_z`, `cov` |
| `TracePoint` | Punto della traccia sul detector: `t`, `r`, `mu_y`, `mu_z`, `cov` |
| `LensConfig` | Chiave di configurazione `(x1, x2)` con confronto a tolleranza 1e-4 mm |
| `LineFitResult` | Risultato completo del fit ODR (vedi sotto) |

#### Funzioni

**`load_psf_database(path)`** — carica `psf_data.root` in memoria come `PSFDatabase`.

**`find_nearest_config(cfg, db)`** — trova la configurazione più vicina nel database (distanza euclidea in `(x1, x2)`).

**`interpolate(r, cfg, db)`** — interpola `mu_y`, `mu_z` e `Sigma` per un raggio `r` arbitrario. Usa interpolazione lineare tra i due punti della griglia adiacenti; clamp agli estremi senza estrapolazione.

**`build_trace(y0, cfg, db, L=10, dt=0.1)`** — costruisce la traccia media sul detector per una traccia ideale a distanza `y0` dall'asse ottico, parametrizzata da `t ∈ [-L/2, +L/2]` con step `dt`.

**`fit_trace(trace, max_iter=20, tol=1e-8)`** — fit lineare pesato ODR della traccia (vedi sezione dedicata).

#### `fit_trace` — Fit ODR iterativo

Esegue il fit della retta `z = a·y + b` sui punti `{(mu_y_i, mu_z_i)}` usando le matrici di covarianza `Σ_i` come peso statistico, tramite Orthogonal Distance Regression (ODR) con schema IRLS.

**Algoritmo:**

1. Stima iniziale `(a, b)` con OLS non pesato.
2. Per ogni iterazione:
   - Calcola il vettore normale `n̂ = (-a, 1)ᵀ / √(1+a²)`
   - Calcola i pesi `w_i = 1 / (n̂ᵀ Σ_i n̂) = 1 / (n_y²·σ²_y + 2·n_y·n_z·σ_yz + n_z²·σ²_z)`
   - Risolve il sistema normale 2×2 pesato in forma chiusa per `(a, b)` — costo O(N)
   - Controlla convergenza su `|Δa| < tol`
3. Calcola χ², residui e pull finali con distanza perpendicolare esatta.

**`LineFitResult` — campi:**

| Campo | Tipo | Descrizione |
|---|---|---|
| `a`, `b` | `double` | Parametri della retta `z = a·y + b` |
| `sigma_a`, `sigma_b` | `double` | Incertezze sui parametri |
| `cov_ab` | `double` | Covarianza tra `a` e `b` |
| `chi2` | `double` | χ² totale |
| `ndof` | `int` | Gradi di libertà = N − 2 |
| `chi2_ndof` | `double` | χ² ridotto = χ²/ndof |
| `n_iter` | `int` | Iterazioni eseguite |
| `converged` | `bool` | `true` se `\|Δa\| < tol` |
| `residuals` | `vector<double>` | Distanza perpendicolare `d_i` [mm] |
| `residual_sig` | `vector<double>` | `σ_{d,i} = √(n̂ᵀ Σ_i n̂)` [mm] |
| `pull` | `vector<double>` | Pull normalizzato `d_i / σ_{d,i}` |

**Utilizzo tipico:**

```cpp
#include "psf_interpolator.hpp"

// Carica il database PSF (una volta sola)
auto db = riptide::load_psf_database("output/psf/psf_data.root");

// Trova la configurazione disponibile più vicina
riptide::LensConfig cfg{94.9, 186.4};
auto actual = riptide::find_nearest_config(cfg, db);

// Costruisce la traccia media per y0 = 5.0 mm
auto trace = riptide::build_trace(5.0, actual, db);

// Fit lineare pesato ODR
auto result = riptide::fit_trace(trace);

std::cout << "a = " << result.a << " ± " << result.sigma_a << "\n";
std::cout << "b = " << result.b << " ± " << result.sigma_b << " mm\n";
std::cout << "chi2/ndof = " << result.chi2_ndof << "\n";
```

### trace\_viewer

Visualizza la traccia media sul detector con ellissi di covarianza, colorate per parametro `t` lungo la traccia (palette Rainbow: blu = inizio, rosso = fine).

```bash
./build/analysis/psf_analysis/Release/trace_viewer \
    --x1 94.9 --x2 186.4 --y0 5.0

# Opzioni complete:
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
| `--output` | auto (con timestamp) | Path immagine PNG |
| `--dt` | `0.1` | Step traccia [mm] |
| `--L` | `10.0` | Lunghezza traccia [mm] |
| `--sigma` | `1.0` | Scala delle ellissi di covarianza (in unità σ) |

---

## Test unitari

### Esecuzione

```bash
# Esegue solo i test PSF
./build/analysis/psf_analysis/Release/test_fit_trace

# Oppure tramite CTest dalla directory di build
cd build && ctest -R fit_trace_unit -V
```

L'eseguibile ritorna `0` se tutti i test passano, `1` altrimenti.

### Casi di test (`test_fit_trace.cpp`)

I test costruiscono `TracePoint` sintetici con soluzione analitica nota — nessun file ROOT necessario.

| Test | Scenario | Verifica |
|---|---|---|
| **T1** | Retta `z = y`, cov isotropa uniforme | `a=1`, `b=0`, `χ²≈0`, pull tutti ≈ 0 |
| **T2** | Retta `z = 0.5y + 3`, cov isotropa | `a=0.5`, `b=3`, `χ²≈0`, `σ_a > 0` |
| **T3** | Retta `z = 2y − 1`, cov anisotropa (`σ_y >> σ_z`) | Parametri invariati: con residui nulli la soluzione non dipende dai pesi |
| **T4** | Retta `z = 0` con outlier `Δz = 0.5 mm` al centro | Pull outlier ≈ `Δz/σ`, pull degli altri < 2, `χ²/ndof` dominato dall'outlier |
| **T5** | Traccia con 0, 2 o 3 punti | Eccezione `std::invalid_argument` per N < 3; nessuna eccezione per N = 3 |
| **T6** | Covarianza degenere (`var = 0`) | Nessun crash grazie al floor `sd²_floor = 1e-6 mm²`; convergenza e parametri corretti |
| **T7** | Retta `z = 0`, cov non diagonale; outlier `Δz` al centro | `σ_{d,outlier} ≈ √(cov_zz)` (verifica formula `n̂ᵀ Σ n̂` con `a ≈ 0`) |

---

## File di configurazione

`config/config.json` — letto da entrambi i programmi di simulazione e da `plot2D`:

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
| `lower_percentile` | Soglia bassa per `plot2D` (esclude il 45% inferiore) |
| `upper_percentile` | Soglia alta per `plot2D` (0 = nessun taglio superiore) |

I vincoli geometrici applicati durante la scansione sono:

```
x1_min = x_min - r1 + h1
x1_max = x_max - h2 - 1 - r1
x2_min = x1 + r1 + r2 + 1       (gap minimo 1 mm tra le lenti)
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
# 1. Build
cmake -S . -B build/ -G "Ninja Multi-Config"
cmake --build build/ --config Release

# 2. Test unitari (verifica che fit_trace funzioni correttamente)
./build/analysis/psf_analysis/Release/test_fit_trace
# oppure: cd build && ctest -R fit_trace_unit -V

# 3a. Scansione efficienza — output locale
./build/Release/optimization_main -g geometry/main.gdml -o
# → output/events.root

# 3b. Mappa 2D efficienza
./build/analysis/Release/plot2D
# → output/efficiency2D.png

# 4a. Beam scan dettagliato — output locale (parallelizzato)
./scripts/run.sh lens local --jobs $(nproc --all)
# → output/lens_simulation/lens_<timestamp>.root

# 4b. Estrazione PSF
./build/analysis/Release/psf_extractor \
    output/lens_simulation/lens.root \
    output/psf/psf_data.root
# → output/psf/psf_data.root

# 4c. Visualizzazione traccia con ellissi di covarianza
./build/analysis/psf_analysis/Release/trace_viewer \
    --x1 94.9 --x2 186.4 --y0 5.0
# → output/psf_analysis/trace_x1_94.9_x2_186.4_y0_5.0.png

# 4d. Fit ODR della traccia (uso programmatico della libreria)
#     Vedere esempio in sezione "Libreria psf_analysis"

# 5. Strumenti di diagnostica
./build/analysis/Release/beam_scan_plot 94.9 186.4
./build/analysis/Release/m_c_creator   94.9 186.4 5.0

# 6. Visualizzazione interattiva geometria
./build/Release/optimization_main    -g geometry/main.gdml -v
./build/Release/lens_simulation_main -g geometry/main.gdml -v
```
