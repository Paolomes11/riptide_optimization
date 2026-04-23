# riptide_optimization

Simulazione Monte Carlo Geant4 per l'ottimizzazione del posizionamento di un sistema ottico a doppia lente (UVFS) davanti a un fotocatodo GaAsP. Il progetto fa parte dell'esperimento **RIPTIDE** e comprende due programmi principali — `optimization` e `lens_simulation` — sei strumenti di analisi ROOT, e una libreria di analisi PSF (`psf_analysis`) con test unitari integrati.

---

## Tool: lens\_cutter

**Eseguibile**: `lens_cutter_main`

**Scopo**: Gestione e generazione di geometrie per lenti commerciali Thorlabs. Permette di elencare le lenti disponibili e generare snippet GDML o solidi Geant4.

### Utilizzo
```bash
# Elenca tutte le lenti Thorlabs disponibili
./build/Release/lens_cutter_main --list

# Genera lo snippet GDML per una lente specifica
./build/Release/lens_cutter_main --id LB4592
```

### Funzionamento

Il tool legge le specifiche da file TSV nel formato Thorlabs. Sono supportati due formati:

- **Biconvesse** (`thorlabs_biconvex.tsv`): colonne `Item # | Diameter | Focal Length | Radius of Curvature | Center Thickness | Edge Thickness | Back Focal Length`. Il tipo viene rilevato automaticamente dall'assenza della colonna `Rotation_deg`.

- **Plano-convesse** (`thorlabs_planoconvex.tsv`): stesso schema con colonna aggiuntiva `Rotation_deg` [deg]. Indica la rotazione attorno all'asse Y applicata al volume fisico: `0` = lato curvo verso la sorgente, `180` = lato piano verso la sorgente. La rotazione è propagata automaticamente al `G4VPhysicalVolume` in `DetectorConstruction`.

Per caricare più cataloghi contemporaneamente:
```bash
./build/Release/lens_cutter_main --list \
    --catalog lens_cutter/lens_data/thorlabs_planoconvex.tsv
```

Le simulazioni `optimization` e `lens_simulation` accettano ID da entrambi i cataloghi tramite `--lens75-id` e `--lens60-id`; la geometria (solido + rotazione + offset centro geometrico) viene impostata automaticamente in base al tipo rilevato.

---

## Indice

1. [Requisiti](#requisiti)
2. [Struttura del progetto](#struttura-del-progetto)
3. [Setup fisico](#setup-fisico)
4. [Build](#build)
5. [Programma: optimization](#programma-optimization)
6. [Programma: lens\_simulation](#programma-lens_simulation)
7. [Output su SSD esterna](#output-su-ssd-esterna)
8. [Esecuzione Parallela e Dashboard](#esecuzione-parallela-e-dashboard)
9. [Strumenti di analisi](#strumenti-di-analisi)
10. [Libreria psf\_analysis](#libreria-psf_analysis)
11. [Test unitari](#test-unitari)
12. [Configurazione raccomandata per `q_map`](#configurazione-raccomandata-per-q_map)
13. [File di configurazione](#file-di-configurazione)
14. [Formato dei file ROOT](#formato-dei-file-root)
15. [Workflow completo](#workflow-completo)
16. [Tool: exp2\_main — Analisi immagini FITS laser](#tool-exp2_main--analisi-immagini-fits-laser)
17. [Riferimenti metodologici](#riferimenti-metodologici)

---

## Requisiti

| Dipendenza | Versione minima | Note |
|---|---|---|
| CMake | 3.26 | |
| C++ | 17 | |
| Geant4 | 11.x | Con supporto `ui_all`, `vis_all`, GDML, ottica |
| ROOT | 6.x | Componenti: `Core`, `Hist`, `RIO`, `Graf`, `Graf3d` |
| spdlog | qualsiasi | Logging strutturato |
| nlohmann/json | (header-only, incluso) | Lettura configurazione |
| lyra | (header-only, incluso in `external/`) | Parsing CLI |

---

## Struttura del progetto
```
riptide_optimization/
├── scripts/
│   ├── run.sh                           # Lancia lens_simulation o optimization (locale o SSD)
│   ├── monitor.sh                       # Monitoraggio progresso chunk paralleli (bash)
│   └── dashboard.py                     # Dashboard Python/rich (prioritario su monitor.sh)
│
├── programs/
│   ├── optimization_main.cpp            # Main per la scansione di efficienza geometrica
│   ├── lens_simulation_main.cpp         # Main per la simulazione beam scan PSF
│   └── lens_cutter_main.cpp             # Main per la gestione delle lenti Thorlabs
│
├── src/
│   ├── common/                          # Codice condiviso
│   │   ├── physics_list.cpp
│   │   ├── primary_generator_action.cpp
│   │   └── stepping_action.cpp
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
│   │   ├── physics_list.hpp
│   │   ├── primary_generator_action.hpp
│   │   └── stepping_action.hpp
│   ├── optimization/
│   │   ├── action_initialization.hpp
│   │   ├── detector_construction.hpp
│   │   ├── event_action.hpp
│   │   ├── optimizer.hpp
│   │   ├── run_action.hpp
│   │   └── sensitive_detector.hpp
│   └── lens_simulation/
│       ├── action_initialization.hpp
│       ├── detector_construction.hpp
│       ├── event_action.hpp
│       ├── lens_scan.hpp
│       ├── run_action.hpp
│       └── sensitive_detector.hpp
│
├── analysis/
│   ├── plot2D.cpp                       # Mappa 2D efficienza geometrica
│   ├── plot3D.cpp                       # Visualizzazione 3D multi-lente
│   ├── beam_scan_plot.cpp               # Posizione fotoni vs posizione sorgente (3D + 2D slices)
│   ├── m_c_creator.cpp                  # Istogramma 2D hit su detector per un run
│   ├── psf_extractor.cpp                # Estrazione media e covarianza PSF per tutti i run
│   ├── CMakeLists.txt
│   └── psf_analysis/                    # Libreria PSF + analisi traccia + ottimizzazione
│       ├── psf_interpolator.hpp         # API pubblica: strutture dati, dichiarazioni
│       ├── psf_interpolator.cpp         # load, interpolate, build_trace_3d, fit_trace, compute_Q
│       ├── trace_viewer.cpp             # Visualizza traccia 2D/3D con ellissi di covarianza
│       ├── q_map.cpp                    # Mappa 2D di Q(x1,x2) oppure copertura geometrica
│       ├── chi2_map.cpp                 # Mappa 2D linearità risposta (fit piano μ_y,z)
│       ├── test_fit_trace.cpp           # Test unitari per fit_trace() e compute_Q()
│       └── CMakeLists.txt
│
├── lens_cutter/
│   ├── include/lens_cutter.hpp
│   ├── src/lens_cutter.cpp
│   └── lens_data/
│       ├── thorlabs_biconvex.tsv
│       └── thorlabs_planoconvex.tsv
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
│   ├── lens_profile.mac
│   ├── run.mac
│   └── vis.mac
│
├── config/
│   ├── config.json
│   └── config_profile.json
│
├── external/
│   └── lyra/
│
└── output/
    ├── lens_simulation/
    ├── mean_covariance_maps/
    ├── psf/
    └── psf_analysis/
```

---

## Setup fisico

Il sistema simulato è composto da tre elementi posizionati lungo l'asse X, all'interno di un volume mondo cubico di 1000×1000×1000 mm³ riempito d'aria:
```
Sorgente fotoni  →  [Lente 75mm]  →  [Lente 60mm]  →  [Fotocatodo GaAsP 16×16mm]
      (GPS)           (UVFS)           (UVFS)              (sensore)
         x=0        x ≈ 14–170mm    x ≈ 45–186mm           x ≈ 180mm
```

**Lente 75 mm** (`lens75`): ellissoide in UV Fused Silica (UVFS), raggio 38.6 mm, spessore 12.5 mm. Sostituibile con qualsiasi lente Thorlabs tramite `--lens75-id`.

**Lente 60 mm** (`lens60`): ellissoide UVFS, raggio 30.9 mm, spessore 16.3 mm. Sostituibile tramite `--lens60-id`.

**Fotocatodo GaAsP**: lastra quadrata 16×16×0.01 mm, indice di rifrazione 3.5–3.8 nel range 2–4 eV, lunghezza di assorbimento ~1 µm.

**Sorgente**: fotoni ottici a 2.5 eV (≈ 496 nm), generati da GPS Geant4 su una sferetta di raggio 0.1 mm con distribuzione angolare isotropa diretta lungo +X (`mintheta=0`, `maxtheta=90 deg`).

**Processi ottici attivi**: rifrazione, riflessione (Fresnel), assorbimento. Scintillazione, Cherenkov, WLS e WLS2 sono esplicitamente disabilitati in `PhysicsList` perché i fotoni sono generati direttamente dalla GPS.

**SteppingAction**: uccide i fotoni che si allontanano nella direzione −X (x < −1 mm con px < 0) o che escono dalla regione utile (|y| > 150 mm o |z| > 150 mm), evitando tracking inutile fuori dall'apertura delle lenti.

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
#   build/Release/lens_cutter_main
#   build/analysis/Release/plot2D
#   build/analysis/Release/plot3D
#   build/analysis/Release/beam_scan_plot
#   build/analysis/Release/m_c_creator
#   build/analysis/Release/psf_extractor
#   build/analysis/psf_analysis/Release/trace_viewer
#   build/analysis/psf_analysis/Release/q_map
#   build/analysis/psf_analysis/Release/chi2_map
#   build/analysis/psf_analysis/Release/test_fit_trace
```

> **Nota**: CMake crea automaticamente le cartelle `output/`, `output/lens_simulation/`, `output/mean_covariance_maps/`, `output/psf/` e `output/psf_analysis/` durante la configurazione.

> **Build Release**: la build Release aggiunge automaticamente `-march=native -O3 -ffast-math` tramite il `CMakeLists.txt` root.

---

## Programma: optimization

**Eseguibile**: `optimization_main`

**Scopo**: scansione grezza dell'efficienza geometrica del sistema a doppia lente. Per ogni coppia di posizioni (x1, x2) delle lenti, esegue una simulazione con fotoni emessi da una sorgente rettangolare e conta quanti raggiungono il fotocatodo.

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
| `--all-lenses` | flag | Scansione di **tutte** le combinazioni di lenti Thorlabs |
| `--lens75-id` | string | ID lente Thorlabs per L1 (es. `LB4592`) |
| `--lens60-id` | string | ID lente Thorlabs per L2 (es. `LB4553`) |
| `--output` | path | File ROOT di output (default: `output/events.root`) |
| `--config` | path | File `config.json` (default: `config/config.json`) |
| `--ssd` | flag | Output sull'SSD esterna con timestamp automatico |
| `--ssd-mount` | path | Mount point dell'SSD (default: `/mnt/external_ssd`) |

### Output (`events.root`)

**`Configurations`**: una riga per ogni configurazione di lenti testata.

| Branch | Tipo | Descrizione |
|---|---|---|
| `config_id` | `Int_t` | Identificativo unico della configurazione |
| `x1` | `Double_t` | Posizione della lente 1 (75mm) [mm] |
| `x2` | `Double_t` | Posizione della lente 2 (60mm) [mm] |
| `lens75_id` | `String` | ID Thorlabs lente 1 |
| `lens60_id` | `String` | ID Thorlabs lente 2 |

**`Efficiency`**: una riga per ogni configurazione, con il conteggio dei fotoni.

| Branch | Tipo | Descrizione |
|---|---|---|
| `config_id` | `Int_t` | Identificativo della configurazione |
| `n_photons` | `Int_t` | Numero di fotoni generati |
| `n_hits` | `Int_t` | Numero di fotoni che hanno raggiunto il detector |

---

## Programma: lens\_simulation

**Eseguibile**: `lens_simulation_main`

**Scopo**: simulazione dettagliata del beam scan per la caratterizzazione della PSF. Per ogni configurazione di lenti e per ogni coppia (x\_source, y\_source) sulla griglia, registra la posizione (y, z) di ogni fotone sul fotocatodo in vettori per-run.

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
| `--lens75-id` | string | ID lente Thorlabs per L1 |
| `--lens60-id` | string | ID lente Thorlabs per L2 |
| `--output` | path | File ROOT di output (default: `output/lens_simulation/lens.root`) |
| `--config` | path | File `config.json` (default: `config/config.json`) |
| `--ssd` | flag | Output sull'SSD esterna |
| `--ssd-mount` | path | Mount point dell'SSD |

### Griglia di campionamento

La sorgente viene spostata su una griglia 2D definita in `config.json`:

- **Asse X**: da `source_x_min` a `source_x_max` con passo `source_dx`
- **Asse Y** (radiale): da `source_y_min` a `source_y_max` con passo `source_dy`

Il limite superiore `source_y_max = 10√2 ≈ 14.14 mm` garantisce che tracce di lunghezza 10 mm a qualsiasi y₀ ≤ 10 mm siano completamente coperte dalla griglia PSF (raggio massimo ≈ 11.18 mm ≪ 14.14 mm).

### Output (`lens.root`)

**`Configurations`**: identico a `optimization`.

**`Runs`**: una riga per ogni esecuzione (configurazione × posizione sorgente).

| Branch | Tipo | Descrizione |
|---|---|---|
| `run_id` | `Int_t` | Indice run globale |
| `config_id` | `Int_t` | Configurazione di appartenenza |
| `x_source` | `Float_t` | Posizione X della sorgente [mm] |
| `y_source` | `Float_t` | Posizione Y della sorgente [mm] |
| `n_hits` | `Int_t` | Fotoni rilevati in questo run |
| `y_hits` | `vector<float>` | Coordinate Y hit sul fotocatodo [mm] |
| `z_hits` | `vector<float>` | Coordinate Z hit sul fotocatodo [mm] |

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

## Esecuzione Parallela e Dashboard

Il framework supporta la parallelizzazione tramite lo script `scripts/run.sh`, che divide il carico di lavoro in chunk indipendenti eseguiti in parallelo e li fonde con `hadd` al termine.

### Avvio
```bash
./scripts/run.sh lens local --jobs $(nproc --all)
./scripts/run.sh lens ssd   --jobs 8
./scripts/run.sh opt  local --jobs 4
./scripts/run.sh opt  local --jobs 4 --all-lenses
```

### Dashboard di Monitoraggio

Lo script tenta prima di avviare `scripts/dashboard.py` (Python + `rich`); se non disponibile, ricade su `scripts/monitor.sh` (bash puro). La dashboard mostra:

- Progresso globale e per-chunk con barre di avanzamento.
- Tempo trascorso e stima ETA calcolata sul rate corrente.
- Stato di ogni job: `[wait]`, `[ >> ]` (in corso), `[DONE]`.

### Gestione processi e cleanup

Premendo `Ctrl+C`, il trap `SIGINT`/`SIGTERM` termina tutti i processi Geant4 in background, rimuove i config temporanei in `/tmp/riptide_chunks_*` e i file ROOT parziali.

### Chunking automatico

`run.sh` delega a uno script Python embedded il calcolo di tutte le coppie (x1, x2) valide (rispettando i vincoli di non-collisione geometrica) e la loro distribuzione bilanciata in N chunk. Ogni chunk riceve un `config_id_offset` e un `run_id_offset` univoci per garantire ID globalmente monotoni dopo il merge con `hadd`.

---

## Strumenti di analisi

### plot2D

Genera una mappa di calore 2D dell'efficienza in funzione delle posizioni (x₁, x₂).
```bash
./build/analysis/Release/plot2D \
    -i output/events.root \
    --lens1 LB4553 --lens2 LB4592 \
    --low 0.05 --high 0.05
# → output/efficiency2D.png
```

| Opzione | Default | Descrizione |
|---|---|---|
| `-i`, `--input` | `output/events.root` | File ROOT di input |
| `-c`, `--config` | `config/config.json` | File di configurazione |
| `-o`, `--output` | `output/efficiency2D.png` | File PNG di output |
| `--lens1`, `--lens2` | prima coppia disponibile | Selezione coppia di lenti |
| `--low`, `--high` | da `config.json` | Percentili da escludere dalla scala colori |

### plot3D

Visualizzazione 3D multi-lente per la scansione `--all-lenses`. Mostra un canvas con un plot 3D per ogni modello di prima lente: asse X = modello seconda lente, asse Y = x₁, asse Z = x₂, colore = efficienza. Una linea rossa traccia il valor medio pesato (x̄₁, x̄₂) per ogni seconda lente.
```bash
./build/analysis/Release/plot3D \
    -i output/events.root \
    --low 0.1 --high 0.1

# Solo per una prima lente, con griglia 2D per ogni seconda lente:
./build/analysis/Release/plot3D \
    -i output/events.root \
    --lens1 LB4553 --2d
```

### beam\_scan\_plot

Per una coppia (x1, x2) fissata, calcola posizione media e deviazione standard dei fotoni al variare della posizione sorgente (x\_src, y\_src). Applica filtro outlier 2σ iterativo (2 iterazioni). Produce un grafico 3D con piano di fit sovrapposto. Con l'opzione `x0_filters`, produce anche sezioni 2D a x\_source fissato.
```bash
./build/analysis/Release/beam_scan_plot <x1> <x2> [file.root] [x0_1,x0_2,...]
./build/analysis/Release/beam_scan_plot 94.9 186.4
./build/analysis/Release/beam_scan_plot 94.9 186.4 output/lens_simulation/lens.root -30,0,30
# → output/lens_simulation/beam_scan_3D_x1_94.90_x2_186.40.png
# → output/lens_simulation/beam_scan_2D_x1_94.90_x2_186.40.png  (se x0_filters forniti)
```

### m\_c\_creator

Istogramma 2D delle hit per una terna (x1, x2, y0) fissata, con filtro ellittico iterativo a 2σ (4 iterazioni).
```bash
./build/analysis/Release/m_c_creator <x1> <x2> <y0>
./build/analysis/Release/m_c_creator 94.9 186.4 5.0
# → output/mean_covariance_maps/detector_hits_config_<id>_y0_5.0.png
```

### psf\_extractor

Estrae media bidimensionale e matrice di covarianza della PSF per tutti i run. Applica filtro outlier ellittico (distanza di Mahalanobis) a 3σ con 4 iterazioni. Marca ogni run con il flag `on_detector` (true se `n_hits_filtered >= min_hits`, default 50).
```bash
./build/analysis/Release/psf_extractor [input.root] [output.root] [min_hits]
./build/analysis/Release/psf_extractor \
    output/lens_simulation/lens.root \
    output/psf/psf_data.root \
    50
```

Output `psf_data.root` — TTree `PSF`:

| Branch | Tipo | Descrizione |
|---|---|---|
| `config_id` | `Int_t` | Indice configurazione |
| `x1`, `x2` | `Double_t` | Posizioni lenti [mm] |
| `x_source` | `Float_t` | Posizione X sorgente [mm] |
| `y_source` | `Float_t` | Posizione Y sorgente [mm] |
| `mean_y`, `mean_z` | `Double_t` | Media PSF [mm] |
| `cov_yy`, `cov_yz`, `cov_zz` | `Double_t` | Matrice di covarianza [mm²] |
| `n_hits_filtered` | `Int_t` | Hit dopo filtro Mahalanobis |
| `n_hits_raw` | `Int_t` | Hit prima del filtro |
| `on_detector` | `Bool_t` | `true` se `n_hits_filtered >= min_hits` |

---

## Libreria psf\_analysis

La libreria `psf_analysis` (in `analysis/psf_analysis/`) implementa la catena analitica PSF-based completa: caricamento database, interpolazione 2D, costruzione tracce 3D, fit ODR pesato, temporal unfolding, calcolo Q e copertura. È compilata come libreria statica e usata da `trace_viewer`, `q_map`, `chi2_map` e `test_fit_trace`.

### API pubblica (`psf_interpolator.hpp`)

#### Strutture dati principali

| Struttura | Descrizione |
|---|---|
| `PSFPoint` | Un punto del database PSF: `x_source`, `y_source`, `mu_y`, `mu_z`, covarianza, `on_detector`, `n_hits` |
| `Cov2` | Matrice di covarianza 2×2: `yy`, `yz`, `zz` |
| `PSFValue` | Risultato interpolato: `mu_y`, `mu_z`, `cov`, `on_detector`, `n_hits_interp` |
| `TracePoint` | Punto della traccia sul detector: `t`, `r`, `x_src`, `y_src`, `z_src`, `mu_y`, `mu_z`, `cov`, `valid`, `n_hits` |
| `LensConfig` | Chiave di configurazione `(x1, x2)` con confronto a tolleranza 1e-4 mm |
| `LineFitResult` | Risultato completo del fit ODR: `a`, `b`, `sigma_a`, `sigma_b`, `cov_ab`, `chi2`, `ndof`, `chi2_ndof`, `n_iter`, `converged`, residui, pull |
| `QConfig` | Parametri per `compute_Q`: dimensioni scintillatore, `n_tracks`, `trace_dt`, soglie, unfolding |
| `QResult` | Risultato di `compute_Q`: `Q`, `n_traces`, `n_failed`, `n_invalid`, `config_valid`, chi2 per-traccia |
| `CoverageResult` | Risultato di `compute_coverage`: `coverage`, `n_y0_evaluated`, `config_valid` |

#### Funzioni

**`load_psf_database(path)`** — carica `psf_data.root` in memoria come `PSFDatabase`. I punti vengono ordinati per (x\_source, y\_source) crescente per facilitare l'interpolazione bilineare.

**`find_nearest_config(cfg, db)`** — trova la configurazione più vicina nel database (distanza euclidea in `(x1, x2)`). Logga un warning se la distanza supera 1e-4 mm.

**`interpolate(x, y, cfg, db)`** — interpolazione bilineare di `mu_y`, `mu_z` e `Σ` per una sorgente in posizione `(x, y)`. Clamp agli estremi; gestisce il caso vicino all'asse (y < y\_min) con interpolazione lineare verso l'origine e covarianza isotropizzata.

**`build_trace_3d(p1, p2, cfg, db, dt)`** — costruisce la traccia media sul detector per un segmento 3D P1→P2. Per ogni punto lungo il segmento, calcola la distanza radiale r dall'asse X, chiama `interpolate(x, r, ...)` e ruota media e covarianza dell'angolo azimutale φ = atan2(z, y) tramite la matrice di rotazione R(φ).

**`build_trace(y0, cfg, db, L, dt)`** — wrapper di compatibilità: costruisce la traccia per un segmento rettilineo lungo X a distanza y₀ dall'asse, con P1 = (−L/2, y₀, 0) e P2 = (+L/2, y₀, 0).

**`is_trace_valid(trace, point_valid_fraction)`** — ritorna `true` se la frazione di `TracePoint` con `valid == true` supera la soglia (default 0.75).

**`fit_trace(trace, min_hits_per_point, max_iter, tol)`** — fit lineare pesato ODR della traccia media (vedi sezione dedicata).

**`compute_Q(cfg, db, qcfg, include_non_converged)`** — calcola la funzione di qualità Q(x1,x2) su tracce casuali 3D nello scintillatore (vedi sezione dedicata).

**`compute_coverage(cfg, db, qcfg)`** — calcola la copertura geometrica media: frazione di `TracePoint` con `n_hits >= min_hits_per_point`, mediata su tutte le tracce campionate. Non esegue nessun fit ODR.

#### `fit_trace` — Fit ODR iterativo

Esegue il fit della retta `z = a·y + b` sui punti `{(mu_y_i, mu_z_i)}` con `valid == true` e `n_hits >= min_hits_per_point`, usando le matrici di covarianza `Σ_i` come peso statistico tramite Orthogonal Distance Regression (IRLS).

**Algoritmo:**

1. Stima iniziale `(a, b)` con OLS non pesato.
2. Per ogni iterazione: calcola il vettore normale `n̂ = (−a, 1)/‖…‖`, poi i pesi `w_i = 1 / (n̂ᵀ Σ_i n̂)` con floor `1e-6 mm²`, risolve il sistema WLS 2×2 in forma chiusa, controlla convergenza su `|Δa| < tol`.
3. Calcola χ², ndof, residui perpendicolari e pull finali.

Output: `LineFitResult` completo con `n_points_used`, `converged`, `n_iter`.

#### `compute_Q` — Funzione di qualità

Genera `n_tracks` tracce casuali nello scintillatore (punti su facce casuali del parallelepipedo `scint_x × scint_y × scint_z`), per ognuna:

1. Costruisce la traccia 3D con `build_trace_3d`.
2. Verifica la validità (≥75% dei punti on-detector).
3. Applica il temporal unfolding (se abilitato).
4. Esegue `fit_trace`.
5. Accumula `chi2_ndof`.

Restituisce la media del Chi-squared ridotto sulle tracce valide:
```
Q(x1, x2) = (1 / n_traces) · Σ chi²_ndof
```

La configurazione è marcata `config_valid = false` se la frazione di tracce valide scende sotto `trace_valid_fraction`.

**`QConfig` — campi principali:**

| Campo | Default | Descrizione |
|---|---|---|
| `scint_x`, `scint_y`, `scint_z` | 60, 20, 20 mm | Dimensioni scintillatore |
| `n_tracks` | 100 | Numero di tracce casuali per configurazione |
| `trace_dt` | 0.1 mm | Passo campionamento traccia |
| `min_hits_per_point` | 10.0 | Hit minime per considerare un punto PSF valido |
| `trace_valid_fraction` | 0.75 | Soglia tracce valide per configurazione valida |
| `apply_temporal_unfolding` | true | Abilita/disabilita il temporal unfolding |
| `z_unfold_step` | 0.0 (auto) | Passo di srotolamento fisso [mm/passo]; 0 = L/(N−1) |

### trace\_viewer

Visualizza la traccia media sul detector con ellissi di covarianza, colorate per parametro `t` lungo la traccia (palette Rainbow: blu = inizio, rosso = fine). Supporta sia la modalità 2D (sorgente a y₀ fissato) che 3D (segmento P1→P2 arbitrario). Opzionalmente esegue il fit ODR e mostra χ²/ndof nel pannello info.
```bash
# Modalità 2D
./build/analysis/psf_analysis/Release/trace_viewer \
    --x1 94.9 --x2 186.4 --y0 5.0

# Modalità 3D con fit e unfolding
./build/analysis/psf_analysis/Release/trace_viewer \
    --x1 108.0 --x2 153.0 \
    --p1x 30.0 --p1y 10.0 --p1z 10.0 \
    --p2x -30.0 --p2y -10.0 --p2z -9.0 \
    --dt 1.0 --fit --unfold \
    --output output/psf_analysis/trace_3d.png
```

| Opzione | Default | Descrizione |
|---|---|---|
| `--x1`, `--x2` | — | **Obbligatori.** Posizioni lenti [mm] |
| `--y0` | — | Distanza radiale sorgente [mm] (modalità 2D) |
| `--p1x/y/z`, `--p2x/y/z` | — | Estremi traccia 3D [mm] (modalità 3D) |
| `--psf` | `output/psf/psf_data.root` | Database PSF |
| `--output` | auto | Path immagine PNG |
| `--dt` | 0.1 | Step traccia [mm] |
| `--L` | 10.0 | Lunghezza traccia [mm] (solo modalità 2D) |
| `--sigma` | 1.0 | Scala delle ellissi di covarianza (in unità σ) |
| `--fit` | off | Esegue il fit ODR e mostra χ²/ndof |
| `--unfold` | off | Applica temporal unfolding prima del fit |

### q\_map

Genera la mappa 2D della funzione di qualità `Q(x1, x2)` su tutte le configurazioni del database PSF, oppure la mappa di copertura geometrica in percentuale (`--coverage`). Individua e marca automaticamente il minimo di Q (o il massimo di copertura) con un marker a stella rosso.
```bash
# Mappa Q (default)
./build/analysis/psf_analysis/Release/q_map \
    --psf    output/psf/psf_data.root \
    --config config/config.json \
    --n-tracks 100 --dt 0.1 \
    --log \
    --tsv output/psf_analysis/q_map.tsv

# Mappa copertura geometrica
./build/analysis/psf_analysis/Release/q_map --coverage \
    --psf output/psf/psf_data.root \
    --output output/psf_analysis/coverage_map.png
```

| Opzione | Default | Descrizione |
|---|---|---|
| `--psf` | `output/psf/psf_data.root` | Database PSF |
| `--config` | `config/config.json` | Parametri griglia lenti |
| `--output` | auto | Immagine PNG di output |
| `--tsv` | (disabilitato) | Esporta valori in formato TSV |
| `--scint-x/y/z` | 60, 20, 20 mm | Dimensioni scintillatore |
| `--n-tracks` | 100 | Tracce casuali per configurazione |
| `--dt` | 0.1 | Step traccia [mm] |
| `--min-hits` | 10.0 | Hit minime per punto PSF valido |
| `--trace-frac` | 0.75 | Soglia tracce valide |
| `--unfold-dz` | 0.0 (auto) | Passo di srotolamento fisso [mm/passo] |
| `--no-unfold` | off | Disabilita il temporal unfolding |
| `--coverage` | off | Modalità mappa copertura geometrica |
| `--log` | off | Scala logaritmica sull'asse Z (colori) |
| `--norm` | off | (non usato internamente, per compatibilità) |

**Output grafico (modalità Q)**: mappa TH2D con palette `kBird`, scala colori ai percentili 0–95, pannello info con parametri e coordinate del minimo.

**Output grafico (modalità coverage)**: mappa TH2D con palette `kViridis`, range [0, 100]%, pannello info con il massimo di copertura.

**Output TSV**: tabella con colonne `x1 \t x2 \t Q \t n_traces \t n_failed \t n_invalid \t config_valid` (modalità Q) oppure `x1 \t x2 \t coverage_pct \t n_y0_evaluated \t config_valid` (modalità coverage).

### chi2\_map

Genera la mappa 2D della linearità della risposta ottica, calcolando per ogni configurazione un fit di piano `μ_{y,z} = a + b·y₀ + c·x₀` sui dati PSF e riportando il χ²/ndof combinato (Y + Z). Complementare a `q_map`: bassa linearità (basso χ²) indica che il sistema si comporta come una lente ideale su tutto il campo.
```bash
./build/analysis/psf_analysis/Release/chi2_map \
    --psf output/psf/psf_data.root \
    --output output/psf_analysis/chi2_map.png \
    --min-hits 10 \
    --log
```

| Opzione | Default | Descrizione |
|---|---|---|
| `--psf` | `output/psf/psf_data.root` | Database PSF |
| `--config` | `config/config.json` | Parametri griglia |
| `--output` | `output/psf_analysis/chi2_map.png` | Immagine PNG |
| `--tsv` | (disabilitato) | Esportazione TSV |
| `--min-hits` | 10.0 | Hit minime per considerare un punto PSF valido |
| `--log` | off | Scala logaritmica |
| `--no-reduced` | off | Usa χ² totale invece di χ²/ndof |

---

## Temporal Unfolding

Il temporal unfolding risolve la degenerazione del fit ODR per configurazioni con forte aberrazione di campo, in cui la traccia media si ripiega su se stessa nel piano (μ_y, μ_z). La trasformazione applica uno spostamento artificiale progressivo lungo la direzione **ortogonale** alla traccia stessa:
```
ũ_i = u_i + i · δS · n̂⊥
```

dove `n̂⊥` è il vettore normale alla traccia calcolato dai punti estremi validi, e `δS` è il passo di srotolamento (auto = L/(N−1) oppure fisso via `z_unfold_step`).

La trasformazione è applicata **esclusivamente** all'interno di `compute_Q` su una copia locale della traccia. Il database PSF, `build_trace` e `trace_viewer` operano sempre sulle coordinate fisiche (μ_y, μ_z).

Le matrici di covarianza Σ_i non vengono modificate perché l'offset è deterministico (non stocastico): `Cov(ũ_i, ũ_j) = Cov(u_i, u_j)`.

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

### Suite T — `fit_trace`

| Test | Scenario | Verifica |
|---|---|---|
| **T1** | Retta `z = y`, cov isotropa | `a=1`, `b=0`, `χ²≈0`, pull ≈ 0, `n_points_used=21` |
| **T2** | Retta `z = 0.5y + 3`, cov isotropa | `a=0.5`, `b=3`, `χ²≈0`, `σ_a > 0`, `n_points_used=31` |
| **T3** | Retta `z = 2y − 1`, cov anisotropa | Parametri invariati rispetto ai pesi |
| **T4** | Outlier singolo `Δz = 0.5 mm` | Pull ≈ `Δz/σ`, χ²/ndof dominato dall'outlier |
| **T5** | Meno di 3 punti validi | Eccezione `std::invalid_argument` |
| **T6** | Covarianza degenere (`var = 0`) | Nessun crash, floor `1e-6 mm²`, convergenza corretta |
| **T7** | Cov non diagonale, outlier | `σ_d ≈ √(cov_zz)` verifica formula `n̂ᵀ Σ n̂` |

### Suite TV — `is_trace_valid`

| Test | Scenario | Verifica |
|---|---|---|
| **TV1** | Tutti `valid=true` | Valida per soglie 0.75 e 1.0 |
| **TV2** | 60% valid | Valida per soglia 50%, invalida per 75% |
| **TV3** | Traccia vuota | `false` |

### Suite TQ — `compute_Q`

| Test | Scenario | Verifica |
|---|---|---|
| **TQ1** | PSF ideale, unfolding OFF | `Q ≈ 0`, `n_traces > 0`, `n_failed = 0` |
| **TQ2** | Scintillatore personalizzato | `n_traces > 0`, `Q ≈ 0` |
| **TQ3** | Config non presente | `std::invalid_argument` |
| **TQ4** | PSF curva, due σ_z diverse | `Q(σ_piccola) > Q(σ_grande)` |
| **TQ5** | Traccia ripiegata vs lineare, unfolding ON/OFF | `χ²(fold) >> χ²(lin)` con ON; `Δχ²(ON) >> Δχ²(OFF)` |

### Note sull'inizializzazione dei dati sintetici

- **`TracePoint::valid = true`** deve essere impostato esplicitamente: la zero-inizializzazione lascia `valid = false`, causando `std::invalid_argument` in `fit_trace`.
- **`PSFPoint::on_detector = true`** deve essere impostato esplicitamente: `build_trace` copia questo flag in `TracePoint::valid`.

---

## Configurazione raccomandata per `q_map`
```bash
./build/analysis/psf_analysis/Release/q_map \
    --psf    output/psf/psf_data.root \
    --n-tracks 100                    \
    --dt 0.1                          \
    --unfold-dz 0.000002              \
    --trace-frac 0.50                 \
    --min-hits 10                     \
    --log
```

| Parametro | Valore | Motivazione |
|---|---|---|
| `--n-tracks 100` | 100 tracce | Compromesso tra stabilità statistica e tempo di calcolo (~22 s per 16 chunk) |
| `--dt 0.1` | Passo 0.1 mm | ~101 punti per traccia di 10 mm; sufficiente per il fit ODR |
| `--unfold-dz 0.000002` | δz = 2×10⁻⁶ mm/passo | Rotto la degenerazione senza distorcere χ² delle tracce lineari (offset totale ≈ 0.0002 mm ≪ σ_PSF) |
| `--trace-frac 0.50` | 50% | Permette di mappare configurazioni borderline senza perdere celle |
| `--min-hits 10` | 10 hit | Soglia conservativa per PSF affidabile |
| `--log` | scala log | Q copre 2–3 ordini di grandezza; la scala log rende leggibile l'intera mappa |

---

## File di configurazione

### `config/config.json`
```json
{
  "x_min": 33.0,
  "x_max": 171.0,
  "dx": 3.0,
  "r1": 38.6,
  "h1": 12.5,
  "r2": 30.9,
  "h2": 16.3,
  "source_width": 10.0,
  "source_height": 5.0,
  "source_thickness": 60.0,
  "n_photons": 10000,
  "lower_percentile": 0.0,
  "upper_percentile": 0.0,
  "source_x_min": -30.0,
  "source_x_max": 30.0,
  "source_dx": 1.0,
  "source_y_min": 0.0,
  "source_y_max": 14.14,
  "source_dy": 1.0
}
```

| Parametro | Descrizione |
|---|---|
| `x_min`, `x_max`, `dx` | Range e passo della scansione lenti [mm] |
| `r1`, `h1` | Raggio e spessore della lente 75mm (solo GDML default) [mm] |
| `r2`, `h2` | Raggio e spessore della lente 60mm (solo GDML default) [mm] |
| `source_width/height/thickness` | Dimensioni sorgente per `optimization` [mm] |
| `n_photons` | Fotoni per run in `optimization` |
| `lower_percentile`, `upper_percentile` | Taglio scala colori per `plot2D` |
| `source_x_min/max`, `source_dx` | Griglia X sorgente per `lens_simulation` [mm] |
| `source_y_min/max`, `source_dy` | Griglia Y sorgente per `lens_simulation` [mm] |

I parametri `config_id_offset` e `run_id_offset` vengono aggiunti dinamicamente da `run.sh` durante la parallelizzazione per garantire ID globalmente univoci.

### `config/config_profile.json`

Configurazione ridotta usata per il profiling rapido (passo `dx = 30 mm`, soglie percentile più strette):
```json
{
  "x_min": 33.0, "x_max": 171.0, "dx": 30.0,
  "r1": 38.6, "h1": 12.5, "r2": 30.9, "h2": 16.3,
  "lower_percentile": 0.45, "upper_percentile": 0.0
}
```

---

## Formato dei file ROOT

I file ROOT utilizzano un formato a **vettori di hit** per massimizzare l'efficienza di archiviazione (riduzione ~40% rispetto al formato riga-per-hit). La compressione è LZ4 livello 4 (`SetCompressionLevel(404)`).

### lens.root (lens\_simulation)
```
lens.root
├── TTree "Configurations"
│   ├── config_id  Int_t
│   ├── x1, x2     Double_t [mm]
│   ├── lens75_id  String
│   └── lens60_id  String
│
└── TTree "Runs"
    ├── run_id      Int_t
    ├── config_id   Int_t
    ├── x_source    Float_t  [mm]
    ├── y_source    Float_t  [mm]
    ├── n_hits      Double_t
    ├── y_hits      vector<float> [mm]
    └── z_hits      vector<float> [mm]
```

### events.root (optimization)
```
events.root
├── TTree "Configurations"
│   ├── config_id  Int_t
│   ├── x1, x2     Double_t
│   ├── lens75_id  String
│   └── lens60_id  String
│
└── TTree "Efficiency"
    ├── config_id  Int_t
    ├── n_photons  Int_t
    └── n_hits     Double_t
```

### psf\_data.root (psf\_extractor)
```
psf_data.root
└── TTree "PSF"
    ├── config_id            Int_t
    ├── x1, x2               Double_t [mm]
    ├── x_source             Float_t  [mm]
    ├── y_source             Float_t  [mm]
    ├── mean_y, mean_z       Double_t [mm]
    ├── cov_yy, cov_yz, cov_zz  Double_t [mm²]
    ├── n_hits_filtered      Double_t
    ├── n_hits_raw           Double_t
    └── on_detector          Bool_t
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
./build/Release/optimization_main -g geometry/main.gdml -b -o
./build/analysis/Release/plot2D
# → output/efficiency2D.png

# ── 4. Beam scan PSF (parallelizzato) ───────────────────────────────────────
./scripts/run.sh lens local --jobs $(nproc --all)
# → output/lens_simulation/lens_<timestamp>.root

# ── 5. Estrazione PSF ────────────────────────────────────────────────────────
./build/analysis/Release/psf_extractor \
    output/lens_simulation/lens_<timestamp>.root \
    output/psf/psf_data.root
# → output/psf/psf_data.root

# ── 6. Mappa linearità (diagnostica veloce, nessun fit traccia) ──────────────
./build/analysis/psf_analysis/Release/chi2_map \
    --psf output/psf/psf_data.root --log
# → output/psf_analysis/chi2_map.png

# ── 7. Mappa Q — ottimizzazione PSF-based ────────────────────────────────────
./build/analysis/psf_analysis/Release/q_map \
    --psf output/psf/psf_data.root \
    --n-tracks 100 --dt 0.1 --unfold-dz 0.000002 \
    --trace-frac 0.50 --log \
    --tsv output/psf_analysis/q_map.tsv
# → output/psf_analysis/q_map.png
# Stampa su stdout: x1*, x2*, Q_min

# ── 8. Mappa copertura geometrica (complementare a Q) ────────────────────────
./build/analysis/psf_analysis/Release/q_map --coverage \
    --psf output/psf/psf_data.root
# → output/psf_analysis/coverage_map.png

# ── 9. Analisi traccia per la configurazione ottimale ───────────────────────
./build/analysis/psf_analysis/Release/trace_viewer \
    --x1 <x1*> --x2 <x2*> --y0 5.0 --fit
# → output/psf_analysis/trace_x1_<x1*>_x2_<x2*>_y0_5.0.png

# ── 10. Diagnostica beam scan ────────────────────────────────────────────────
./build/analysis/Release/beam_scan_plot <x1*> <x2*>
./build/analysis/Release/m_c_creator   <x1*> <x2*> 5.0

# ── 11. Visualizzazione interattiva geometria ────────────────────────────────
./build/Release/lens_simulation_main -g geometry/main.gdml -v
```

---

## Tool: exp2\_main — Analisi immagini FITS laser

**Eseguibile**: `build/analysis/exp2/Release/exp2_main`

**Scopo**: Confronto quantitativo di quattro configurazioni ottiche (good/bad × focus/nofocus) su stack di immagini FITS 16-bit acquisite con un laser. Per ogni configurazione estrae il profilo trasversale della traccia laser, misura la larghezza PSF e la linearità del centroide.

### Pipeline

```
FITS frames (signal + background)
  → sigma-clipping stack (SigmaClip / Mean / Median)
    → differenza signal − background
      → stima angolo traccia (PCA + momento di inerzia)
        → raffinamento angolo iterativo (fino a 2 iterazioni ODR)
          → estrazione profilo a slice (centroide pesato ISO 11146)
            → fit ODR lineare sul centroide
              → sigma_minor / sigma_mean / chi²/ndof
```

### Metriche principali

| Metrica | Descrizione |
|---|---|
| `sigma_minor` | Semiasse minore della distribuzione 2D (invariante per rotazione; confrontabile tra blob e streak) |
| `sigma_mean` | Media pesata inversa-varianza delle larghezze Gaussiane delle slice (valida solo per tracce lineari) |
| `aspect_ratio` | `sigma_major / sigma_minor`; se < `min_aspect_ratio` (default 2.0), il blob è troppo circolare e la traccia non viene estratta |
| `chi²/ndof` | Qualità del fit ODR lineare sul centroide (≈1 per traccia lineare ideale) |

### Algoritmi di estrazione traccia

#### Centroide pesato (Priorità 1 — metodo principale)

Il centroide di ogni slice trasversale viene calcolato con il **primo momento pesato per intensità** (standard ISO 11146 per beam profiler):

```
w(s)       = max(0, I(s) − B)         B = mediana dei bordi della slice
centroid   = Σ(s · w(s)) / Σw
σ_centroid = σ_dist / √N_eff          σ_dist = √(secondo momento)
N_eff      = (Σw)² / Σ(w²)            pixel count effettivo
```

Questo metodo non richiede convergenza, tollera profili troncati ai bordi del FOV e profili non-Gaussiani. Il **fit Gaussiano 1D** (Levenberg-Marquardt) viene mantenuto solo per la misura della larghezza PSF `sigma`, con fallback al secondo momento in caso di non convergenza.

#### Raffinamento iterativo dell'angolo (Priorità 2)

L'angolo iniziale (PCA + momento di inerzia) viene raffinato con un loop di massimo 2 iterazioni:

```
δ (YvsZ axis) = atan(a_ODR)      con  center = a·t + b
δ (ZvsY axis) = atan(1 / a_ODR)  con  t = a·center + b
angle_new = angle_old + δ   (interrotto se |δ| < 0.05° o |δ| > 10°)
```

La deriva lineare del centroide nel fit ODR è geometricamente equivalente all'errore angolare residuo; il pendio `a` fornisce quindi una stima diretta della correzione necessaria.

#### Trimming basato su `center_err` (Priorità 3)

Il segmento valido della traccia viene identificato come il più lungo blocco contiguo di slice con `center_err ≤ trace_trim_max_center_err` (default 5.0 px), invece che con la soglia SNR precedente. Questo preserva slice a basso SNR (PSF larga, bordi FOV) fintanto che la posizione del centroide è stimabile con precisione sufficiente.

### Utilizzo

```bash
./build/analysis/exp2/Release/exp2_main \
  --data-dir /mnt/external_ssd/riptide/exp2/ \
  --no-root \
  --snr-min 1.0 \
  --slice-width 15 --slice-step 1 \
  --center-err-floor 0.5 --center-err-scale 2.0 \
  --sigma-err-floor 0.5 --sigma-err-scale 1.0 \
  --trim-frac 0.10 --trim-min-snr 3.0 \
  --trim-pad-slices 25 --trim-max-center-err 5.0 \
  --min-aspect-ratio 2.0 \
  --x1-good 110 --x2-good 130 \
  --x1-bad 60  --x2-bad 125 \
  --xdet-good-opt 165 --xdet-bad-opt 150
```

| Parametro | Default | Descrizione |
|---|---|---|
| `--min-aspect-ratio` | 2.0 | Soglia aspect_ratio per rilevare tracce analizzabili |
| `--trim-max-center-err` | 5.0 px | Soglia max `center_err` per il trimming |
| `--center-err-scale` | 1.0 | Fattore moltiplicativo per l'incertezza del centroide |
| `--sigma-err-floor` | 0.2 px | Incertezza minima per `sigma` |
| `--trim-pad-slices` | 10 | Slice di padding intorno alla regione valida |

---

## Riferimenti metodologici

### Estrazione traccia e centroide

- **ISO 11146-1:2005** — *Lasers and laser-related equipment: Test methods for laser beam widths, divergence angles and beam propagation ratios*. Definisce il metodo dei momenti del secondo ordine per la misura della larghezza del fascio (usato per `sigma_minor`, `sigma_dist`).

- **Zhang, C. & Couloigner, I. (2007)**. *Accurate Centerline Detection and Line Width Estimation of Thick Lines Using the Radon Transform*. IEEE Transactions on Image Processing, 16(2), 310–316. DOI: [10.1109/TIP.2006.887731](https://doi.org/10.1109/TIP.2006.887731). Riferimento per la robustezza della trasformata di Radon pesata nella localizzazione di linee diffuse; base per il metodo del centroide pesato applicato alle slice.

- **Vyas, A. et al. (2009)**. *Centroid Detection by Gaussian Pattern Matching in Adaptive Optics*. arXiv:0910.3386. Confronto tra fit Gaussiano e centroide pesato per la stima della posizione: il centroide pesato è più robusto a basso SNR e profili troncati.

- **Thomas, S. et al. (2006)**. *Comparison of centroid computation algorithms in a Shack–Hartmann sensor*. Monthly Notices of the Royal Astronomical Society, 371(1), 323–336. DOI: [10.1111/j.1365-2966.2006.10661.x](https://doi.org/10.1111/j.1365-2966.2006.10661.x). Dimostra che il centroide pesato per intensità (WCoG) ha bias sistematico ridotto rispetto al fit Gaussiano quando il profilo è non ideale.

### Rilevazione di streak in immagini astronomiche

- **Nir, G. et al.** *pyradon: Python tools for streak detection in astronomical images using the Fast Radon Transform*. GitHub: [guynir42/pyradon](https://github.com/guynir42/pyradon). Implementazione di riferimento per la trasformata di Radon veloce applicata a immagini con streak diffuse.

- **Yanagisawa, T. et al. (2015)**. *Streak Detection and Analysis Pipeline for Space-debris Optical Images*. ResearchGate. Base per la pipeline di estrazione di features (centroide, larghezza, flusso) da immagini ottiche con streak lineari a basso SNR.

### Fit robusto (RANSAC / ODR)

- **Fischler, M.A. & Bolles, R.C. (1981)**. *Random Sample Consensus: A Paradigm for Model Fitting with Applications to Image Analysis and Automated Cartography*. Communications of the ACM, 24(6), 381–395. Riferimento classico per RANSAC applicato al fit di modelli in presenza di outlier (base concettuale per il loop ODR con rejection in `fit_centroid_line`).
```