# riptide_optimization

Simulazione Monte Carlo Geant4 per l'ottimizzazione del posizionamento di un sistema ottico a doppia lente (UVFS) davanti a un fotocatodo GaAsP. Il progetto fa parte dell'esperimento **RIPTIDE** e comprende quattro programmi di simulazione, una suite di strumenti di analisi (mappe di efficienza/PSF/DoF/risoluzione, selezione Pareto, calibrazione sperimentale FITS), una libreria PSF con fit ODR, test unitari e di regressione, e script di gestione simulazioni/trasferimento CNAF.

---

## Indice

1. [Requisiti](#requisiti)
2. [Struttura del progetto](#struttura-del-progetto)
3. [Setup fisico](#setup-fisico)
4. [Build](#build)
5. [Programmi di simulazione](#programmi-di-simulazione)
   - [optimization\_main](#optimization_main)
   - [lens\_simulation\_main](#lens_simulation_main)
   - [dof\_simulation\_main](#dof_simulation_main)
   - [psf\_dof\_scan\_main](#psf_dof_scan_main)
   - [lens\_cutter\_main](#lens_cutter_main)
6. [Esecuzione Parallela e Dashboard](#esecuzione-parallela-e-dashboard)
7. [Output su SSD esterna](#output-su-ssd-esterna)
8. [Log delle simulazioni](#log-delle-simulazioni)
9. [Script di gestione e CNAF](#script-di-gestione-e-cnaf)
10. [Ottimizzazione autonoma](#ottimizzazione-autonoma)
11. [Simulazione singola coppia di lenti](#simulazione-singola-coppia-di-lenti)
12. [Strumenti di analisi](#strumenti-di-analisi)
    - [plot2D](#plot2d)
    - [plot3D](#plot3d)
    - [beam\_scan\_plot](#beam_scan_plot)
    - [m\_c\_creator](#m_c_creator)
    - [psf\_extractor](#psf_extractor)
    - [dof\_map](#dof_map)
    - [resolution\_map](#resolution_map)
    - [pareto\_selector](#pareto_selector)
    - [exp1\_main](#exp1_main)
    - [exp2\_main](#exp2_main)
    - [exp3\_main](#exp3_main)
13. [Libreria psf\_analysis](#libreria-psf_analysis)
    - [trace\_viewer](#trace_viewer)
    - [q\_map](#q_map)
    - [chi2\_map](#chi2_map)
    - [psf\_gaussianity](#psf_gaussianity)
14. [Temporal Unfolding](#temporal-unfolding)
15. [Importance Sampling](#importance-sampling)
16. [Test unitari](#test-unitari)
    - [Regression testing (pytest)](#regression-testing-pytest)
17. [Bibliografia](#bibliografia)
18. [Configurazione raccomandata per `q_map`](#configurazione-raccomandata-per-q_map)
19. [File di configurazione](#file-di-configurazione)
20. [Formato dei file ROOT](#formato-dei-file-root)
21. [Workflow completo](#workflow-completo)
22. [Riferimenti metodologici](#riferimenti-metodologici)
23. [References — Pipeline di ottimizzazione RIPTIDE (Tecniche A/C/E)](#references--pipeline-di-ottimizzazione-riptide-tecniche-ace)

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
│   ├── autonomous_optimizer.py              # Driver autonomo per sweep parametrico completo (two-phase, Pareto)
│   ├── lens_runner.py                       # Pipeline completa per una singola coppia di lenti (no Pareto)
│   ├── pareto_runner.py                     # Post-processing Pareto multi-peso su simulazione già eseguita
│   ├── simulation_log.csv                   # Registro persistente di tutte le simulazioni eseguite
│   ├── run.sh                               # Lancia lens_simulation, optimization o dof_simulation (locale o SSD)
│   ├── monitor.sh                           # Monitoraggio progresso chunk paralleli (bash)
│   ├── dashboard.py                         # Dashboard Python/rich (prioritario su monitor.sh)
│   ├── sim_manager.py                       # Gestione simulation_log.csv e cartelle output/lens_simulations/ (list/delete/sync)
│   ├── verify_chos.sh                       # Verifica rapida risultati di una coppia già simulata
│   ├── sync_to_cnaf.sh                      # Push del progetto su CNAF (rsync, esclude build/output/ROOT)
│   ├── fetch_from_cnaf.sh                   # Pull output/sweep da CNAF (esclude ROOT di default)
│   ├── fetch_lens_from_cnaf.sh              # Pull output/lens_simulations da CNAF (opzioni --pair/--with-root/--log-only)
│   ├── fetch_any_from_cnaf.sh               # Pull generico da cartella CNAF arbitraria (--filter)
│   ├── cnaf_sbatch.sh                       # Job SLURM per autonomous_optimizer.py su CNAF
│   ├── validate_tecnica_c.py                # Validazione Tecnica C: confronto momenti PSF+DoF a 10k vs 5k fotoni
│   └── validate_tecnica_e.py                # Validazione Tecnica E: confronto PsfDofRuns embedded vs standalone
│
├── programs/
│   ├── optimization_main.cpp                # Scansione efficienza geometrica
│   ├── lens_simulation_main.cpp             # Beam scan PSF
│   ├── dof_simulation_main.cpp              # Scansione profondità di campo
│   ├── psf_dof_scan_main.cpp                # Beam scan PSF + acquisizione raggi DoF
│   └── lens_cutter_main.cpp                 # Gestione lenti Thorlabs
│
├── src/
│   ├── common/                              # Codice condiviso
│   │   ├── physics_list.cpp
│   │   ├── primary_generator_action.cpp
│   │   ├── stepping_action.cpp
│   │   └── importance_sampling.cpp          # Campionamento geometrico preferenziale
│   ├── optimization/                        # Sorgenti optimization
│   ├── lens_simulation/                     # Sorgenti lens_simulation
│   ├── dof_simulation/                      # Sorgenti dof_simulation
│   └── psf_dof_scan/                        # Sorgenti psf_dof_scan
│
├── include/
│   ├── common/
│   │   ├── physics_list.hpp
│   │   ├── primary_generator_action.hpp
│   │   ├── stepping_action.hpp
│   │   ├── importance_sampling.hpp          # API campionamento
│   │   └── focus_map.hpp                    # Caricamento TSV piano di fuoco
│   ├── optimization/
│   ├── lens_simulation/
│   ├── dof_simulation/
│   └── psf_dof_scan/
│
├── analysis/
│   ├── plot2D.cpp                           # Mappa 2D efficienza geometrica
│   ├── plot3D.cpp                           # Visualizzazione 3D multi-lente
│   ├── beam_scan_plot.cpp                   # Posizione fotoni vs posizione sorgente
│   ├── m_c_creator.cpp                      # Istogramma 2D hit su detector
│   ├── psf_extractor.cpp                    # Estrazione media e covarianza PSF
│   ├── lens_ranking.hpp                     # Ranking robusto coppie di lenti (score mean-variance), usato da plot2D/plot3D --ranking-only
│   ├── psf_analysis/                        # Libreria PSF + analisi traccia + ottimizzazione
│   │   ├── psf_interpolator.hpp/cpp         # API: load, interpolate, build_trace, fit_trace, compute_Q
│   │   ├── trace_viewer.cpp                 # Traccia 2D/3D con ellissi di covarianza
│   │   ├── q_map.cpp                        # Mappa 2D Q(x1,x2) o copertura geometrica
│   │   ├── chi2_map.cpp                     # Mappa 2D linearità risposta
│   │   ├── psf_gaussianity.cpp              # Test di gaussianità multivariata nucleo/contaminazione (Mardia, Henze-Zirkler)
│   │   ├── test_fit_trace.cpp               # Test unitari: T1–T7, TV1–TV3, TQ1–TQ5
│   │   ├── test_chi2_diagnostics.cpp        # Test unitari diagnostica χ²/ndof: TD1–TD5, TD_REG, TD_ROT
│   │   └── CMakeLists.txt
│   ├── dof_analysis/                        # Profondità di campo e magnificazione
│   │   ├── dof_map.cpp                      # Piano di fuoco, DoF, M, EE80 da focal.root
│   │   ├── dof_plot.cpp                     # Diagnostica estesa
│   │   ├── magnification_map.cpp            # Solo mappa M(x1,x2)
│   │   └── CMakeLists.txt
│   ├── resolution_analysis/                 # Metriche ottiche al piano di fuoco
│   │   ├── resolution_map.cpp               # EE80, Δy_min, DoF da psf_dof.root
│   │   └── CMakeLists.txt
│   ├── pareto_analysis/                     # Selezione ottimale multi-criterio
│   │   ├── pareto_core.hpp                  # Logica fronte di Pareto e Mtot
│   │   ├── pareto_selector.cpp              # Tool principale
│   │   ├── test_pareto_selector.cpp         # Test unitari
│   │   └── CMakeLists.txt
│   ├── exp1/                                # Analisi frame-by-frame immagini FITS
│   │   ├── exp1_main.cpp
│   │   ├── include/
│   │   └── src/
│   ├── exp2/                                # Analisi tracce laser su FITS 16-bit
│   │   ├── exp2_main.cpp
│   │   ├── include/
│   │   └── src/
│   ├── exp3/                                # Calibrazione lineare + confronto Q_exp/Q_sim + M(r) + preview FITS
│   │   ├── exp3_main.cpp
│   │   ├── include/
│   │   ├── src/
│   │   └── test_exp3_calib.cpp              # Test unitari: C1–C6 su immagini FITS sintetiche
│   ├── exp_common/                          # Moduli FITS condivisi
│   │   ├── stacking.hpp                     # Sigma-clipping, mean, median stack
│   │   └── fits_io.hpp                      # Lettura/scrittura FITS
│   └── CMakeLists.txt
│
├── lens_cutter/
│   ├── include/lens_cutter.hpp
│   ├── src/lens_cutter.cpp
│   └── lens_data/
│       ├── thorlabs_biconvex.tsv
│       └── thorlabs_planoconvex.tsv
│
├── geometry/
│   ├── main.gdml                            # Geometria standard per optimization e lens_simulation
│   ├── dof_geometry.gdml                    # Geometria per dof_simulation e psf_dof_scan
│   ├── define.xml
│   ├── materials.xml
│   ├── solids.xml
│   └── structure.xml
│
├── macros/
├── tests/                                   # Regression testing (pytest) — verifica 10 bug storici
│
├── config/
│   ├── config.json
│   ├── config_profile.json
│   ├── analysis_params.json                 # Parametri di tuning per autonomous_optimizer
│   ├── config_{detailed,intermediate,fast,screening,test,psf_gaussianity,validate_c_10k,validate_c_5k,validate_e}.json
│   │                                         # Varianti velocità/validazione, vedi tabella in File di configurazione
│   └── exp3/                                # Configurazione exp3
│
├── external/
│   └── lyra/
│
└── output/
    ├── sweep/                   # autonomous_optimizer.py — sweep --all-lenses su tutte le coppie
    │   └── <tag>/               #   plots/, optimizer_*.log, report_*.md
    └── lens_simulations/        # lens_runner.py — simulazione di una coppia specifica
        └── {l1_id}_{l2_id}/    #   plots/, data/, pareto/, pipeline.log
```

---

## Setup fisico

Il sistema simulato è composto da tre elementi posizionati lungo l'asse X, all'interno di un volume mondo cubico di 1000×1000×1000 mm³ riempito d'aria:
```
Sorgente fotoni  →  [Lente L1]  →  [Lente L2]  →  [Fotocatodo GaAsP 16×16mm]
      (GPS)           (UVFS)           (UVFS)              (sensore)
         x=0        x ≈ 99–300mm    x ≈ 130–340mm           variabile
```

**Lente L1** (`l1`): ellissoide in UV Fused Silica (UVFS), raggio 38.6 mm, spessore 12.5 mm. Sostituibile con qualsiasi lente Thorlabs tramite `--l1-id`.

**Lente L2** (`l2`): ellissoide UVFS, raggio 30.9 mm, spessore 16.3 mm. Sostituibile tramite `--l2-id`.

**Fotocatodo GaAsP**: lastra quadrata 16×16×0.01 mm, indice di rifrazione 3.5–3.8 nel range 2–4 eV, lunghezza di assorbimento ~1 µm.

**Sorgente**: fotoni ottici a 2.5 eV (≈ 496 nm), generati da GPS Geant4 su una sferetta di raggio 0.1 mm con distribuzione angolare isotropa diretta lungo +X (`mintheta=0`, `maxtheta=90 deg`).

**Processi ottici attivi**: rifrazione, riflessione (Fresnel), assorbimento. Scintillazione, Cherenkov, WLS e WLS2 sono esplicitamente disabilitati in `PhysicsList` perché i fotoni sono generati direttamente dalla GPS.

**SteppingAction**: uccide i fotoni che si allontanano nella direzione −X (x < −1 mm con px < 0) o che escono dalla regione utile (|y| > 150 mm o |z| > 150 mm), evitando tracking inutile.

**Geometria DoF** (`dof_geometry.gdml`): stessa configurazione ottica di `main.gdml`, con l'aggiunta di un piano virtuale di acquisizione a `x_virtual = x2 + dof_x_virtual_offset` (default 30 mm dopo la seconda lente). Usata da `dof_simulation_main` e `psf_dof_scan_main` per catturare posizione e direzione `(y₀, z₀, dy, dz)` di ogni raggio dopo la seconda lente, senza elementi ottici intermedî.

---

## Build
```bash
cmake -S . -B build/ -G "Ninja Multi-Config"
cmake --build build/ --config Release
```

Eseguibili prodotti:
```
build/Release/optimization_main
build/Release/lens_simulation_main
build/Release/dof_simulation_main
build/Release/psf_dof_scan_main
build/Release/lens_cutter_main

build/analysis/Release/plot2D
build/analysis/Release/plot3D
build/analysis/Release/beam_scan_plot
build/analysis/Release/m_c_creator
build/analysis/Release/psf_extractor

build/analysis/psf_analysis/Release/trace_viewer
build/analysis/psf_analysis/Release/q_map
build/analysis/psf_analysis/Release/chi2_map
build/analysis/psf_analysis/Release/psf_gaussianity
build/analysis/psf_analysis/Release/test_fit_trace
build/analysis/psf_analysis/Release/test_chi2_diagnostics

build/analysis/dof_analysis/Release/dof_map
build/analysis/dof_analysis/Release/dof_plot
build/analysis/dof_analysis/Release/magnification_map

build/analysis/resolution_analysis/Release/resolution_map

build/analysis/pareto_analysis/Release/pareto_selector
build/analysis/pareto_analysis/Release/test_pareto_selector

build/analysis/exp1/Release/exp1_main
build/analysis/exp2/Release/exp2_main
build/analysis/exp3/Release/exp3_main
build/analysis/exp3/Release/test_exp3_calib
```

> **Build Release**: aggiunge automaticamente `-march=native -O3 -ffast-math`.

> **Cartelle output**: CMake crea automaticamente tutte le sottocartelle `output/` durante la configurazione.

---

## Programmi di simulazione

### optimization\_main

**Scopo**: scansione grezza dell'efficienza geometrica. Per ogni coppia (x1, x2), conta quanti fotoni emessi da una sorgente rettangolare raggiungono il fotocatodo.

```bash
./build/Release/optimization_main -g geometry/main.gdml -b -o
./build/Release/optimization_main -g geometry/main.gdml -b -o --all-lenses
```

| Flag | Tipo | Descrizione |
|---|---|---|
| `-g`, `--geometry` | path | **Obbligatorio.** Percorso al file GDML |
| `-m`, `--macro` | path | Macro Geant4 (default: `macros/optimization.mac`) |
| `-v`, `--visualize` | flag | Visualizzazione interattiva |
| `-b`, `--batch` | flag | Modalità batch |
| `-o`, `--optimize` | flag | Avvia la scansione |
| `--all-lenses` | flag | Scansione di tutte le combinazioni Thorlabs |
| `--l1-id` | string | ID lente Thorlabs per L1 |
| `--l2-id` | string | ID lente Thorlabs per L2 |
| `--output` | path | File ROOT di output (default: `output/events.root`) |
| `--config` | path | File `config.json` (default: `config/config.json`) |
| `--ssd` | flag | Output su SSD esterna con timestamp |
| `--ssd-mount` | path | Mount point SSD (default: `/mnt/external_ssd`) |

**Output**: `events.root` — TTree `Configurations` + `Efficiency` (vedi [Formato dei file ROOT](#formato-dei-file-root)).

> **`--all-lenses` e cataloghi pre-definiti**: il flag legge a runtime i cataloghi TSV in `lens_cutter/lens_data/` (`thorlabs_biconvex.tsv`, 9 lenti; `thorlabs_planoconvex.tsv`, 7 lenti) e genera il **prodotto cartesiano completo** di tutte le coppie disponibili. Il catalogo è fisso al contenuto dei TSV: non è limitabile a un sottoinsieme senza modificare quei file. La stessa logica vale per `dof_simulation_main`.

---

### lens\_simulation\_main

**Scopo**: simulazione dettagliata del beam scan per la caratterizzazione della PSF. Per ogni configurazione e per ogni coppia (x\_source, y\_source), registra le posizioni (y, z) dei fotoni sul fotocatodo.

```bash
./build/Release/lens_simulation_main -g geometry/main.gdml -b -l
./build/Release/lens_simulation_main -g geometry/main.gdml -b -l \
    --focus-tsv output/dof_analysis/dof_map.tsv
```

| Flag | Tipo | Descrizione |
|---|---|---|
| `-g`, `--geometry` | path | **Obbligatorio.** |
| `-m`, `--macro` | path | Macro (default: `macros/lens_simulation.mac`) |
| `-v`, `--visualize` | flag | Visualizzazione interattiva |
| `-b`, `--batch` | flag | Modalità batch |
| `-l`, `--lens-sim` | flag | Avvia il beam scan |
| `--l1-id` | string | ID lente Thorlabs per L1 |
| `--l2-id` | string | ID lente Thorlabs per L2 |
| `--output` | path | File ROOT (default: `output/lens_simulation/lens.root`) |
| `--config` | path | File `config.json` |
| `--focus-tsv` | path | TSV da `dof_map` con colonne `x1,x2,x_focus`: sposta il detector al piano di fuoco ottimale per ogni (x1,x2), purché `x_focus ≥ x2 + lens_det_gap` |
| `--ssd` | flag | Output su SSD |
| `--ssd-mount` | path | Mount point SSD |

**Griglia di campionamento**: la sorgente scorre su una griglia 2D definita in `config.json` (`source_x_min/max/dx`, `source_y_min/max/dy`). Il limite `source_y_max = 10√2 ≈ 14.14 mm` garantisce copertura completa per tracce di lunghezza 10 mm a qualsiasi y₀ ≤ 10 mm.

**Output**: `lens.root` — TTree `Configurations` + `Runs`.

---

### dof\_simulation\_main

**Scopo**: scansione della profondità di campo. Per ogni coppia (x1, x2), registra al piano virtuale `x_virtual = x2 + dof_x_virtual_offset` la posizione trasversa `(y₀, z₀)` e la direzione `(dy = py/px, dz = pz/px)` di ogni fotone dopo la seconda lente. Questi dati permettono di propagare i raggi linearmente verso qualsiasi piano di rilevazione senza ulteriori simulazioni Geant4.

```bash
./build/Release/dof_simulation_main -g geometry/dof_geometry.gdml -b -d
./build/Release/dof_simulation_main -g geometry/dof_geometry.gdml -b -d \
    --all-lenses --ssd
```

| Flag | Tipo | Descrizione |
|---|---|---|
| `-g`, `--geometry` | path | **Obbligatorio.** (default: `geometry/dof_geometry.gdml`) |
| `-m`, `--macro` | path | Macro Geant4 |
| `-v`, `--visualize` | flag | Visualizzazione interattiva |
| `-b`, `--batch` | flag | Modalità batch |
| `-d`, `--dof` | flag | Avvia la scansione DoF |
| `--all-lenses` | flag | Scansione tutte le combinazioni Thorlabs |
| `--l1-id` | string | ID lente Thorlabs per L1 |
| `--l2-id` | string | ID lente Thorlabs per L2 |
| `--output` | path | File ROOT (default: `output/dof_simulation/focal.root`) |
| `--config` | path | File `config.json` |
| `--ssd` | flag | Output su SSD con timestamp |
| `--ssd-mount` | path | Mount point SSD |

**Parametri config**: `dof_n_photons` (default 100000), `dof_source_halfy/halfx` (semi-estensione sorgente), `dof_source_y_centre`, `dof_x_virtual_offset`, `use_importance_sampling`.

**Output**: `focal.root` — TTree `FocalConfigurations` + `FocalRays` (vedi [Formato dei file ROOT](#formato-dei-file-root)).

---

### psf\_dof\_scan\_main

**Scopo**: beam scan PSF + acquisizione raggi DoF in un'unica run. Per ogni configurazione e per ogni posizione sorgente (x\_src, y\_src) della griglia, scansiona simultaneamente il profilo PSF e registra i raggi al piano virtuale. Produce un unico file `psf_dof.root` usato da `resolution_map` per calcolare le metriche al piano di fuoco.

```bash
./build/Release/psf_dof_scan_main -g geometry/dof_geometry.gdml -b -p
./build/Release/psf_dof_scan_main -g geometry/dof_geometry.gdml -b -p \
    --l1-id LB4592 --l2-id LB4553
```

| Flag | Tipo | Descrizione |
|---|---|---|
| `-g`, `--geometry` | path | **Obbligatorio.** (default: `geometry/dof_geometry.gdml`) |
| `-m`, `--macro` | path | Macro Geant4 |
| `-v`, `--visualize` | flag | Visualizzazione |
| `-b`, `--batch` | flag | Batch mode |
| `-p`, `--psf-dof` | flag | Avvia la modalità PSF+DoF |
| `--l1-id` | string | ID lente L1 |
| `--l2-id` | string | ID lente L2 |
| `--output` | path | File ROOT (default: `output/psf_dof_simulation/psf_dof.root`) |
| `--config` | path | File `config.json` |
| `--ssd` | flag | Output su SSD |
| `--ssd-mount` | path | Mount point SSD |

**Parametri config**: `psf_dof_n_photons` (default 10000), `psf_dof_x_virtual_offset`, `psf_dof_save_hits`.

---

### lens\_cutter\_main

**Scopo**: gestione e generazione di geometrie per lenti commerciali Thorlabs. Permette di elencare le lenti disponibili e generare snippet GDML o solidi Geant4.

```bash
./build/Release/lens_cutter_main --list
./build/Release/lens_cutter_main --id LB4592
./build/Release/lens_cutter_main --list \
    --catalog lens_cutter/lens_data/thorlabs_planoconvex.tsv
```

Il tool legge le specifiche da file TSV nel formato Thorlabs:

- **Biconvesse** (`thorlabs_biconvex.tsv`): `Item # | Diameter | Focal Length | Radius of Curvature | Center Thickness | Edge Thickness | Back Focal Length`. Il tipo è rilevato automaticamente.
- **Plano-convesse** (`thorlabs_planoconvex.tsv`): schema identico con colonna aggiuntiva `Rotation_deg` [deg]. `0` = lato curvo verso la sorgente, `180` = lato piano verso la sorgente. La rotazione è propagata automaticamente al `G4VPhysicalVolume`.

Le simulazioni `optimization` e `lens_simulation` accettano ID da entrambi i cataloghi tramite `--l1-id` e `--l2-id`; la geometria (solido + rotazione + offset centro geometrico) è impostata automaticamente in base al tipo rilevato.

---

## Esecuzione Parallela e Dashboard

Il framework supporta la parallelizzazione tramite `scripts/run.sh`, che divide il carico di lavoro in chunk indipendenti eseguiti in parallelo e li fonde con `hadd` al termine.

```bash
./scripts/run.sh lens local --jobs $(nproc --all)
./scripts/run.sh lens ssd   --jobs 8
./scripts/run.sh opt  local --jobs 4
./scripts/run.sh opt  local --jobs 4 --all-lenses
./scripts/run.sh dof  local --jobs $(nproc --all)
```

Lo script tenta prima di avviare `scripts/dashboard.py` (Python + `rich`); se non disponibile, ricade su `scripts/monitor.sh` (bash puro). La dashboard mostra progresso globale e per-chunk, tempo trascorso, ETA e stato di ogni job.

**Chunking automatico**: uno script Python embedded calcola tutte le coppie (x1, x2) valide (rispettando i vincoli di non-collisione geometrica) e le distribuisce bilanciata in N chunk. Ogni chunk riceve `config_id_offset` e `run_id_offset` univoci per garantire ID globalmente monotoni dopo il merge.

**Modalità super-batch** (invocata da `lens_runner.py --disk-guard-gb`): se il config chunk JSON contiene la chiave `pairs`, `run.sh` usa quelle coppie pre-computate invece di ricalcolarle, e legge `config_id_offset_base`/`run_id_offset_base` per calcolare gli offset corretti tra batch sequenziali. Questo permette a `lens_runner.py` di eseguire la simulazione in blocchi reattivi al disco senza perdere la continuità degli ID.

**Cleanup**: `Ctrl+C` termina tutti i processi Geant4 in background, rimuove i config temporanei in `/tmp/riptide_chunks_*` e i file ROOT parziali.

---

## Output su SSD esterna

```bash
lsblk
sudo mkdir -p /mnt/external_ssd
sudo mount -o noatime,nodiratime,discard /dev/nvme1n1p2 /mnt/external_ssd
mountpoint -q /mnt/external_ssd && echo "OK"

./build/Release/lens_simulation_main -g geometry/main.gdml -b -l --ssd
./build/Release/dof_simulation_main  -g geometry/dof_geometry.gdml -b -d --ssd
```

Il path di output viene generato automaticamente con timestamp:
```
/mnt/external_ssd/riptide/runs/run_20260315_094512/lens.root
/mnt/external_ssd/riptide/runs/run_20260315_094512/focal.root
```

---

## Log delle simulazioni

**Simulazioni Geant4 via `run.sh` (modalità parallela):** spdlog scrive esclusivamente su stdout/stderr. Durante l'esecuzione i log per-chunk vengono catturati in:
```
/tmp/riptide_chunks_<timestamp>/chunk_N.log
```
Questi file vengono rimossi al termine della run (o su Ctrl+C insieme ai ROOT parziali). **Non esiste un log persistente per simulazioni già concluse**: l'artefatto durevole è il file ROOT prodotto (`lens.root`, `events.root`, `focal.root`, `psf_dof.root`).

**Autonomous optimizer:** crea un log strutturato persistente per ogni variante sweep in:
```
output/sweep/<tag>/optimizer_<timestamp>.log
```
Il log include stdout/stderr di ogni step della pipeline, con timestamp e livello spdlog.

**Registro CSV persistente:** `scripts/simulation_log.csv` traccia tutte le simulazioni eseguite tramite `lens_runner.py`, con colonne:

| Colonna | Descrizione |
|---------|-------------|
| `date` | Data e ora di esecuzione |
| `l1_id`, `l2_id` | ID Thorlabs delle due lenti |
| `detector_mode` | `fixed` (detector a `x_det`) o `mobile` (piano focale) |
| `x_min`, `x_max`, `dx` | Range griglia di scansione |
| `n_photons`, `lens_n_photons`, `dof_n_photons`, `psf_dof_n_photons` | Fotoni per step |
| `lens_gap_margin` | Margine di sicurezza collisione lenti (mm) |
| `storage_type` | `local` o `ssd` |
| `storage_path` | Path dei file ROOT pesanti |
| `output_dir` | Cartella di output grafici e TSV |
| `notes` | Note libere (aggiungibili con `--notes "..."`) |

Il file è pre-popolato con le simulazioni storiche estratte da `exp.txt`. Le nuove simulazioni vengono aggiunte automaticamente al termine di ogni run `lens_runner.py`.

---

## Script di gestione e CNAF

### Gestione simulazioni locali

`scripts/sim_manager.py` opera su `simulation_log.csv` e sulle cartelle `output/lens_simulations/`:

```bash
python3 scripts/sim_manager.py list [--json]                    # elenca le coppie simulate
python3 scripts/sim_manager.py delete <l1_id> <l2_id> [--yes]   # rimuove una coppia (log + output)
python3 scripts/sim_manager.py sync [--dry-run]                 # riallinea log e cartelle su disco
```

### Trasferimento CNAF

Il cluster remoto CNAF (`cnaf-ntof`, path remoto `/home/NTOF/mesini/riptide_optimization`) è usato per eseguire scan pesanti via SLURM. Script di trasferimento:

| Script | Direzione | Note |
|---|---|---|
| `sync_to_cnaf.sh` | locale → CNAF | rsync del progetto, esclude `build/`, `output/`, `.git/`, `*.root` |
| `fetch_from_cnaf.sh` | CNAF → locale | pull di `output/sweep`, esclude i ROOT di default (`--with-root` per includerli) |
| `fetch_lens_from_cnaf.sh` | CNAF → locale | pull di `output/lens_simulations`; opzioni `--with-root`, `--pair L1ID_L2ID`, `--dry-run`, `--log-only` |
| `fetch_any_from_cnaf.sh` | CNAF → locale | pull generico: `bash scripts/fetch_any_from_cnaf.sh <remote_dir> [local_dest] [opzioni]`, opzioni `--with-root`, `--filter GLOB` (ripetibile), `--dry-run` |
| `cnaf_sbatch.sh` | — | job SLURM per `autonomous_optimizer.py` (`--time=52:00:00 --cpus-per-task=16 --mem=32G`); adattare partition/cpus/mem/module prima del submit |

### Validazione e diagnostica

- `validate_tecnica_c.py` — confronta i momenti PSF+DoF calcolati a `lens_n_photons=10000` vs `5000`; approva la riduzione a 5k fotoni (speedup ×2) se l'errore relativo medio su `sigma_y`, `sigma_dz`, `cov_z_dz` è sotto soglia (default 5%):
  ```bash
  python3 scripts/validate_tecnica_c.py psf_dof_10k.root psf_dof_5k.root [--threshold 0.05] [--csv report.csv]
  ```
- `validate_tecnica_e.py` — confronta spot-per-spot i `PsfDofRuns` prodotti embedded (`lens_simulation_main`) vs standalone (`psf_dof_scan_main`):
  ```bash
  python3 scripts/validate_tecnica_e.py embedded_lens.root standalone_psf_dof.root [--threshold 3] [--csv report.csv]
  ```
- `verify_chos.sh` — verifica rapida dei risultati per una coppia già simulata (default `L1=CHOS60R`, `L2=CHOS25R`); richiede che `simulation_summary.json` esista già (eseguire prima `lens_runner.py`).

---

## Ottimizzazione autonoma

`scripts/autonomous_optimizer.py` è il driver autonomo per lo sweep parametrico completo del design space RIPTIDE. Esegue la pipeline di analisi in modalità **two-phase** su un insieme di varianti geometriche, con recovery da fallimento, validazione post-run e ranking Pareto finale.

### Varianti sweep

Di default l'ottimizzatore genera **9 varianti** (3 geometrie × 3 valori di margine):

| Geometria | `source_dx` [mm] | Note |
|-----------|-----------------|------|
| `nominal` | 0.5 | Griglia fine standard |
| `coarse` | 1.0 | Griglia rada |
| `extended` | 1.0 | Range sorgente esteso |

| `lens_gap_margin` [mm] | 0.5 | 1.0 | 2.0 |

### Strategia two-phase

**Fase 1 (fast):** tutte le 9 varianti con `x_max = 150 mm`. Al termine, ogni variante viene classificata per Mtot medio (colonna `Mtot` di `pareto_results.tsv`).

**Fase 2 (full):** le top-k varianti (default k = 3) vengono ri-eseguite con `x_max` pieno. In modalità fuoco fisso, se `x_max > x_det - lens_det_gap`, `x_max` viene clampato al limite geometrico. Dopo ogni run viene applicata una validazione di 3 invarianti; in caso di fallimento viene eseguito un auto-fix (ricompilazione selettiva del target CMake coinvolto) e la run viene ripetuta fino a 3 volte.

### Modalità detector

**Fuoco fisso (default):** il rivelatore rimane a `x_det` per tutte le configurazioni. Adatto a sistemi in cui la distanza detector–ottica è fissa.

**Fuoco mobile (`--mobile-focus`):** il rivelatore si sposta sul piano focale di ogni coppia (x₁, x₂). Richiede `--all-lenses` con `--l1-id`/`--l2-id` come lenti di riferimento per la stima thin-lens iniziale. La pipeline di screening è:

1. Stima thin-lens del piano focale per ogni (x₁, x₂) → `thin_focus.tsv`
2. `opt --all-lenses` con `thin_focus.tsv` → classifica preliminare
3. `plot2D --ranking-only` → CSV ranking senza immagini
4. `dof` per le top-N coppie → `focus_accurate.tsv` via simulazione reale
5. `opt` con `focus_accurate.tsv` (fotoni da `--refine-photons` o `config.json`) → classifica precisa con immagini

### Pipeline per variante

| Step | Binario | Input → Output | Note |
|------|---------|----------------|------|
| 1 | `run.sh opt` | config → `events.root` + `plots/efficiency2D.png` | |
| 2 | `run.sh lens` | config → `lens.root` | con `--focus-tsv` se fuoco mobile |
| 3 | `psf_extractor` | `lens.root` → `psf_data.root` | min_hits ≥ 50 |
| 4 | `q_map` | `psf_data.root` → `q_map.tsv` | |
| 5 | `chi2_map` | `psf_data.root` → `chi2_map.tsv` | |
| 6 | `run.sh dof` | config → `focal.root` | solo fuoco fisso; fuoco mobile usa screening |
| 7 | `dof_map` | `focal.root` → `dof_map.tsv` | solo fuoco fisso |
| 8 | `run.sh psf-dof` | config → `psf_dof.root` | |
| 9 | `resolution_map` | `psf_dof.root` → `resolution_map.tsv` | |
| 10 | `pareto_selector` | tutti i TSV → `pareto_results.tsv` | |

I parametri di ogni step di analisi (min_hits, n_tracks, pesi Pareto, ecc.) sono letti da `config/analysis_params.json` (vedi [File di configurazione](#file-di-configurazione)).

### Recovery

Un file registry JSON in `output/sweep/` registra i path dei file ROOT prodotti. Se il processo viene interrotto e riavviato, i passi già completati vengono saltati.

### Output

```
output/sweep/<tag>/optimizer_<timestamp>.log    # log strutturato per variante
output/sweep/<tag>/plots/efficiency2D.png       # mappa efficienza per variante
output/sweep/lens_screening/plots/              # plot screening (fuoco mobile)
output/sweep/report_<timestamp>.md              # report Markdown finale con top-10 Pareto
```

> `output/sweep/` è la cartella usata dallo sweep `--all-lenses`; per una coppia singola vedi [`lens_runner.py`](#simulazione-singola-coppia-di-lenti) → `output/lens_simulations/`.

### Utilizzo

```bash
python3 scripts/autonomous_optimizer.py [opzioni]
```

| Flag | Default | Descrizione |
|------|---------|-------------|
| `--fast` | off | Solo fase 1 (tutte le varianti, `x_max=150 mm`) |
| `--no-two-phase` | off | Disabilita fase 2 |
| `--top-k N` | 3 | Varianti promosse a fase 2 |
| `--geom` | tutte | Filtra geometria: `nominal`/`coarse`/`extended` |
| `--margin` | tutti | Filtra margine: `0.5`/`1.0`/`2.0` |
| `--l1-id`, `--l2-id` | — | ID lenti specifici |
| `--all-lenses` | off | Fase 0: screening su tutte le coppie del catalogo |
| `--lens-subset` | — | IDs comma-separated per limitare lo screening |
| `--top-n-lenses` | 3 | Top-N coppie selezionate dallo screening |
| `--screening-photons` | 1000 | Fotoni per la fase di screening |
| `--mobile-focus` | off | Fuoco mobile: thin-lens + DOF accurato per posizionare il detector |
| `--refine-photons` | 0 | Fotoni per il re-run opt con fuoco accurato (0 = usa `config.json`) |
| `--local` | off | Output in `output/sweep/` locale anziché SSD esterna |
| `--ssd-mount` | `/mnt/external_ssd` | Mount point SSD |
| `--analysis-params` | `config/analysis_params.json` | Override file parametri analisi |
| `--max-hours` | 48 | Budget tempo totale |
| `--jobs` | nproc | Job Geant4 paralleli per step |
| `--dry-run` | off | Stampa comandi senza eseguirli |
| `--no-commit` | off | Disabilita git commit automatici dei risultati |

---

## Simulazione singola coppia di lenti

`scripts/lens_runner.py` esegue la pipeline completa per una coppia di lenti specificata, senza screening iniziale né Pareto finale. È il modo più diretto per analizzare una coppia specifica già nota.

> Per lo sweep su tutte le coppie del catalogo (`--all-lenses`) usa invece `autonomous_optimizer.py` → `output/sweep/`.

### Output

```
output/lens_simulations/{l1_id}_{l2_id}/
  plots/                    ← tutti i grafici PNG
    efficiency2D.png
    q_map.png  chi2_map.png  coverage_map.png
    dof_focus_map.png  dof_dof_map.png  dof_M_map.png  dof_EE80_map.png
    dof_stripe_map.png  dof_M_error_map.png  dof_within_photocathode.png
    resolution_focus_map.png  resolution_dof_map.png  ...
  data/                     ← TSV + ROOT di analisi
    q_map.tsv  chi2_map.tsv  dof_map.tsv  resolution_map.tsv
    psf_data.root  psf_dof_data.root
    events.root → (symlink a SSD)   lens.root → (symlink a SSD)
    focal.root  → (symlink a SSD)
  pipeline.log              ← log strutturato di tutti gli step
  simulation_summary.json   ← tutti i parametri usati nella run
  focus_accurate.tsv        ← (solo --mobile-focus) piano focale da simulazione DOF
  pareto/                   ← output di pareto_runner.py
    w_M_dominant/
      pareto_results.tsv  pareto_plot.png
    balanced/
      pareto_results.tsv  pareto_plot.png
```

I file ROOT pesanti (`events.root`, `lens.root`, `focal.root`, `psf_dof.root`) vengono salvati su SSD esterna (default) o in locale con `--local`. Un symlink in `data/` mantiene il riferimento senza duplicare i dati.

### Utilizzo

```bash
# Run completo su SSD esterna (ROOT pesanti → /mnt/external_ssd)
python3 scripts/lens_runner.py --l1-id LA4464 --l2-id LA4464R

# Run locale con detector mobile sul piano focale
python3 scripts/lens_runner.py --l1-id LA4464 --l2-id LA4464R --local --mobile-focus

# Dry-run: mostra tutti i comandi senza eseguire nulla
python3 scripts/lens_runner.py --l1-id LA4464 --l2-id LA4464R --local --dry-run

# Resume: salta opt e lens già completati
python3 scripts/lens_runner.py --l1-id LA4464 --l2-id LA4464R --local --skip-opt --skip-lens
```

### Flag principali

| Flag | Default | Descrizione |
|------|---------|-------------|
| `--l1-id`, `--l2-id` | — | ID lenti Thorlabs (**obbligatori**) |
| `--mobile-focus` | off | Fuoco mobile: esegue dof+dof_map prima della pipeline per ottenere il piano focale reale |
| `--ssd-mount PATH` | `/mnt/external_ssd` | Mountpoint SSD per i file ROOT pesanti (es. `/mnt/external_ssd`) |
| `--local` | off | Salva ROOT files in locale (no SSD mount) |
| `--config PATH` | `config/config.json` | Override configurazione |
| `--analysis-params PATH` | `config/analysis_params.json` | Override parametri analisi |
| `--jobs N` | nproc | Job paralleli per run.sh |
| `--dry-run` | off | Stampa comandi senza eseguire |
| `--keep-lens` | off | Non eliminare `lens.root` dopo `psf_extractor` |
| `--keep-psf-dof` | off | Non eliminare `psf_dof.root` dopo `psf_dof_extractor` |
| `--skip-opt` | off | Salta opt (richiede `data/events.root` esistente) |
| `--skip-lens` | off | Salta lens+psf_extractor (richiede `data/psf_data.root`) |
| `--skip-dof` | off | Salta dof+dof_map (richiede `data/dof_map.tsv`) |
| `--skip-psf-dof` | off | Salta psf-dof+psf_dof_extractor (richiede `data/psf_dof_data.root`) |
| `--disk-guard-gb N` | 0 (off) | Super-batch reattivo: esegue `psf_extractor` e cancella `lens.root` a batch, mantenendo almeno N GB liberi su disco. Raccomandato su CNAF con `--disk-guard-gb 400` |
| `--notes "..."` | — | Note aggiuntive nel CSV di log |

### Modalità detector

**Fisso (default):** il detector rimane a `x_det` da `config.json` per tutte le configurazioni.

**Mobile (`--mobile-focus`):** esegue `dof`+`dof_map` come pre-pipeline prima di `opt`, ricavando il piano focale reale per ogni (x₁, x₂) via simulazione Geant4 completa (file `focus_accurate.tsv`). Il risultato viene passato come `--focus-tsv` agli step `opt`, `lens` e `psf-dof`. I passi `dof`/`dof_map` normali (step 7–8) vengono saltati automaticamente. Adatto a sistemi in cui il detector può essere posizionato sul piano focale.

### Post-processing Pareto (`pareto_runner.py`)

Dopo aver completato una run con `lens_runner.py`, è possibile eseguire il Pareto selector con diversi set di pesi senza rieseguire le simulazioni:

```bash
# Configura config/pareto_runs.json con l1_id, l2_id e i set di pesi desiderati
# poi lancia:
python3 scripts/pareto_runner.py

# Override della coppia di lenti:
python3 scripts/pareto_runner.py --l1-id LA4078 --l2-id LA4464R
```

Formato di `config/pareto_runs.json`:

```json
{
  "l1_id": "LA4464",
  "l2_id": "LA4464R",
  "data_dir": null,
  "runs": [
    { "label": "w_M_dominant", "ee80_max": 10.0, "w_eta": 0.1, "w_Q": 0.1, "w_dof": 0.3, "w_M": 0.5 },
    { "label": "balanced",     "ee80_max": 10.0, "w_eta": 0.25, "w_Q": 0.25, "w_dof": 0.25, "w_M": 0.25 }
  ]
}
```

`data_dir` (opzionale): path alternativo ai file di input se diverso da `output/lens_simulations/{l1_id}_{l2_id}/data/`.

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
| `-c`, `--config` | `config/config.json` | Configurazione |
| `-o`, `--output` | `output/efficiency2D.png` | PNG di output |
| `--lens1`, `--lens2` | prima coppia disponibile | Selezione coppia di lenti |
| `--low`, `--high` | da `config.json` | Percentili da escludere dalla scala colori |

---

### plot3D

Visualizzazione 3D multi-lente per la scansione `--all-lenses`. Mostra un canvas con un plot 3D per ogni modello di prima lente: asse X = modello seconda lente, asse Y = x₁, asse Z = x₂, colore = efficienza. Una linea rossa traccia il valor medio pesato (x̄₁, x̄₂) per ogni seconda lente.
```bash
./build/analysis/Release/plot3D -i output/events.root --low 0.1 --high 0.1
./build/analysis/Release/plot3D -i output/events.root --lens1 LB4553
```

---

### beam\_scan\_plot

Per una coppia (x1, x2) fissata, calcola posizione media e deviazione standard dei fotoni al variare della posizione sorgente (x\_src, y\_src). Applica filtro outlier 2σ iterativo (2 iterazioni). Produce un grafico 3D con piano di fit sovrapposto.
```bash
./build/analysis/Release/beam_scan_plot 94.9 186.4
./build/analysis/Release/beam_scan_plot 94.9 186.4 output/lens_simulation/lens.root -30,0,30
# → output/lens_simulation/beam_scan_3D_x1_94.90_x2_186.40.png
```

---

### m\_c\_creator

Istogramma 2D delle hit per una terna (x1, x2, y0) fissata, con filtro ellittico iterativo a 2σ (4 iterazioni).
```bash
./build/analysis/Release/m_c_creator 94.9 186.4 5.0
# → output/mean_covariance_maps/detector_hits_config_<id>_y0_5.0.png
```

---

### psf\_extractor

Estrae media bidimensionale e matrice di covarianza della PSF per tutti i run. Applica filtro outlier ellittico (distanza di Mahalanobis) a 3σ con 4 iterazioni. Marca ogni run con il flag `on_detector` (true se `n_hits_filtered >= min_hits`, default 50).
```bash
./build/analysis/Release/psf_extractor \
    output/lens_simulation/lens.root \
    output/psf/psf_data.root \
    50
# → output/psf/psf_data.root (TTree "PSF")
```

---

### dof\_map

Elabora il file `focal.root` prodotto da `dof_simulation_main` e calcola per ogni configurazione (x1, x2):

- **Piano di fuoco** `x_focus`: posizione che minimizza la deviazione standard trasversa σ_z(x) dei raggi propagati linearmente. Affinamento sub-step con interpolazione parabolica del vertice di σ_z²(x) (tre punti).
- **Profondità di campo** `DoF`: range assiale dove σ_z(x) ≤ k·σ_z_min (default k = √2).
- **Magnificazione** `M`: pendenza della regressione lineare pesata di `y_focus` (posizione trasversa al piano di fuoco) rispetto a `y_source_hits` (posizione trasversa alla sorgente), per fotone. Non un rapporto di deviazioni standard: quella formulazione confonderebbe la magnificazione geometrica con lo spot-blur da aberrazione, sovrastimandola sempre (vedi [Riferimenti](#riferimenti-metodologici)).
- **EE80**: diametro che raccoglie l'80% dell'energia (vedi [Riferimenti](#riferimenti-metodologici)).
- **Larghezza striscia** `stripe_width`: estensione del fascio sul fotocatodo 16 mm.
- **Flag** `focus_before_lens2`: fuoco fisicamente non raggiungibile (x_focus < x2); i bin corrispondenti sono delimitati da una linea di contorno tratteggiata arancione nelle mappe PNG e sintetizzati in una mappa dedicata `dof_focus_zone_map.png` (palette blu→rosso).

**Algoritmo propagazione**: `z(xᵢ) = z₀ + dz·(xᵢ − x_virtual)` — valida nell'ottica geometrica perché nel tratto x_virtual → detector non esistono elementi ottici. Le direzioni `(dy, dz)` codificano già la rifrazione calcolata da Geant4, aberrazioni di ordine superiore incluse.

Il range di scan è limitato inferiormente a `max(scan_min, x2 + lens_det_gap)` per rispettare il vincolo fisico lente→detector.

```bash
./build/analysis/dof_analysis/Release/dof_map \
    --input  output/dof_simulation/focal.root \
    --config config/config.json \
    --k      1.414 \
    --tsv    output/dof_analysis/dof_map.tsv
# → output/dof_analysis/dof_map.png   (quattro heatmap: x_focus, DoF, M, margine fotocatodo)
# → output/dof_analysis/dof_map.tsv
```

| Opzione | Default | Descrizione |
|---|---|---|
| `-i`, `--input` | `output/dof_simulation/focal.root` | File ROOT di input |
| `-c`, `--config` | `config/config.json` | Configurazione |
| `-o`, `--output` | `output/dof_analysis/` | Cartella di output |
| `--scan-min` | `dof_x_scan_min` da config | Inizio scan x_det [mm] |
| `--scan-max` | `dof_x_scan_max` da config | Fine scan x_det [mm] |
| `--scan-step` | `dof_x_scan_step` da config | Passo scan [mm] |
| `--k` | `dof_k_threshold` da config | Soglia DoF: k·σ_z_min |
| `--core-fraction` | 1.0 | Frazione core PSF selezionata tramite distanza di Mahalanobis |
| `--m-target` | `m_target` da config | Magnificazione target per calcolo M_abs_err |
| `--tsv` | (disabilitato) | Esporta risultati in formato TSV |

**Colonne TSV**: `x1 x2 x_focus x_focus_scan dof M m_target M_abs_err EE80 stripe_width within_photocathode n_rays n_rays_core core_fraction config_id focus_before_lens2`

---

### resolution\_map

Calcola le metriche ottiche standardizzate **al piano di fuoco** da `psf_dof.root` (prodotto da `psf_dof_scan_main`):

- **EE80**: `EE80 = 2·1.7941·σ_rms` dove `σ_rms = √((σ_y² + σ_z²)/2)` al fuoco. Approssimazione gaussiana 2D isotropa (vedi [Riferimenti](#riferimenti-metodologici)).
- **Δy\_min**: `Δy_min = k·σ_z_min / |M|` — errore trasversale minimo sull'asse Y della traccia ricondotto alle coordinate dello scintillatore.
- **DoF**: profondità di campo al piano di fuoco (stessa definizione di `dof_map`).

```bash
./build/analysis/resolution_analysis/Release/resolution_map \
    --input  output/psf_dof_simulation/psf_dof.root \
    --config config/config.json \
    --k      1.414 \
    --tsv    output/resolution_analysis/resolution_map.tsv
# → output/resolution_analysis/resolution_DoF_mean_map.png
# → output/resolution_analysis/resolution_EE80_mean_map.png
# → output/resolution_analysis/resolution_map.tsv
```

| Opzione | Default | Descrizione |
|---|---|---|
| `-i`, `--input` | `output/psf_dof_simulation/psf_dof.root` | File ROOT di input |
| `-c`, `--config` | `config/config.json` | Configurazione |
| `-o`, `--output` | `output/resolution_analysis/` | Cartella di output |
| `--k` | `dof_k_threshold` da config | Costante per calcolo DoF |
| `--scan-min`, `--scan-max`, `--scan-step` | da config | Range scansione x_det [mm] |
| `--low`, `--high` | 0.0 | Percentili per scala colori PNG |
| `--tsv` | (disabilitato) | Esportazione TSV |
| `--dump-csv` | (disabilitato) | Dump dati raw |
| `--max-entries`, `--entry-stride`, `--entry-offset` | — | Limitatori dataset (debug) |

**Colonne TSV**: `x1 x2 dof_mean EE80_mean config_id n_runs_dof n_runs_EE80`

---

### pareto\_selector

Seleziona la configurazione ottica ottimale aggregando i risultati di tutti gli strumenti di analisi tramite analisi del **fronte di Pareto** e una **funzione di merito scalare pesata**.

```bash
./build/analysis/pareto_analysis/Release/pareto_selector \
    --events   output/optimization/events.root \
    --qmap     output/psf_analysis/q_map.tsv \
    --chi2map  output/psf_analysis/chi2_map.tsv \
    --dofmap   output/dof_analysis/dof_map.tsv \
    --resolution output/resolution_analysis/resolution_map.tsv \
    --output   output/pareto_analysis/pareto_plot.png \
    --tsv      output/pareto_analysis/pareto_results.tsv \
    --eta-frac 0.75 --focus-tol 15.0 \
    --w-eta 0.35 --w-Q 0.40 --w-dof 0.15 --w-M 0.10
```

**Pipeline interna:**
1. Join inner di tutti i file su `(x1, x2)` con tolleranza 1e-3 mm
2. Filtri hard: `η ≥ eta_frac·η_max`, `|x_focus − x_det| ≤ focus_tol`, `DoF ≥ dof_min`, `EE80 ≤ ee80_max` (opzionale)
3. Calcolo `Mtot` pesato e normalizzato internamente in [0, 1]:
   ```
   Mtot = w_η·(η/η_max) + w_Q·(1−Q/Q_max) + w_dof·(DoF/DoF_max) + w_M·(1−M_abs_err/M_abs_err_max)
   ```
   I pesi sono normalizzati automaticamente se non sommano a 1 (con warning su stderr).
4. Fronte di Pareto su `(η, Q)`: A domina B se `A.η ≥ B.η AND A.Q ≤ B.Q` con almeno una disuguaglianza stretta
5. Ranking dei punti sul fronte per `Mtot` decrescente

| Opzione | Default | Descrizione |
|---|---|---|
| `--events` | `output/optimization/events.root` | Efficienza geometrica |
| `--qmap` | `output/psf_analysis/q_map.tsv` | Qualità Q |
| `--chi2map` | `output/psf_analysis/chi2_map.tsv` | Linearità χ² |
| `--dofmap` | `output/dof_analysis/dof_map.tsv` | DoF e magnificazione |
| `--resolution` | (disabilitato) | EE80 da resolution_map (opzionale) |
| `--output` | `output/pareto_analysis/pareto_plot.png` | PNG output |
| `--tsv` | `output/pareto_analysis/pareto_results.tsv` | Risultati TSV |
| `--eta-frac` | `0.75` | Soglia η relativa |
| `--x-det` | `180.0` | Posizione detector nominale [mm] |
| `--focus-tol` | `15.0` | Tolleranza sul piano di fuoco [mm] |
| `--dof-min` | `0.0` | DoF minima [mm] (0 = disabilitato) |
| `--ee80-max` | `0.0` | EE80 massimo [mm] (0 = disabilitato) |
| `--w-eta` | `0.35` | Peso efficienza in Mtot |
| `--w-Q` | `0.40` | Peso qualità in Mtot |
| `--w-dof` | `0.15` | Peso DoF in Mtot |
| `--w-M` | `0.10` | Peso magnificazione in Mtot |
| `--l1-id`, `--l2-id` | `""` | Filtra per ID lente |

**Output PNG**: scatter plot η/η_max vs Q_max/Q, punti colorati per |M−m\_target|, dimensione ∝ DoF. Fronte di Pareto con bordo rosso; configurazione raccomandata con stella rossa. Tabella top-5 nel pad inferiore.

**Colonne TSV**: `x1 x2 eta eta_norm Q chi2 DoF M M_abs_err x_focus EE80 on_pareto Mtot pareto_rank`

**Test unitari:**
```bash
./build/analysis/pareto_analysis/Release/test_pareto_selector
cd build && ctest -R pareto_selector_unit -V
```

#### Modalità `--weight-sweep` (mappa dei pesi)

Invece del singolo peso `(w_eta, w_Q, w_dof, w_M)`, calcola quale configurazione vince (Mtot massimo) su una **griglia densa di combinazioni di pesi** (vedi [Bibliografia](#bibliografia) per la tecnica). Il join, i pre-filtri (`--eta-frac`, `--focus-tol`, `--dof-min`, `--ee80-max`) e il fronte di Pareto restano **fissi** per tutto lo sweep — solo `Mtot` viene ricalcolato per ogni combinazione di pesi, il che rende lo sweep economico anche su griglie con migliaia di punti:

```bash
./build/analysis/pareto_analysis/Release/pareto_selector \
    --events   output/optimization/events.root \
    --qmap     output/psf_analysis/q_map.tsv \
    --chi2map  output/psf_analysis/chi2_map.tsv \
    --dofmap   output/dof_analysis/dof_map.tsv \
    --resolution output/resolution_analysis/resolution_map.tsv \
    --ee80-max 10.0 \
    --weight-sweep --weight-step 0.05 \
    --weight-sweep-tsv output/pareto_analysis/weight_sweep_results.tsv \
    --output   output/pareto_analysis/pareto_ternary.png
```

| Opzione | Default | Descrizione |
|---|---|---|
| `--weight-sweep` | off | Attiva la modalità sweep pesi invece del plot singolo-peso |
| `--weight-step` | `0.05` | Passo della griglia baricentrica dei 4 pesi (per `0.05` → 1771 combinazioni) |
| `--weight-sweep-tsv` | `output/pareto_analysis/weight_sweep_results.tsv` | TSV aggregato con un vincitore per combinazione di pesi |

**Colonne TSV**: `w_eta w_Q w_dof w_M x1_winner x2_winner Mtot_winner category_id` (una riga per punto della griglia dei pesi).

**Output**: due PNG puliti (tassellazione triangolare esatta, nessuna intersezione/buco per costruzione):
- **Ternario a sviluppo del tetraedro** (`--output`): le 4 facce del tetraedro dei pesi `(w_eta, w_Q, w_dof, w_M)` disegnate come **net** (sviluppo) — un triangolo centrale (faccia `w_M=0`, coordinate baricentriche `(w_eta, w_Q, w_dof)`) più 3 lembi ribaltati verso l'esterno sui lati condivisi, uno per ciascuna delle altre 3 facce (`w_eta=0`, `w_Q=0`, `w_dof=0`); i 3 apici esterni sono lo stesso vertice `w_M=1` del tetraedro, ripiegato su ciascuna faccia adiacente. Tecnica standard per diagrammi di composizione quaternari (4 componenti), vedi Shimura & Kemp (2015) in [Bibliografia](#bibliografia). Ogni regione colorata è la proiezione della configurazione vincente per quella combinazione di pesi; il colore non è enumerato in legenda (vedi la mappa fisica per la corrispondenza colore↔(x1,x2)).
- **Mappa fisica (x1,x2)** (stesso path con suffisso `_map`): stessi colori per categoria, ma nel dominio nativo `(x1,x2)` di `q_map`/`chi2_map`/`dof_map`, con l'intero dominio simulato riempito (grigio chiaro = valido, nero = invalido/fuori vincoli geometrici) e i punti del fronte mai vincenti come cerchi vuoti.

**Caveat**: essendo `Mtot` una scalarizzazione lineare, la mappa mostra solo configurazioni sull'involucro convesso del fronte di Pareto (vedi Das & Dennis 1997 in [Bibliografia](#bibliografia)).

Invocabile anche per l'intera pipeline via `scripts/pareto_runner.py`, aggiungendo la sezione `"weight_sweep"` a `config/pareto_runs.json`:
```json
"weight_sweep": {"step": 0.05, "ee80_max": 10.0}
```
Output in `output/lens_simulations/{l1_id}_{l2_id}/pareto/weight_sweep/`.

---

### exp1\_main

**Scopo**: analisi statistica frame-by-frame di immagini FITS 16-bit per il confronto di configurazioni ottiche in laboratorio. Confronta una configurazione "buona" e una "cattiva" con un fondo comune; produce mappe di intensità stackate, profili integrati e significatività statistica del segnale.

```bash
./build/analysis/exp1/Release/exp1_main \
    --data-dir analysis/exp1/data \
    --good good --bad bad --bg background \
    --method sigma --sigma 3.0 --iter 3 \
    --output output/exp1/
```

| Opzione | Default | Descrizione |
|---|---|---|
| `--data-dir` | `analysis/exp1/data` | Radice cartella dati |
| `--good` | `good` | Sottocartella configurazione buona |
| `--bad` | `bad` | Sottocartella configurazione cattiva |
| `--bg` | `background` | Sottocartella fondo |
| `--output` | `output/exp1/` | Cartella di output |
| `--frames N` | 0 (tutti) | Limita il numero di frame |
| `--roi x0 y0 x1 y1` | intero frame | Region of interest |
| `--method` | `sigma` | Metodo di stacking: `sigma`, `mean`, `median` |
| `--sigma N` | 3.0 | Soglia sigma-clipping |
| `--iter N` | 3 | Iterazioni sigma-clipping |
| `--no-root` | off | Disabilita salvataggio ROOT |
| `--no-png` | off | Disabilita salvataggio PNG |

**Output**: `stacked_means.png`, `sigma_maps.png`, `diff_maps.png`, `profiles.png`, `summary.png`.

---

### exp2\_main

**Scopo**: confronto quantitativo di quattro configurazioni ottiche (good/bad × focus/nofocus) su stack di immagini FITS 16-bit acquisite con un laser. Estrae il profilo trasversale della traccia laser, misura la larghezza PSF e la linearità del centroide.

**Pipeline:**
```
FITS frames (signal + background)
  → sigma-clipping stack
    → differenza signal − background
      → stima angolo traccia (PCA + momento di inerzia)
        → raffinamento angolo iterativo (≤2 iter ODR)
          → estrazione profilo a slice (centroide pesato ISO 11146)
            → fit ODR lineare sul centroide
              → σ_minor / σ_mean / χ²/ndof
```

**Metriche principali:**

| Metrica | Descrizione |
|---|---|
| `σ_minor` | Semiasse minore della distribuzione 2D (invariante per rotazione; confrontabile tra blob e streak) |
| `σ_mean` | Media pesata inversa-varianza delle larghezze Gaussiane delle slice (valida solo per tracce lineari) |
| `aspect_ratio` | `σ_major / σ_minor`; se < `--min-aspect-ratio` (default 2.0), il blob è circolare e la traccia non viene estratta |
| `χ²/ndof` | Qualità del fit ODR lineare sul centroide (≈1 per traccia lineare ideale) |

**Algoritmi:**

**Centroide pesato (ISO 11146)** — metodo principale:
```
w(s)       = max(0, I(s) − B)         B = mediana dei bordi della slice
centroid   = Σ(s · w(s)) / Σw
σ_centroid = σ_dist / √N_eff
N_eff      = (Σw)² / Σ(w²)
```
Robusto a profili non-gaussiani e a profili troncati ai bordi del FOV.

**Raffinamento iterativo dell'angolo**: l'angolo iniziale (PCA) viene corretto fino a 2 iterazioni tramite il pendio del fit ODR (`δ = atan(a_ODR)`), interrotto se |δ| < 0.05° o |δ| > 10°.

**Trimming basato su `center_err`**: il segmento valido è il blocco contiguo più lungo dove `center_err ≤ trim_max_center_err` (default 5.0 px). Sostituisce il trimming SNR precedente, preservando slice a basso SNR con centroide ben stimabile.

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
    --x1-bad  60  --x2-bad  125 \
    --xdet-good-opt 165 --xdet-bad-opt 150
```

| Parametro | Default | Descrizione |
|---|---|---|
| `--min-aspect-ratio` | 2.0 | Soglia aspect_ratio per rilevare tracce analizzabili |
| `--trim-max-center-err` | 5.0 px | Soglia max `center_err` per il trimming |
| `--center-err-floor` | 0.2 px | Incertezza minima centroide |
| `--center-err-scale` | 1.0 | Fattore moltiplicativo incertezza centroide |
| `--sigma-err-floor` | 0.2 px | Incertezza minima sigma |
| `--sigma-err-scale` | 1.0 | Fattore moltiplicativo incertezza sigma |
| `--trim-pad-slices` | 10 | Padding slice intorno alla regione valida |
| `--xdet-good-opt` | — | Posizione detector ottimale config. buona [mm] |
| `--xdet-bad-opt` | — | Posizione detector ottimale config. cattiva [mm] |

---

### exp3\_main

**Scopo**: confronto sperimentale su banco ottico (display smartphone) tra qualità
sperimentale Q_exp e qualità simulata Q_sim, più stima del fattore di
non-linearità M(r) (Exp3b). Modalità separate, eseguite in sequenza
(`calib` → `measure` → `exp3b` → `plots`); `fits-preview` è indipendente e
serve solo per ispezione visiva dei frame FITS grezzi.

**Pipeline:**
```
Linea di calibrazione FITS → sottrazione background
  → scale_d{dist}mm.json (mm/px sensore, L_sensore)

Tracce laser parallele FITS (Exp3) → sottrazione background
  → centroide pesato per slice → fit ODR lineare (outlier rejection Mahalanobis)
    → Q_exp = χ²/ndof per traccia → q_exp_map.tsv
      → confronto con Q_sim (fisso da config o nearest-neighbor su q_map.tsv)
        → q_comparison.png, q_vs_dax_r{idx}.png

Colonna di dot FITS (Exp3b) → M(r) locale/globale → M_summary.tsv
  → M_vs_dax.png, M_profile_d{dist}_{config}.png, M_nonlinearity_d{dist}_{config}.png

FITS grezzi (qualsiasi cartella) → fits-preview → PNG 1:1 per ispezione
```

```bash
./build/analysis/exp3/Release/exp3_main --mode calib   --config config/exp3/exp3_config.json --data-dir data/exp3/ --output output/exp3/
./build/analysis/exp3/Release/exp3_main --mode measure --config config/exp3/exp3_config.json --data-dir data/exp3/ --output output/exp3/
./build/analysis/exp3/Release/exp3_main --mode exp3b   --config config/exp3/exp3_config.json --data-dir-b data/exp3b/ --output-b output/exp3b/
./build/analysis/exp3/Release/exp3_main --mode plots   --config config/exp3/exp3_config.json --output output/exp3/ --output-b output/exp3b/
```

| Opzione | Default | Descrizione |
|---|---|---|
| `--mode` | *(obbligatoria)* | `calib`, `measure`, `exp3b`, `plots`, `fits-preview` |
| `--config` | `config/exp3/exp3_config.json` | Parametri exp3 (distanze, soglie estrazione traccia, config Q_sim) |
| `--data-dir` | `data/exp3/` | Radice dati Exp3 (`calib`, `measure`, `fits-preview`) |
| `--data-dir-b` | `data/exp3b/` | Radice dati Exp3b (`exp3b`) |
| `--output` | `output/exp3/` | Radice output Exp3 |
| `--output-b` | `output/exp3b/` | Radice output Exp3b |
| `--lens-config` | `all` | `good`, `bad`, o `all` (per `measure`/`exp3b`) |
| `--no-png` | off | Salta la generazione PNG in `plots` |
| `--verbose` | off | Output dettagliato |

**Struttura dati attesa:**
```
data/exp3/
├── calib/d{dist}/{signal,background}/       # FITS linea di calibrazione a distanza d [mm]
└── {good,bad}/d{dist}/{signal,background}/  # FITS tracce laser parallele

data/exp3b/
└── {good,bad}/d{dist}/{signal,background}/  # FITS colonna di dot
```

**Output:**
- `output/exp3/calib/scale_d{dist}mm.json` — scala mm/px e lunghezza sensore per distanza
- `output/exp3/{good,bad}/q_exp_map.tsv` — Q_exp per traccia (d_ax, r_idx, χ²/ndof, n_valid_slices, warning)
- `output/exp3/q_comparison.png` — Q_exp(d_ax) good/bad vs Q_sim di riferimento (linee tratteggiate; valori fissi impostabili in `q_comparison.q_sim_good`/`q_sim_bad` nel config, altrimenti nearest-neighbor da `q_map_tsv`)
- `output/exp3/q_vs_dax_r{idx}.png` — Q_exp vs distanza assiale per ogni offset radiale
- `output/exp3b/{good,bad}/M_summary.tsv` — M globale e residuo RMS per distanza
- `output/exp3b/M_vs_dax.png`, `M_profile_d{dist}_{config}.png`, `M_nonlinearity_d{dist}_{config}.png`
- `output/exp3/fits_preview/<stessa struttura di --data-dir>.png` — anteprima 1:1 di ogni FITS grezzo (titolo con config/d_ax/categoria/esposizione/frame, assi in px, colorbar in ADU)

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
| `TracePoint` | Punto della traccia: `t`, `r`, `x_src`, `y_src`, `z_src`, `mu_y`, `mu_z`, `cov`, `valid`, `n_hits` |
| `LensConfig` | Chiave `(x1, x2)` con confronto a tolleranza 1e-4 mm |
| `LineFitResult` | Risultato ODR: `a`, `b`, `sigma_a`, `sigma_b`, `cov_ab`, `chi2`, `ndof`, `chi2_ndof`, `n_iter`, `converged`, residui, pull |
| `QConfig` | Parametri per `compute_Q`: dimensioni scintillatore, `n_tracks`, `trace_dt`, soglie, unfolding |
| `QResult` | Risultato `compute_Q`: `Q`, `n_traces`, `n_failed`, `n_invalid`, `config_valid`, chi2 per-traccia |
| `CoverageResult` | Risultato `compute_coverage`: `coverage`, `n_y0_evaluated`, `config_valid` |

#### Funzioni

**`load_psf_database(path)`** — carica `psf_data.root` in memoria come `PSFDatabase`. I punti sono ordinati per (x\_source, y\_source) crescente.

**`find_nearest_config(cfg, db)`** — trova la configurazione più vicina (distanza euclidea in (x1,x2)). Warning se distanza > 1e-4 mm.

**`interpolate(x, y, cfg, db)`** — interpolazione bilineare di `mu_y`, `mu_z` e `Σ`. Clamp agli estremi; per y < y\_min usa interpolazione lineare verso l'origine con covarianza isotropizzata.

**`build_trace_3d(p1, p2, cfg, db, dt)`** — traccia media per un segmento 3D P1→P2. Calcola la distanza radiale r dall'asse X, chiama `interpolate(x, r, ...)` e ruota media e covarianza dell'angolo azimutale φ = atan2(z, y).

**`build_trace(y0, cfg, db, L, dt)`** — wrapper: segmento rettilineo lungo X a distanza y₀, con P1 = (−L/2, y₀, 0) e P2 = (+L/2, y₀, 0).

**`is_trace_valid(trace, point_valid_fraction)`** — ritorna `true` se la frazione di `TracePoint` con `valid == true` supera la soglia (default 0.75).

**`fit_trace(trace, min_hits_per_point, max_iter, tol)`** — fit ODR pesato della traccia.

**`compute_Q(cfg, db, qcfg, include_non_converged)`** — calcola Q(x1,x2) su tracce casuali 3D.

**`compute_coverage(cfg, db, qcfg)`** — copertura geometrica media (nessun fit ODR).

#### `fit_trace` — Fit ODR iterativo

Fit della retta `z = a·y + b` sui punti `{(mu_y_i, mu_z_i)}` con `valid == true` e `n_hits >= min_hits_per_point`, usando le matrici di covarianza `Σ_i` come peso via IRLS (Boggs & Rogers, 1990).

**Algoritmo:**
1. Stima iniziale `(a, b)` con OLS non pesato.
2. Per ogni iterazione: calcola il vettore normale `n̂ = (−a, 1)/‖…‖`, poi i pesi `w_i = 1 / (n̂ᵀ Σ_i n̂)` con floor `1e-6 mm²`, risolve il sistema WLS 2×2 in forma chiusa, controlla convergenza su `|Δa| < tol`.
3. Calcola χ², ndof, residui perpendicolari e pull finali.

#### `compute_Q` — Funzione di qualità

Genera `n_tracks` tracce casuali nello scintillatore (punti su facce casuali del parallelepipedo `scint_x × scint_y × scint_z`), per ognuna:
1. Costruisce la traccia 3D con `build_trace_3d`.
2. Verifica la validità (≥75% dei punti on-detector).
3. Applica il temporal unfolding (se abilitato).
4. Esegue `fit_trace`.
5. Accumula `chi2_ndof`.

```
Q(x1, x2) = (1 / n_traces) · Σ chi²_ndof
```

**`QConfig` — campi principali:**

| Campo | Default | Descrizione |
|---|---|---|
| `scint_x`, `scint_y`, `scint_z` | 60, 20, 20 mm | Dimensioni scintillatore |
| `n_tracks` | 100 | Tracce casuali per configurazione |
| `trace_dt` | 0.1 mm | Passo campionamento traccia |
| `min_hits_per_point` | 10.0 | Hit minime per punto PSF valido |
| `trace_valid_fraction` | 0.75 | Soglia tracce valide per configurazione valida |
| `apply_temporal_unfolding` | true | Abilita/disabilita temporal unfolding |
| `z_unfold_step` | 0.0 (auto) | Passo di srotolamento fisso [mm/passo]; 0 = L/(N−1) |

### Strumenti basati su psf\_analysis

#### trace\_viewer

Visualizza la traccia media sul detector con ellissi di covarianza (palette Rainbow: blu = inizio, rosso = fine). Supporta modalità 2D (sorgente a y₀ fissato) e 3D (segmento P1→P2 arbitrario).
```bash
# Modalità 2D
./build/analysis/psf_analysis/Release/trace_viewer \
    --x1 94.9 --x2 186.4 --y0 5.0 --fit

# Modalità 3D con unfolding
./build/analysis/psf_analysis/Release/trace_viewer \
    --x1 108.0 --x2 153.0 \
    --p1x 30.0 --p1y 10.0 --p1z 10.0 \
    --p2x -30.0 --p2y -10.0 --p2z -9.0 \
    --dt 1.0 --fit --unfold
```

| Opzione | Default | Descrizione |
|---|---|---|
| `--x1`, `--x2` | — | **Obbligatori.** Posizioni lenti [mm] |
| `--y0` | — | Distanza radiale [mm] (modalità 2D) |
| `--p1x/y/z`, `--p2x/y/z` | — | Estremi traccia 3D [mm] |
| `--psf` | `output/psf/psf_data.root` | Database PSF |
| `--output` | auto | PNG di output |
| `--dt` | 0.1 | Step traccia [mm] |
| `--L` | 10.0 | Lunghezza traccia [mm] (solo modalità 2D) |
| `--sigma` | 1.0 | Scala ellissi di covarianza |
| `--fit` | off | Esegue fit ODR e mostra χ²/ndof |
| `--unfold` | off | Applica temporal unfolding prima del fit |

#### q\_map

Genera la mappa 2D di Q(x1,x2) o la mappa di copertura geometrica. Marca automaticamente il minimo di Q (o il massimo di copertura) con un marker a stella rosso.
```bash
./build/analysis/psf_analysis/Release/q_map \
    --psf output/psf/psf_data.root \
    --n-tracks 100 --dt 0.1 --log \
    --tsv output/psf_analysis/q_map.tsv

./build/analysis/psf_analysis/Release/q_map --coverage \
    --psf output/psf/psf_data.root
```

| Opzione | Default | Descrizione |
|---|---|---|
| `--psf` | `output/psf/psf_data.root` | Database PSF |
| `--config` | `config/config.json` | Parametri griglia lenti |
| `--output` | auto | PNG di output |
| `--tsv` | (disabilitato) | Esporta in TSV |
| `--scint-x/y/z` | 60, 20, 20 mm | Dimensioni scintillatore |
| `--n-tracks` | 100 | Tracce casuali per configurazione |
| `--dt` | 0.1 | Step traccia [mm] |
| `--min-hits` | 10.0 | Hit minime per punto PSF valido |
| `--trace-frac` | 0.75 | Soglia tracce valide |
| `--unfold-dz` | 0.0 (auto) | Passo di srotolamento fisso [mm/passo] |
| `--no-unfold` | off | Disabilita temporal unfolding |
| `--coverage` | off | Modalità mappa copertura geometrica |
| `--log` | off | Scala logaritmica |

#### chi2\_map

Genera la mappa 2D della linearità della risposta ottica. Per ogni configurazione, esegue un fit di piano `μ_{y,z} = a + b·y₀ + c·x₀` e riporta il χ²/ndof combinato (Y + Z).
```bash
./build/analysis/psf_analysis/Release/chi2_map \
    --psf output/psf/psf_data.root --log
```

#### psf\_gaussianity

Verifica se la PSF (hit non troncati dal piano virtuale, Tecnica E — `PsfDofRuns/y_virtual_hits,z_virtual_hits`) è ben approssimata da una gaussiana bivariata, separando formalmente il **nucleo** dalla **contaminazione** (luce spuria da riflessione totale interna/scattering ai bordi lente, tipicamente <1% degli hit) prima di applicare il test di curtosi multivariata di Mardia. Il test di Mardia sull'intero campione grezzo è estremamente sensibile anche a contaminazioni minime e da solo non distingue "nucleo non gaussiano" da "nucleo gaussiano + coda di contaminazione separata".

```bash
./build/analysis/psf_analysis/Release/psf_gaussianity \
    --input output/lens_simulation/<coppia>/data/lens.root \
    --output output/psf_analysis/psf_gaussianity_analysis.png
```

Genera tre immagini separate, una per grafico (titoli/assi/legenda a piena larghezza): `psf_gaussianity_analysis_scatter.png` (hit sul detector + ellissi 1/2/3σ), `psf_gaussianity_analysis_marginal.png` (profilo marginale + fit gaussiano), `psf_gaussianity_analysis_qq.png` (Q-Q plot vs χ²₂). I valori numerici (cutoff Bonferroni, contaminazione, Mardia γ₂, Henze-Zirkler) sono riportati solo su console/log, non nei grafici.

| Opzione | Default | Descrizione |
|---|---|---|
| `--input` | — | **Obbligatorio.** File ROOT con `Configurations`+`Runs`/`PsfDofRuns` |
| `--output` | auto | Percorso base per i 3 PNG di output (suffissi `_scatter`/`_marginal`/`_qq`) |
| `--sigma-cut` | 3.0 | Cutoff σ del filtro Mahalanobis iterativo (`compute_psf`) |
| `--n-iter` | 4 | Iterazioni del filtro |
| `--min-hits` | — | Hit minime per run valido |
| `--contamination-alpha` | 0.05 | Errore family-wise α per il cutoff Bonferroni-χ²₂ nucleo/contaminazione |
| `--hz-subsample-size` | 2000 | Dimensione sottocampione per il test Henze-Zirkler (seed fisso, riproducibile) |

**Metodo**: sotto H0 (nucleo gaussiano bivariato) la distanza di Mahalanobis quadra D² di ogni hit segue χ²₂, la cui coda ha forma chiusa esatta `P(D²>x) = exp(-x/2)`. Con N hit pooled ed errore family-wise α_fw, il cutoff Bonferroni-corretto è `x_cutoff = -2·ln(α_fw/N)`; gli hit con D² oltre il cutoff sono classificati come contaminazione e la curtosi di Mardia γ₂ = mean(D⁴)/8 − 1 viene ricalcolata solo sul nucleo (D² ≤ cutoff), atteso vicino a 0 per una gaussiana bivariata. A titolo di diagnostica secondaria (più potente di Mardia ma O(n²), quindi limitata a un sottocampione del nucleo) viene riportato anche il test di Henze-Zirkler. Nota: a N dell'ordine di 10⁷ hit, l'errore standard asintotico di γ₂ è talmente piccolo che anche deviazioni di piccola entità pratica risultano formalmente "significative" — l'entità della riduzione di |γ₂| dopo il taglio (non il solo p-value) è l'indicatore rilevante. Dettagli e citazioni in [Riferimenti metodologici § Test di gaussianità multivariata](#test-di-gaussianità-multivariata-psf_gaussianity).

**Origine fisica della contaminazione**: la frazione di hit oltre il cutoff (tipicamente ~0.01–0.1% del campione) non è un artefatto del fit ma un effetto ottico reale del sistema. La fisica ottica di Geant4 (`G4OpticalPhysics` in `src/common/physics_list.cpp`, mai disabilitata) modella la riflessione parziale di Fresnel a ogni interfaccia dielettrica lente/aria (nessuna `G4OpticalSurface` custom è definita, quindi vale il comportamento di default con ~4% di riflettività per superficie in UVFS, n≈1.5). Una piccola frazione di fotoni subisce una riflessione anziché una rifrazione a una delle interfacce (una probabilità di ordine 4%×4%≈0.16% per una doppia riflessione tra le due lenti, coerente in ordine di grandezza con la contaminazione osservata) ed esce con un angolo del tutto diverso da quello atteso, colpendo il detector lontano dal fuoco primario. Questo produce visivamente due firme distinte e coerenti tra loro: (1) nello scatter dei hit sul detector, un nucleo gaussiano stretto circondato da un alone sparso a raggi maggiori; (2) nel Q-Q plot, un andamento sub-diagonale concavo (tipo logaritmo) per il nucleo seguito da un salto quasi verticale oltre il cutoff, dove i pochi punti di contaminazione hanno D² enormi rispetto al resto del campione — la firma tipica di una mistura `(1-ε)·N(μ,Σ) + ε·coda estranea` piuttosto che di code pesanti di un'unica gaussiana.

---

## Temporal Unfolding

Il temporal unfolding risolve la degenerazione del fit ODR per configurazioni con forte aberrazione di campo, in cui la traccia media si ripiega su se stessa nel piano (μ_y, μ_z). La trasformazione applica uno spostamento artificiale progressivo lungo la direzione **ortogonale** alla traccia:
```
ũ_i = u_i + i · δS · n̂⊥
```
dove `n̂⊥` è il vettore normale calcolato dai punti estremi validi, e `δS` è il passo di srotolamento (auto = L/(N−1) o fisso via `z_unfold_step`).

La trasformazione è applicata **esclusivamente** all'interno di `compute_Q` su una copia locale della traccia. Le matrici di covarianza Σ_i non vengono modificate perché l'offset è deterministico: `Cov(ũ_i, ũ_j) = Cov(u_i, u_j)`.

---

## Importance Sampling

L'importance sampling geometrico riduce il numero di fotoni "sprecati" nelle simulazioni Geant4, dirigendo preferenzialmente i fotoni verso l'apertura della prima lente.

**Algoritmo:**
1. Calcola il cono di minima apertura che copre l'intera apertura della lente da qualunque punto del rettangolo sorgente (tangenza sferica sulle superfici R1/R2 + tangenza cilindrica sui bordi).
2. L'asse del cono è il vettore dalla posizione media della sorgente verso il centro della lente.
3. L'angolo di apertura massimo garantisce copertura completa anche per sorgenti al bordo del rettangolo.
4. Geant4 genera fotoni uniformemente nel cono; un peso `w = 1 / (solid_angle_cone / 2π)` compensa la riduzione di angolo solido.

**Attivazione**: `"use_importance_sampling": true` in `config.json`.

**Usato in**: `dof_simulation_main`, `psf_dof_scan_main`, `lens_simulation_main`.

**File**: `src/common/importance_sampling.cpp`, `include/common/importance_sampling.hpp`.

---

## Test unitari

```bash
# Test PSF analysis
./build/analysis/psf_analysis/Release/test_fit_trace
cd build && ctest -R fit_trace_unit -V

# Test Pareto
./build/analysis/pareto_analysis/Release/test_pareto_selector
cd build && ctest -R pareto_selector_unit -V

# Test diagnostica chi2
./build/analysis/psf_analysis/Release/test_chi2_diagnostics

# Test calibrazione exp3
./build/analysis/exp3/Release/test_exp3_calib
```

### Suite T — `fit_trace`

| Test | Scenario | Verifica |
|---|---|---|
| **T1** | Retta `z = y`, cov isotropa | `a=1`, `b=0`, χ²≈0, pull≈0, `n_points_used=21` |
| **T2** | Retta `z = 0.5y + 3`, cov isotropa | `a=0.5`, `b=3`, χ²≈0 |
| **T3** | Retta `z = 2y − 1`, cov anisotropa | Parametri invariati |
| **T4** | Outlier singolo `Δz = 0.5 mm` | Pull ≈ `Δz/σ`, χ²/ndof dominato dall'outlier |
| **T5** | Meno di 3 punti validi | `std::invalid_argument` |
| **T6** | Covarianza degenere (`var = 0`) | Nessun crash, floor `1e-6 mm²`, convergenza |
| **T7** | Cov non diagonale, outlier | `σ_d ≈ √(cov_zz)`, verifica formula `n̂ᵀ Σ n̂` |

### Suite TV — `is_trace_valid`

| Test | Scenario | Verifica |
|---|---|---|
| **TV1** | Tutti `valid=true` | Valida per soglie 0.75 e 1.0 |
| **TV2** | 60% valid | Valida per soglia 50%, invalida per 75% |
| **TV3** | Traccia vuota | `false` |

### Suite TQ — `compute_Q`

| Test | Scenario | Verifica |
|---|---|---|
| **TQ1** | PSF ideale, unfolding OFF | `Q ≈ 0`, `n_failed = 0` |
| **TQ2** | Scintillatore personalizzato | `Q ≈ 0` |
| **TQ3** | Config non presente | `std::invalid_argument` |
| **TQ4** | Due σ_z diverse | `Q(σ_piccola) > Q(σ_grande)` |
| **TQ5** | Traccia ripiegata vs lineare | `χ²(fold) >> χ²(lin)` con unfolding ON |

### Suite Pareto — `test_pareto_selector`

Test unitari per `pareto_core.hpp`: filtri hard (η, DoF, EE80), calcolo Mtot con normalizzazione, fronte di Pareto (dominanza), join su (x1,x2) con tolleranza.

### Suite TD — `test_chi2_diagnostics`

Verifica la diagnostica χ²/ndof di `fit_trace` su dati sintetici con proprietà statistiche note:

| Test | Scenario | Verifica |
|---|---|---|
| **TD_REG** | Regressione: `build_trace_3d` con `n_hits_count` variabile | La covarianza è scalata correttamente per 1/n_hits |
| **TD_ROT** | PSF ruotata (angolo azimutale φ non banale) | `mu_y`/`mu_z`/`cov` ruotati correttamente, `trace(cov)` invariante |
| **TD1** | Baseline, dati perfetti (cov = errore sulla media) | χ²≈0, `ndof=N-2`, fit convergente |
| **TD2** | Scaling della covarianza vs n_hits, rumore i.i.d. | Confronto scenari "over"/"correct"/"under" scaling su 200 realizzazioni Monte Carlo |
| **TD3** | Residui correlati AR(1), covarianza marginale corretta | χ²/ndof confrontato con formula chiusa `expected_chi2_ndof_ar1_linear_fit` per ρ={0, 0.5, 0.9} |
| **TD4** | Combinato: covarianza di distribuzione + rumore AR(1) sulla media | Coerenza tra i due effetti |
| **TD5** | Monotonia χ²/ndof vs densità di campionamento N, lunghezza di correlazione fisica fissa | Andamento atteso al crescere di N |

### Suite C — `test_exp3_calib`

Test unitari per il modulo di calibrazione sperimentale `exp3` (calibrazione via linea, tracce parallele, `exp3b`), interamente su immagini sintetiche costruite in memoria (nessun file FITS reale richiesto):

| Test | Scenario | Verifica |
|---|---|---|
| **C1** | `calibrate_line()` su profilo rettangolare sintetico | Rilevamento soglia, calcolo scala mm/px |
| **C2** | `segment_parallel_lines()` su 3 strisce gaussiane a offset noti | Segmentazione corretta delle linee |
| **C3** | `run_exp3b()` fit lineare su colonna di punti sintetica (`M_true=0.5`, `q_true=5.0`) | Recupero dei parametri M, q |
| **C4** | `load_exp3_config()` su vecchio schema JSON piatto | Retrocompatibilità |
| **C5** | `LineCalib` save/load JSON | Round-trip senza perdita di dati |
| **C6** | `write_q_results_tsv()` | Header contiene `r_idx`, non contiene `theta_deg` |

### Note sui dati sintetici

- **`TracePoint::valid = true`** deve essere impostato esplicitamente (zero-inizializzazione lascia `valid = false`, causando `std::invalid_argument` in `fit_trace`).
- **`PSFPoint::on_detector = true`** deve essere impostato esplicitamente.

---

### Regression testing (pytest)

Suite Python in `tests/` che verifica 10 bug storici regressi. Richiede i binari compilati in `build/Release/` e `build/analysis/`.

```bash
cd tests && pytest -v            # tutta la suite
cd tests && pytest -v -k scan_min  # singolo test per nome
```

| Test | Bug verificato |
|------|----------------|
| `test_nhits_type` | `n_hits` deve essere `Double_t`, non `Int_t` (mismatch writer/reader) |
| `test_mtot_normalization` | `Mtot` ∈ [0, 1] dopo normalizzazione Pareto |
| `test_scan_min_bug` | `scan_min` non deve escludere il punto di fuoco (`max(scan_min, x_virtual)` rimosso) |
| `test_ee80_mean_column` | Colonna `EE80_mean` presente in `resolution_map.tsv` |
| `test_delta_m_column` | Colonna `M_abs_err` presente in `dof_map.tsv` |
| `test_margin_consistency` | `lens_gap_margin` letto da config, non hardcodato |
| `test_all_lenses_catalog` | `--all-lenses` legge entrambi i cataloghi TSV |
| `test_registry_recovery` | Il registry JSON permette di saltare step già completati |
| `test_focus_tsv_passthrough` | `--focus-tsv` viene propagato da `run.sh` a `lens_simulation_main` |
| `test_pareto_join_tolerance` | Join (x1, x2) con tolleranza 1e-3 mm |

CI: `.github/workflows/regression.yml` esegue la suite ad ogni push su `main`.

---

## Bibliografia

Riferimenti usati per progettare la modalità `--weight-sweep` di `pareto_selector` (mappa dei pesi + plot ternario, vedi [Modalità `--weight-sweep`](#pareto_selector) sopra):

- Zintgraf, Kanters, Roijers, Oliehoek, Beau (2015), *"Quality Assessment of Multi-Objective Reinforcement Learning Algorithms"*, weight-space maps — https://ai.vub.ac.be/sites/default/files/efficicent_weights.pdf — tecnica diretta per mappare le regioni dello spazio dei pesi in cui una soluzione è ottima sotto scalarizzazione lineare (usata qui per il plot ternario: ogni regione colorata è la proiezione di questa mappa sul simplex a 4 pesi).
- Dächert, Klamroth, Lacour, Vanderpooten (2023), *"A weight-set decomposition algorithm for the biobjective/multiobjective case"*, Journal of Global Optimization — https://link.springer.com/article/10.1007/s10898-023-01284-x — decomposizione dello spazio dei pesi in regioni associate a soluzioni Pareto-ottime, generalizzazione dell'idea precedente a più obiettivi.
- Das & Dennis (1997), *"A closer look at drawbacks of minimizing weighted sums of objectives for Pareto set generation in multicriteria optimization problems"*, Structural Optimization — limite noto della scalarizzazione lineare (Mtot): può selezionare solo soluzioni sull'involucro convesso del fronte di Pareto, mai regioni concave. **Caveat interpretativo**: il ternario e la mappa (x1,x2) mostrano quindi solo i "vincitori" raggiungibili da una combinazione lineare dei pesi — configurazioni Pareto-ottime ma in una regione concava del fronte non compariranno mai come vincitrici, indipendentemente dal peso usato.
- Shimura & Kemp (2015), *"Tetrahedral plot diagram: A geometrical solution for quaternary systems"*, American Mineralogist 100, DOI 10.2138/am-2015-5371 — sviluppo (net) del tetraedro regolare come tassellazione piana per diagrammi di composizione a 4 componenti; tecnica usata qui per `pareto_ternary.png`, che disegna le 4 facce del tetraedro dei pesi `(w_eta, w_Q, w_dof, w_M)` come un triangolo centrale più 3 lembi ribaltati sui lati condivisi, invece di fette 2D a `w_M` fisso.

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
| `--n-tracks 100` | 100 tracce | Compromesso stabilità statistica / tempo (~22 s per 16 chunk) |
| `--dt 0.1` | 0.1 mm | ~101 punti per traccia da 10 mm; sufficiente per il fit ODR |
| `--unfold-dz 0.000002` | 2×10⁻⁶ mm/passo | Rompe la degenerazione senza distorcere χ² (offset totale ≈ 0.0002 mm ≪ σ_PSF) |
| `--trace-frac 0.50` | 50% | Permette di mappare configurazioni borderline |
| `--min-hits 10` | 10 hit | Soglia conservativa per PSF affidabile |
| `--log` | log | Q copre 2–3 ordini di grandezza |

---

## File di configurazione

### `config/config.json`
```json
{
  "x_det": 579.8,
  "lens_det_gap": 10.0,
  "lens_gap_margin": 1.0,
  "x_min": 45.0,
  "x_max": 300.0,
  "dx": 1.0,

  "r1": 38.6,
  "h1": 12.5,
  "r2": 30.9,
  "h2": 16.3,

  "source_width": 10.0,
  "source_height": 5.0,
  "source_thickness": 60.0,
  "source_x_min": -30.0,
  "source_x_max": 30.0,
  "source_dx": 0.5,
  "source_y_min": 0.0,
  "source_y_max": 14.14,
  "source_dy": 0.5,

  "n_photons": 10000,
  "lens_n_photons": 10000,
  "dof_n_photons": 100000,
  "psf_dof_n_photons": 10000,

  "use_importance_sampling": true,
  "lower_percentile": 0.0,
  "upper_percentile": 0.0,

  "dof_x_virtual_offset": 30.0,
  "dof_x_scan_min": 60.0,
  "dof_x_scan_max": 700.0,
  "dof_x_scan_step": 0.5,
  "dof_k_threshold": 1.414,
  "dof_source_halfy": 5.0,
  "dof_source_halfx": 30.0,
  "dof_source_y_centre": 5.0,

  "psf_dof_x_virtual_offset": 30.0,
  "psf_dof_save_hits": false,

  "m_target": 0.1333
}
```

| Parametro | Descrizione |
|---|---|
| `x_det` | Posizione nominale del detector [mm] |
| `lens_det_gap` | Gap minimo lente→detector [mm]; vincolo per `x_focus` valido |
| `lens_gap_margin` | Margine aggiuntivo tra le lenti [mm] |
| `x_min`, `x_max`, `dx` | Range e passo della scansione lenti [mm] |
| `r1`, `h1` / `r2`, `h2` | Raggio e spessore lenti default (solo GDML) [mm] |
| `source_width/height/thickness` | Dimensioni sorgente per `optimization` [mm] |
| `n_photons` | Fotoni per run in `optimization` |
| `lens_n_photons` | Fotoni per run in `lens_simulation` |
| `dof_n_photons` | Fotoni per run in `dof_simulation` |
| `psf_dof_n_photons` | Fotoni per run in `psf_dof_scan` |
| `use_importance_sampling` | Abilita campionamento geometrico preferenziale |
| `source_x_min/max`, `source_dx` | Griglia X sorgente per `lens_simulation` [mm] |
| `source_y_min/max`, `source_dy` | Griglia Y sorgente per `lens_simulation` [mm] |
| `lower_percentile`, `upper_percentile` | Taglio scala colori per `plot2D` |
| `dof_x_virtual_offset` | Offset piano virtuale da x2 [mm] |
| `dof_x_scan_min/max/step` | Range scan x_det per `dof_map` [mm] |
| `dof_k_threshold` | k per DoF = k·σ_z_min (default √2) |
| `dof_source_halfy/halfx` | Semi-estensione sorgente DoF [mm] |
| `dof_source_y_centre` | Centro sorgente Y per DoF [mm] |
| `psf_dof_x_virtual_offset` | Offset piano virtuale per `psf_dof_scan` [mm] |
| `psf_dof_save_hits` | Salva hit individuali in `psf_dof.root` (pesante) |
| `m_target` | Magnificazione target per validazione in `dof_map` |

I parametri `config_id_offset` e `run_id_offset` vengono aggiunti dinamicamente da `run.sh` per garantire ID univoci dopo il merge. In modalità super-batch (invocata da `lens_runner.py --disk-guard-gb`), i campi `config_id_offset_base` e `run_id_offset_base` nel JSON del chunk specificano l'offset globale accumulato tra batch sequenziali.

### `config/config_profile.json`

Configurazione ridotta per profiling rapido (`dx = 30 mm`):
```json
{
  "x_min": 33.0, "x_max": 171.0, "dx": 30.0,
  "r1": 38.6, "h1": 12.5, "r2": 30.9, "h2": 16.3,
  "lower_percentile": 0.45, "upper_percentile": 0.0
}
```

### `config/analysis_params.json`

Parametri di tuning usati esclusivamente da `autonomous_optimizer.py` per configurare ogni step della pipeline di analisi. Tracciato da git con i valori di riferimento; modificare localmente senza commit per campagne specifiche. Valori di riferimento:

```json
{
  "psf_extractor": {
    "min_hits": 50
  },
  "q_map": {
    "n_tracks": 1000,
    "dt": 0.5,
    "min_hits": 50,
    "trace_frac": 0.75,
    "dist_to_target": true
  },
  "chi2_map": {
    "min_hits": 50,
    "p_low": 0.0,
    "p_high": 100.0,
    "adaptive_target": true
  },
  "dof_map": {
    "core_fraction": 1.0,
    "m_target": 0.2
  },
  "resolution_map": {},
  "pareto_selector": {
    "ee80_max": 10.0,
    "w_eta": 0.1,
    "w_Q": 0.1,
    "w_dof": 0.3,
    "w_M": 0.5
  }
}
```

| Chiave | Descrizione |
|--------|-------------|
| `psf_extractor.min_hits` | Hit minime per marcare un run `on_detector = true` |
| `q_map.n_tracks` | Tracce casuali per configurazione nel calcolo di Q |
| `q_map.dist_to_target` | Se `true`, Q è la distanza al target invece di χ²/ndof grezzo |
| `chi2_map.adaptive_target` | Target χ² adattivo per range configurabile con `p_low`/`p_high` |
| `dof_map.m_target` | Magnificazione target per calcolo `M_abs_err` |
| `pareto_selector.w_*` | Pesi per Mtot: `w_eta` (η), `w_Q` (qualità), `w_dof` (DoF), `w_M` (magnificazione) |

### Varianti di config per velocità/validazione

Oltre a `config.json` (default), `config/` contiene varianti con passo di scan e numero di fotoni diversi, per bilanciare velocità e precisione o per le validazioni delle Tecniche C/E:

| File | dx | lens/psf_dof photons | Uso |
|---|---|---|---|
| `config_detailed.json` | 1.0 | 10000 | Scan fine, pipeline completa |
| `config_intermediate.json` | 2.0 | 5000 | Compromesso velocità/precisione |
| `config_fast.json` | 3.0 | 1000 | Iterazione rapida |
| `config_screening.json` | 3.0 | 1000 (n_photons 1000) | Screening `--all-lenses` |
| `config_test.json` | 10.0 | 100 | Smoke test pipeline |
| `config_psf_gaussianity.json` | 1.0, range 100-150 | lens_n_photons 3M | Statistica per test di gaussianità |
| `config_validate_c_10k.json` / `_5k.json` | 50.0, range 100-200 | 10000 / 5000 | Validazione Tecnica C |
| `config_validate_e.json` | 50.0, range 100-200 | 2000 | Validazione Tecnica E |

---

## Formato dei file ROOT

I file ROOT usano **vettori di hit** per massimizzare l'efficienza di archiviazione (riduzione ~40% rispetto al formato riga-per-hit). Compressione LZ4 livello 4 (`SetCompressionLevel(404)`).

### `lens.root` (lens\_simulation)
```
TTree "Configurations"
  config_id  Int_t
  x1, x2     Double_t [mm]
  l1_id  String
  l2_id  String

TTree "Runs"
  run_id     Int_t
  config_id  Int_t
  x_source   Float_t  [mm]
  y_source   Float_t  [mm]
  n_hits     Double_t
  y_hits     vector<float> [mm]
  z_hits     vector<float> [mm]
```

### `events.root` (optimization)
```
TTree "Configurations"
  config_id  Int_t
  x1, x2     Double_t
  l1_id  String
  l2_id  String

TTree "Efficiency"
  config_id  Int_t
  n_photons  Int_t
  n_hits     Double_t
```

### `psf_data.root` (psf\_extractor)
```
TTree "PSF"
  config_id               Int_t
  x1, x2                  Double_t [mm]
  x_source                Float_t  [mm]
  y_source                Float_t  [mm]
  mean_y, mean_z          Double_t [mm]
  cov_yy, cov_yz, cov_zz  Double_t [mm²]
  n_hits_filtered         Double_t
  n_hits_raw              Double_t
  on_detector             Bool_t
```

### `focal.root` (dof\_simulation)
```
TTree "FocalConfigurations"
  config_id   Int_t
  x1, x2      Double_t [mm]
  x_virtual   Double_t [mm]
  l1_id   String
  l2_id   String

TTree "FocalRays"
  config_id      Int_t
  n_rays         Int_t
  y_hits         vector<float>  posizione Y al piano virtuale [mm]
  z_hits         vector<float>  posizione Z al piano virtuale [mm]
  dy_hits        vector<float>  direzione py/px
  dz_hits        vector<float>  direzione pz/px
  weight_hits    vector<float>  peso importance sampling
  y_source_hits  vector<float>  posizione Y sorgente per ogni raggio [mm]
```

### `psf_dof.root` (psf\_dof\_scan)

Estende il formato PSF con i dati di raggio al piano virtuale per ogni punto della griglia (config\_id, x\_source, y\_source). Contiene sia i TTree PSF-style che i vettori di raggi DoF.

---

## Workflow completo

Il workflow primario per una coppia di lenti è gestito da `lens_runner.py`, che esegue automaticamente tutti i passi della pipeline (opt → lens → PSF → DoF → mappe → Pareto).

```bash
# ── 1. Build ────────────────────────────────────────────────────────────────
cmake -S . -B build/ -G "Ninja Multi-Config"
cmake --build build/ --config Release

# ── 2. Test unitari ─────────────────────────────────────────────────────────
./build/analysis/psf_analysis/Release/test_fit_trace
./build/analysis/pareto_analysis/Release/test_pareto_selector
./build/analysis/exp3/Release/test_exp3_calib
# oppure: cd build && ctest -V

# ── 3. Pipeline per singola coppia ──────────────────────────────────────────
python3 scripts/lens_runner.py --l1-id LA4464 --l2-id LA4464R --local --mobile-focus
# → output/lens_simulations/LA4464_LA4464R/{plots,data,pipeline.log}

# ── 4. Post-processing Pareto con pesi multipli ──────────────────────────────
python3 scripts/pareto_runner.py --l1-id LA4464 --l2-id LA4464R
# → output/lens_simulations/LA4464_LA4464R/pareto/{label}/

# ── 5. Diagnostica manuale configurazione ottimale ───────────────────────────
./build/analysis/psf_analysis/Release/trace_viewer \
    --x1 <x1*> --x2 <x2*> --y0 5.0 --fit
```

---

## Riferimenti metodologici

### Metodo di propagazione raggi (`dof_map`)

**Propagazione lineare tra piani**. Nel tratto `x_virtual → detector` non esistono elementi ottici, quindi i raggi si propagano in linea retta: `z(xᵢ) = z₀ + dz·(xᵢ − x_virtual)`. Le direzioni `(dy, dz)` codificano già la rifrazione di entrambe le lenti calcolata da Geant4, aberrazioni di ordine superiore incluse.

**Ricerca del fuoco**: il piano di fuoco `x*` minimizza σ_z(x). Affinamento sub-step con interpolazione parabolica del vertice di σ_z²(x) (fit di Lagrange sui tre punti più vicini). Equivalente al criterio "Quick Focus" di Zemax per minimizzazione RMS dello spot.

**Fuoco prima della seconda lente**: quando il minimo di σ_z cade prima di x2, la propagazione lineare all'indietro attraverserebbe la seconda lente, dove i raggi verrebbero rifratti: il risultato è fisicamente non valido. Questi bin vengono delimitati nelle mappe PNG da una linea di contorno tratteggiata arancione (`CONT3`, kOrange+7, kDashed) sovrapposta a ciascuna heatmap, e sono sintetizzati in una mappa dedicata `dof_focus_zone_map.png` con palette a due fermate (blu = fuoco valido, rosso = fuoco prima della lente 2). I bin sono anche marcati con `focus_before_lens2 = 1` nel TSV.

- **Ray transfer matrix analysis (ABCD)** — Usata come termine di confronto: appropriata per sistemi con parametri gaussiani e approssimazione parAssiale; non usata qui perché il campionamento MC include aberrazioni di ordine superiore. URL: https://en.wikipedia.org/wiki/Ray_transfer_matrix_analysis
- **DevOptical Part 19: A Quick Focus Algorithm** — The Pulsar. Criterio di best focus per minimizzazione RMS spot. URL: https://www.thepulsar.be/article/-devoptical-part-19--a-quick-focus-algorithm
- **Zemax OpticStudio Manual** — Ansys. Conferma che il piano di best focus minimizza l'RMS radiale dei raggi su scansione assiale. URL: https://neurophysics.ucsd.edu/Manuals/Zemax/ZemaxManual.pdf

---

### Magnificazione ottica (`M`)

`M` è calcolato come pendenza di regressione lineare pesata di `y_focus` (posizione al piano di fuoco) rispetto a `y_source_hits` (posizione sorgente), per fotone — non come rapporto di deviazioni standard: quella formulazione confondeva la magnificazione geometrica vera con lo spot-blur da aberrazione (dimostrato: correlazione di Pearson 0.976 tra il vecchio `M` ed `EE80` su 27730 punti griglia LA4464/LA4464R), sovrastimandola sistematicamente.

- **UCSD Physics 1C, Lecture 9 — Multiple-lens systems**. Derivazione standard `M_totale = M1 × M2 = (−v1/u1)×(−v2/u2)` per sistemi a lenti sottili in cascata (l'immagine della prima lente è l'oggetto della seconda). URL: https://courses.physics.ucsd.edu/2011/Summer/session1/physics1c/lecture9.pdf
- **Weighted Least Squares** — Wikipedia. Stimatore ai minimi quadrati pesati della pendenza (usato qui per isolare la relazione lineare sistematica `y_focus/y_source_hits` trattando il blur da aberrazione come rumore indipendente attorno alla retta, anziché sommarlo in quadratura nella σ). URL: https://en.wikipedia.org/wiki/Weighted_least_squares

---

### Fit robusto (ODR / IRLS)

- **Boggs, P.T. & Rogers, J.E. (1990)**. *Orthogonal Distance Regression*. Contemporary Mathematics, 112, 183–194. Formulazione teorica dell'IRLS per ODR pesato con matrici di covarianza generali. Implementato in `fit_trace` (`analysis/psf_analysis/psf_interpolator.cpp`).

- **Fischler, M.A. & Bolles, R.C. (1981)**. *Random Sample Consensus: A Paradigm for Model Fitting with Applications to Image Analysis and Automated Cartography*. Communications of the ACM, 24(6), 381–395. Riferimento classico per il fit di modelli in presenza di outlier (base concettuale per il loop ODR con rejection in `fit_centroid_line`).

---

### Estrazione traccia e centroide

- **ISO 11146-1:2005** — *Lasers and laser-related equipment: Test methods for laser beam widths, divergence angles and beam propagation ratios*. Definisce il metodo dei momenti del secondo ordine per la misura della larghezza del fascio (usato per σ_minor, σ_dist). Implementato come centroide pesato in `exp2/` e `exp3/`.

- **Thomas, S. et al. (2006)**. *Comparison of centroid computation algorithms in a Shack–Hartmann sensor*. Monthly Notices of the Royal Astronomical Society, 371(1), 323–336. DOI: [10.1111/j.1365-2966.2006.10661.x](https://doi.org/10.1111/j.1365-2966.2006.10661.x). Dimostra che il centroide pesato per intensità (WCoG) ha bias sistematico ridotto rispetto al fit Gaussiano quando il profilo è non ideale.

- **Vyas, A. et al. (2009)**. *Centroid Detection by Gaussian Pattern Matching in Adaptive Optics*. arXiv:0910.3386. Confronto tra fit Gaussiano e centroide pesato a basso SNR.

- **Zhang, C. & Couloigner, I. (2007)**. *Accurate Centerline Detection and Line Width Estimation of Thick Lines Using the Radon Transform*. IEEE Transactions on Image Processing, 16(2), 310–316. DOI: [10.1109/TIP.2006.887731](https://doi.org/10.1109/TIP.2006.887731). Riferimento per robustezza del centroide pesato nella localizzazione di linee diffuse.

- **PCA per stima orientazione**: la stima dell'angolo della traccia via decomposizione spettrale della matrice di covarianza 2D del profilo di intensità è usata in `exp2/` come inizializzazione prima del raffinamento ODR iterativo.

---

### Metriche di risoluzione ottica (EE80)

- **EE80 (Encircled Energy 80%)**: metrica standard per la caratterizzazione di PSF in strumenti ottici spaziali (HST, JWST, Euclid). Per una distribuzione gaussiana 2D isotropa: `EE80 = 2 · 1.7941 · σ_rms` dove `σ_rms = √((σ_y² + σ_z²)/2)`. Il fattore 1.7941 è la radice del 80° percentile della distribuzione χ²(2) (raggio da cui il 80% dell'energia è contenuta). Implementata in `analysis/resolution_analysis/resolution_map.cpp` e `analysis/dof_analysis/dof_map.cpp`.

---

### Selezione multi-criterio (Pareto)

- **Pareto optimality** — Concetto classico dell'ottimizzazione multi-obiettivo: una soluzione è non-dominata se non esiste nessuna altra soluzione migliore su tutti gli obiettivi simultaneamente. Il fronte di Pareto su `(η, Q)` raccoglie le configurazioni dove migliorare η richiede peggiorare Q e viceversa. Implementato in `analysis/pareto_analysis/pareto_core.hpp`.

---

### Ranking robusto coppie di lenti (mean-variance / LCB)

- **Markowitz, H. (1952)**. *Portfolio Selection*. The Journal of Finance, 7(1), 77–91. DOI: [10.2307/2975974](https://doi.org/10.2307/2975974). Formulazione classica del trade-off media/varianza; base concettuale dello score `score = mean_eta − k·std_eta` usato per il ranking robusto delle coppie di lenti (penalizza coppie con efficienza media alta ma varianza eccessiva sulla griglia x1,x2). Implementato in `analysis/lens_ranking.hpp` (`compute_pair_ranking`).

- **Auer, P., Cesa-Bianchi, N., & Fischer, P. (2002)**. *Finite-time Analysis of the Multiarmed Bandit Problem*. Machine Learning, 47(2-3), 235–256. DOI: [10.1023/A:1013689704352](https://doi.org/10.1023/A:1013689704352). Origine del criterio Lower/Upper Confidence Bound per selezione "safe" sotto incertezza statistica; qui k=0.5 (calibrato sulla distribuzione empirica del coefficiente di variazione, mediana 0.71 su 225 coppie, per non sacrificare eccessivamente l'efficienza assoluta a favore della sola consistenza).

---

### Rilevazione di streak in immagini astronomiche

- **Nir, G. et al.** *pyradon: Python tools for streak detection in astronomical images using the Fast Radon Transform*. GitHub: [guynir42/pyradon](https://github.com/guynir42/pyradon). Riferimento per la trasformata di Radon veloce applicata a immagini con streak diffuse.

- **Yanagisawa, T. et al. (2015)**. *Streak Detection and Analysis Pipeline for Space-debris Optical Images*. ResearchGate. Base per la pipeline di estrazione di features (centroide, larghezza, flusso) da immagini ottiche con streak lineari a basso SNR.

---

### Test di gaussianità multivariata (`psf_gaussianity`)

- **Mardia, K.V. (1970)**. *Measures of multivariate skewness and kurtosis with applications*. Biometrika, 57(3), 519–530. DOI: [10.1093/biomet/57.3.519](https://doi.org/10.1093/biomet/57.3.519). Origine del test di curtosi multivariata γ₂ = mean(D⁴)/8 − 1 già implementato in `psf_gaussianity`; noto per l'estrema sensibilità a qualunque contaminazione, motivo per cui va applicato al nucleo e non al campione grezzo.

- **Rousseeuw, P.J. & Van Driessen, K. (1999)**. *A Fast Algorithm for the Minimum Covariance Determinant Estimator*. Technometrics, 41(3), 212–223. DOI: [10.1080/00401706.1999.10485670](https://doi.org/10.1080/00401706.1999.10485670). Origine concettuale dello stimatore robusto iterativo (MCD) a cui si ispira il filtro Mahalanobis 3σ di `compute_psf` (`analysis/psf_common.hpp`).

- **Hardin, J. & Rocke, D.M. (2005)**. *The Distribution of Robust Distances*. Journal of Computational and Graphical Statistics, 14(4), 928–946. DOI: [10.1198/106186005X77685](https://doi.org/10.1198/106186005X77685). Correzione finite-sample F-scaled per distanze robuste MCD, pensata per quando il sottoinsieme robusto h è una frazione consistente (es. 50%) di n; non applicata qui perché irrilevante quando il filtro scarta <1% degli hit (h/n>99%), regime in cui la coda esatta di χ²₂ è già un'approssimazione adeguata.

- **Cerioli, A. (2010)**. *Multivariate Outlier Detection with High-Breakdown Estimators*. Journal of the American Statistical Association, 105(489), 147–156. DOI: [10.1198/jasa.2009.tm09147](https://doi.org/10.1198/jasa.2009.tm09147). Framework iterato (reweighted-MCD) con controllo formale dell'errore family-wise, di cui il cutoff Bonferroni-χ²₂ usato in `bonferroni_chi2_2_cutoff` è la versione semplificata valida nel regime h/n≈1 di questa analisi.

- **Henze, N. & Zirkler, B. (1990)**. *A class of invariant consistent tests for multivariate normality*. Communications in Statistics - Theory and Methods, 19(10), 3595–3617. DOI: [10.1080/03610929008830400](https://doi.org/10.1080/03610929008830400). Test diagnostico secondario (più potente di Mardia contro varie alternative, ma O(n²)) implementato in `psf_gaussianity` su un sottocampione casuale del nucleo (`--hz-subsample-size`, seed fisso per riproducibilità).

- **Moffat, A.F.J. (1969)**. *A Theoretical Investigation of Focal Stellar Images in the Photographic and Television Fields, and a Comparison with Existing Models*. Astronomy & Astrophysics, 3, 455. Contesto ottico/astronomico: code non gaussiane nelle PSF reali (profili Moffat/King) sono un fenomeno noto e atteso, non necessariamente un artefatto numerico — motiva la separazione nucleo/contaminazione piuttosto che un taglio spaziale arbitrario.

- **Hecht, E. (2017)**. *Optics*, 5th ed., Pearson, cap. 4 (Fresnel equations). Riferimento standard per la riflettività di Fresnel alle interfacce dielettriche (~4% per superficie UVFS/aria, n≈1.5); è il meccanismo fisico — non un artefatto numerico — dietro la coda di contaminazione (hit oltre il cutoff Bonferroni) osservata in `psf_gaussianity`: fotoni riflessi anziché rifratti a una o più interfacce lente-aria (`G4OpticalPhysics`, mai disabilitata in `src/common/physics_list.cpp`) escono con angolo anomalo e colpiscono il detector lontano dal fuoco primario.

---

## References — Pipeline di ottimizzazione RIPTIDE (Tecniche A/C/E)

- **Blyth, S. et al. (2025)**. *Opticks: GPU Accelerated Optical Photon Simulation with Geant4*. arXiv:2502.13215. Contesto per l'accelerazione GPU delle simulazioni ottiche con Geant4; motiva l'approccio alternativo tramite ottimizzazione della voxelizzazione (Tecnica A).

- **Ai, X. et al. (2025)**. *Simulation of the EIC RICH detector with Opticks*. arXiv:2512.06061. Caso d'uso di riferimento per la simulazione di fotoni ottici in rivelatori Cherenkov con Opticks/Geant4.

- **Geant4 Collaboration**. *Geant4 User's Guide for Application Developers — Parallelism and Geometry Optimization*. Documentazione ufficiale: https://geant4-userdoc.web.cern.ch. Base per `G4GeometryManager::RequestParallelOptimisation(true)` (Tecnica A) e la voxelizzazione MT.

- **Chan, T.F., Golub, G.H., LeVeque, R.J. (1979)**. *Updating Formulae and a Pairwise Algorithm for Computing Sample Variances*. Technical Report STAN-CS-79-773, Stanford University. Formula di Welford/Chan per l'aggiornamento online dei momenti del secondo ordine, implementata in `compute_virtual_stats()` in `src/lens_simulation/lens_scan.cpp` (Tecnica E).

- **Hartig, G.F. et al. (2026)**. *Surrogate-based Bayesian Monte Carlo for optical system optimization*. arXiv:2603.13646. Metodologia di riferimento per la riduzione del costo MC tramite surrogate model; motiva la Tecnica C (riduzione N_photons con controllo dell'errore relativo).

- **Chen, M. et al. (2026)**. *GPU-accelerated optical Monte Carlo for large-scale photon transport*. arXiv:2606.05385. Benchmark per lo speedup atteso nella simulazione MC di fotoni ottici su GPU; confronto con l'approccio CPU ottimizzato di RIPTIDE.
