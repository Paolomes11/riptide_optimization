# riptide_optimization

Simulazione Monte Carlo Geant4 per l'ottimizzazione del posizionamento di un sistema ottico a doppia lente (UVFS) davanti a un fotocatodo GaAsP. Il progetto fa parte dell'esperimento **RIPTIDE** e comprende due programmi principali — `optimization` e `lens_simulation` — più tre strumenti di analisi ROOT.

---

## Indice

1. [Requisiti](#requisiti)
2. [Struttura del progetto](#struttura-del-progetto)
3. [Setup fisico](#setup-fisico)
4. [Build](#build)
5. [Programma: optimization](#programma-optimization)
6. [Programma: lens\_simulation](#programma-lens_simulation)
7. [Strumenti di analisi](#strumenti-di-analisi)
8. [File di configurazione](#file-di-configurazione)
9. [Formato dei file ROOT](#formato-dei-file-root)
10. [Workflow completo](#workflow-completo)

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
├── programs/                        # Entry point dei due eseguibili
│   ├── optimization_main.cpp        # Main per la scansione di efficienza
│   └── lens_simulation_main.cpp     # Main per la simulazione beam scan
│
├── src/
│   ├── common/                      # Codice condiviso tra i due programmi
│   │   ├── physics_list.cpp         # Lista di fisica (ottica + standard EM)
│   │   └── primary_generator_action.cpp  # Generatore particelle (GPS)
│   ├── optimization/                # Sorgenti del programma optimization
│   │   ├── detector_construction.cpp
│   │   ├── action_initialization.cpp
│   │   ├── event_action.cpp         # Raccolta hit, scrittura ROOT
│   │   ├── run_action.cpp
│   │   ├── sensitive_detector.cpp   # Fotocatodo sensibile
│   │   └── optimizer.cpp            # Loop di scansione, output events.root
│   └── lens_simulation/             # Sorgenti del programma lens_simulation
│       ├── detector_construction.cpp
│       ├── action_initialization.cpp
│       ├── event_action.cpp         # Raccolta hit, scrittura ROOT
│       ├── run_action.cpp
│       ├── sensitive_detector.cpp
│       └── lens_scan.cpp            # Loop beam scan, output lens.root
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
│   └── CMakeLists.txt
│
├── geometry/                        # Geometria GDML
│   ├── main.gdml                    # File principale (include gli altri)
│   ├── define.xml                   # Costanti, posizioni, proprietà ottiche
│   ├── materials.xml                # Materiali (aria, UVFS, GaAsP, borosilicato)
│   ├── solids.xml                   # Solidi (ellissoidi lenti, fotocatodo)
│   └── structure.xml                # Volumi fisici e loro posizionamento
│
├── macros/
│   ├── optimization.mac             # Sorgente per optimization (10000 fotoni, r=10mm)
│   ├── lens_simulation.mac          # Sorgente per lens_simulation (1000 fotoni, r=0.1mm)
│   ├── run.mac                      # Macro alternativa
│   └── vis.mac                      # Macro per visualizzazione interattiva
│
├── config/
│   └── config.json                  # Parametri di scansione e geometria lenti
│
├── external/
│   └── lyra/                        # Libreria header-only per parsing CLI
│
└── output/                          # Creata automaticamente da CMake
    ├── events.root                  # Output di optimization
    ├── lens_simulation/
    │   └── lens.root                # Output di lens_simulation
    └── mean_covariance_maps/        # Output di m_c_creator
```

---

## Setup fisico

Il sistema simulato è composto da tre elementi posizionati lungo l'asse X, all'interno di un volume mondo cubico di 1000×1000×1000 mm³ riempito d'aria:

```
Sorgente fotoni  →  [Lente 75mm]  →  [Lente 60mm]  →  [Fotocatodo GaAsP 16×16mm]
      (GPS)           (UVFS)           (UVFS)              (sensore)
         x=0        x ≈ 14–170mm    x ≈ 45–186mm           x ≈ 180mm
```

**Lente 75 mm** (`lens75`): ellissoide in UV Fused Silica (UVFS), raggio 38.6 mm, spessore 12.5 mm. Modella una lente piano-convessa tagliata da una sfera.

**Lente 60 mm** (`lens60`): ellissoide UVFS, raggio 30.9 mm, spessore 16.3 mm. Orientata con la faccia convessa verso la sorgente.

**Fotocatodo GaAsP**: lastra quadrata 16×16×0.01 mm, indice di rifrazione 3.5–3.8 nel range 2–4 eV, lunghezza di assorbimento ~1 µm. Registra i fotoni ottici che attraversano la sua superficie frontale (`fGeomBoundary`).

**Sorgente**: fotoni ottici a 2.5 eV (≈ 496 nm), generati da GPS Geant4 su un disco (Annulus) con distribuzione angolare isotropa diretta lungo +X.

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
```

> **Nota**: CMake crea automaticamente le cartelle `output/`, `output/lens_simulation/` e `output/mean_covariance_maps/` durante la configurazione.

---

## Programma: optimization

**Eseguibile**: `optimization_main`

**Scopo**: scansione grezza dell'efficienza geometrica del sistema a doppia lente. Per ogni coppia di posizioni (x1, x2) delle lenti, esegue una simulazione con 10000 fotoni emessi da una sorgente circolare a r=10 mm e conta quanti raggiungono il fotocatodo. Il risultato è una mappa 2D di efficienza.

### Utilizzo

```bash
# Modalità ottimizzazione (scansione completa)
./build/Release/optimization_main -g geometry/main.gdml -o

# Modalità visualizzazione interattiva
./build/Release/optimization_main -g geometry/main.gdml -v

# Modalità batch con macro personalizzata
./build/Release/optimization_main -g geometry/main.gdml -b -m macros/mia_macro.mac

# Help
./build/Release/optimization_main --help
```

### Opzioni CLI

| Flag | Tipo | Descrizione |
|---|---|---|
| `-g`, `--geometry` | path | **Obbligatorio.** Percorso al file GDML della geometria |
| `-m`, `--macro` | path | Macro Geant4 da eseguire (default: `macros/optimization.mac`) |
| `-v`, `--visualize` | flag | Abilita la visualizzazione interattiva OpenGL/Qt |
| `-b`, `--batch` | flag | Modalità batch senza UI |
| `-o`, `--optimize` | flag | Avvia la scansione di ottimizzazione |

### Output

File ROOT: `output/events.root`

Contiene un singolo TTree `events` con una riga per ogni fotone rilevato:

| Branch | Tipo | Descrizione |
|---|---|---|
| `x1` | `Double_t` | Posizione della lente 75mm [mm] |
| `x2` | `Double_t` | Posizione della lente 60mm [mm] |
| `config_id` | `Int_t` | Indice della configurazione (progressivo) |

L'efficienza per configurazione si calcola come `conteggio_hit / N_generati` (N=10000 di default).

### Comportamento interno

Il loop di scansione in `optimizer.cpp` itera su tutte le coppie (x1, x2) fisicamente valide secondo i vincoli geometrici delle lenti (vedi [config.json](#file-di-configurazione)). Per ogni configurazione:
1. Sposta le lenti con `DetectorConstruction::SetLensPositions(x1, x2)`
2. Imposta un seed casuale basato sul clock ad alta risoluzione
3. Esegue la macro Geant4 (`/run/beamOn 10000`)
4. `EventAction::EndOfEventAction` raccoglie le hit dal `SensitivePhotocathode` e le scrive nel TTree

---

## Programma: lens\_simulation

**Eseguibile**: `lens_simulation_main`

**Scopo**: simulazione dettagliata del beam scan. Per ogni configurazione di lenti e per ogni posizione della sorgente lungo Y (da 0 a 10 mm a passi di 0.1 mm), esegue 1000 fotoni con una sorgente puntiforme (r=0.1 mm) e registra la posizione (y, z) di ogni fotone sul fotocatodo. Permette di ricostruire la funzione di risposta del sistema ottico (PSF).

### Utilizzo

```bash
# Esecuzione della simulazione completa
./build/Release/lens_simulation_main -g geometry/main.gdml -l

# Visualizzazione interattiva
./build/Release/lens_simulation_main -g geometry/main.gdml -v

# Con macro personalizzata
./build/Release/lens_simulation_main -g geometry/main.gdml -l -m macros/mia_macro.mac

# Help
./build/Release/lens_simulation_main --help
```

### Opzioni CLI

| Flag | Tipo | Descrizione |
|---|---|---|
| `-g`, `--geometry` | path | **Obbligatorio.** Percorso al file GDML |
| `-m`, `--macro` | path | Macro Geant4 (default: `macros/lens_simulation.mac`) |
| `-v`, `--visualize` | flag | Abilita visualizzazione interattiva |
| `-b`, `--batch` | flag | Modalità batch |
| `-l`, `--lens-sim` | flag | Avvia il beam scan completo |

### Output

File ROOT: `output/lens_simulation/lens.root`

Contiene tre TTree:

**`Configurations`** — una riga per coppia di posizioni lenti:

| Branch | Tipo | Descrizione |
|---|---|---|
| `config_id` | `Int_t` | Indice configurazione (progressivo) |
| `x1` | `Double_t` | Posizione lente 75mm [mm] |
| `x2` | `Double_t` | Posizione lente 60mm [mm] |

**`Runs`** — una riga per ogni esecuzione (configurazione × posizione sorgente):

| Branch | Tipo | Descrizione |
|---|---|---|
| `run_id` | `Int_t` | Indice run (progressivo globale) |
| `config_id` | `Int_t` | Configurazione di appartenenza |
| `x_source` | `Float_t` | Posizione Y della sorgente [mm] |
| `n_hits` | `Int_t` | Numero di fotoni rilevati in questo run |

**`Hits`** — una riga per ogni fotone rilevato:

| Branch | Tipo | Descrizione |
|---|---|---|
| `y_hit` | `Float_t` | Coordinata Y di arrivo sul fotocatodo [mm] |
| `z_hit` | `Float_t` | Coordinata Z di arrivo sul fotocatodo [mm] |

> **Mapping hits → run**: la colonna `n_hits` in `Runs` permette di ricostruire l'associazione senza un `run_id` in `Hits`. Le hit sono scritte in ordine sequenziale; la i-esima riga di `Runs` corrisponde alle `n_hits[i]` righe di `Hits` a partire dall'offset cumulativo `sum(n_hits[0..i-1])`.

### Comportamento interno

Il loop in `lens_scan.cpp` ha tre livelli annidati:
1. **x1** — posizioni fisicamente valide della lente 75mm
2. **x2** — posizioni fisicamente valide della lente 60mm (con gap minimo tra le lenti)
3. **y_source** — da 0 a 10 mm con passo 0.1 mm (101 posizioni)

Per ogni run:
- Si aggiorna la posizione GPS con `/gps/pos/centre 0 y_source 0 mm`
- Si esegue la macro (`/run/beamOn 1000`)
- `EventAction` accumula le hit evento per evento; a fine run, `lens_scan` legge il conteggio tramite `GetLastRunHitCount()` e lo salva in `Runs.n_hits`

---

## Strumenti di analisi

Tutti e tre gli eseguibili di analisi si compilano insieme al progetto principale e si trovano in `build/analysis/`.

### plot2D

**Legge**: `output/events.root`

Genera una mappa 2D di efficienza geometrica in funzione delle posizioni delle lenti x1 e x2, filtrando i valori tra i percentili configurati in `config.json` (`lower_percentile`, `upper_percentile`).

```bash
./build/analysis/Release/plot2D
# Output: output/efficiency2D.png
```

### beam\_scan\_plot

**Legge**: `output/lens_simulation/lens.root`

Per una coppia (x1, x2) fissata, calcola la posizione media e la deviazione standard dei fotoni sul fotocatodo al variare della posizione della sorgente. Produce un grafico Y e uno Z con barre d'errore. Applica un filtro outlier a 2σ iterativo.

```bash
./build/analysis/Release/beam_scan_plot <x1> <x2> [percorso_file.root]

# Esempio:
./build/analysis/Release/beam_scan_plot 94.9 186.4
# Output: output/lens_simulation/beam_x1_94.90_x2_186.40.png
```

### m\_c\_creator

**Legge**: `output/lens_simulation/lens.root`

Per una terna (x1, x2, y0) fissata, raccoglie tutte le hit del run corrispondente e genera un istogramma 2D (y vs z) sul piano del fotocatodo. Applica un filtro outlier a 2σ iterativo per rimuovere fotoni diffusi.

```bash
./build/analysis/Release/m_c_creator <x1_target> <x2_target> <y0_target>

# Esempio:
./build/analysis/Release/m_c_creator 94.9 186.4 5.0
# Output: output/mean_covariance_maps/detector_hits_config_<id>_y0_5.0.png
```

---

## File di configurazione

`config/config.json` — letto da entrambi i programmi di simulazione e da `plot2D`:

```json
{
  "x_min": 33.0,       // Posizione minima ammessa per x1 [mm]
  "x_max": 171.0,      // Posizione massima ammessa per x2 [mm]
  "dx": 3.0,           // Passo della scansione [mm]
  "r1": 38.6,          // Raggio della lente 75mm [mm]
  "h1": 12.5,          // Spessore della lente 75mm [mm]
  "r2": 30.9,          // Raggio della lente 60mm [mm]
  "h2": 16.3,          // Spessore della lente 60mm [mm]
  "lower_percentile":  0.45,  // Soglia bassa per plot2D (esclude il 45% inferiore)
  "upper_percentile":  0.0    // Soglia alta per plot2D (0 = nessun taglio superiore)
}
```

I vincoli geometrici applicati durante la scansione sono:

```
x1_min = x_min - r1 + h1          (lente 75mm non esce a sinistra)
x1_max = x_max - h2 - 1 - r1      (lente 75mm non collide con lente 60mm)
x2_min = x1 + r1 + r2 + 1         (gap minimo di 1mm tra le lenti)
x2_max = x_max + r2 - h2          (lente 60mm non esce a destra)
```

Con i valori di default (dx=3 mm), la scansione genera **~703 configurazioni** e **~71.000 run** totali per `lens_simulation`.

---

## Formato dei file ROOT

### events.root (optimization)

```
events.root
└── TTree "events"    (~14M righe stimate con efficienza 20%)
    ├── x1        Double_t   posizione lente 75mm [mm]
    ├── x2        Double_t   posizione lente 60mm [mm]
    └── config_id Int_t      indice configurazione
```

### lens.root (lens_simulation)

```
lens.root
├── TTree "Configurations"    (703 righe)
│   ├── config_id  Int_t
│   ├── x1         Double_t  [mm]
│   └── x2         Double_t  [mm]
│
├── TTree "Runs"              (~71.000 righe)
│   ├── run_id     Int_t
│   ├── config_id  Int_t
│   ├── x_source   Float_t   posizione Y sorgente [mm]
│   └── n_hits     Int_t     fotoni rilevati in questo run
│
└── TTree "Hits"              (~14M righe stimate)
    ├── y_hit      Float_t   coordinata Y sul fotocatodo [mm]
    └── z_hit      Float_t   coordinata Z sul fotocatodo [mm]
```

La compressione è impostata con `SetCompressionLevel(404)` (LZ4 livello 4, se supportato dalla versione di Geant4/ROOT installata).

---

## Workflow completo

```bash
# 1. Build
cmake -S . -B build/ -G "Ninja Multi-Config"
cmake --build build/ --config Release

# 2a. Scansione efficienza grezza (optimization)
./build/Release/optimization_main -g geometry/main.gdml -o
#    → output/events.root

# 2b. Analisi mappa 2D efficienza
./build/analysis/Release/plot2D
#    → output/efficiency2D.png

# 3a. Beam scan dettagliato sulle configurazioni di interesse (lens_simulation)
./build/Release/lens_simulation_main -g geometry/main.gdml -l
#    → output/lens_simulation/lens.root

# 3b. Grafico risposta sistema ottico per una configurazione
./build/analysis/Release/beam_scan_plot 94.9 186.4
#    → output/lens_simulation/beam_x1_94.90_x2_186.40.png

# 3c. Istogramma 2D hit per una configurazione e posizione sorgente specifica
./build/analysis/Release/m_c_creator 94.9 186.4 5.0
#    → output/mean_covariance_maps/detector_hits_config_<id>_y0_5.0.png

# 4. Visualizzazione interattiva della geometria
./build/Release/optimization_main -g geometry/main.gdml -v
./build/Release/lens_simulation_main -g geometry/main.gdml -v
```