#!/bin/bash
# =============================================================================
# run.sh — Lancia lens_simulation o optimization, locale o su SSD esterna
#          Supporta la parallelizzazione a processi con chunking automatico
#
# Uso:
#   ./run.sh [lens|opt] [local|ssd] [--jobs N] [opzioni extra passate al binario]
#
# Esempi:
#   ./run.sh lens local
#   ./run.sh lens ssd
#   ./run.sh lens ssd --jobs 4
#   ./run.sh opt  ssd --jobs 4 --ssd-mount /mnt/myusb
#   ./run.sh lens local --output /tmp/test.root --config config/config.json
# =============================================================================

set -e

# Root del progetto (cartella padre di scripts/)
PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# configurazione di default
BUILD_DIR="$PROJECT_ROOT/build"
GEOMETRY="$PROJECT_ROOT/geometry/main.gdml"
SSD_MOUNT="/mnt/external_ssd"
BINARY_LENS="$BUILD_DIR/Release/lens_simulation_main"
BINARY_OPT="$BUILD_DIR/Release/optimization_main"
N_JOBS=1
CONFIG_FILE="$PROJECT_ROOT/config/config.json"
LENS75_ID=""
LENS60_ID=""
ALL_LENSES=0

# parsing argomenti posizionali
MODE="${1:-lens}"    # lens | opt
TARGET="${2:-local}" # local | ssd
shift 2 || true

# parsing --jobs e --ssd-mount prima di passare il resto al binario
EXTRA_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --jobs)
      N_JOBS="$2"; shift 2 ;;
    --ssd-mount)
      SSD_MOUNT="$2"; shift 2 ;;
    --all-lenses)
      ALL_LENSES=1; EXTRA_ARGS+=("--all-lenses"); shift ;;
    --config)
      CONFIG_FILE="$2"; EXTRA_ARGS+=("--config" "$2"); shift 2 ;;
    --lens75-id)
      LENS75_ID="$2"; EXTRA_ARGS+=("--lens75-id" "$2"); shift 2 ;;
    --lens60-id)
      LENS60_ID="$2"; EXTRA_ARGS+=("--lens60-id" "$2"); shift 2 ;;
    *)
      EXTRA_ARGS+=("$1"); shift ;;
  esac
done

# verifica opzioni
if [[ "$MODE" == "lens" && "$ALL_LENSES" -eq 1 ]]; then
  echo "[ERROR] L'opzione --all-lenses non è supportata per lens_simulation." >&2
  echo "        Usala solo con 'opt' per l'ottimizzazione di tutte le combinazioni." >&2
  exit 1
fi

# verifica build
BINARY="$BINARY_LENS"
[[ "$MODE" == "opt" ]] && BINARY="$BINARY_OPT"
if [[ ! -f "$BINARY" ]]; then
  echo "[ERROR] Binario non trovato: $BINARY — esegui cmake + make prima." >&2
  exit 1
fi

# verifica SSD
SSD_ARGS=""
if [[ "$TARGET" == "ssd" ]]; then
  if ! mountpoint -q "$SSD_MOUNT"; then
    echo "[ERROR] $SSD_MOUNT non è un mountpoint attivo." >&2
    echo "        sudo mount -o noatime,nodiratime,discard /dev/nvme1n1 $SSD_MOUNT" >&2
    exit 1
  fi
  echo "[INFO]  SSD rilevata su $SSD_MOUNT"
  SSD_ARGS="--ssd --ssd-mount $SSD_MOUNT"
fi

# run singolo (N_JOBS=1)
if [[ "$N_JOBS" -eq 1 ]]; then
  echo "[INFO]  Avvio $MODE ($TARGET), processo singolo"
  case "$MODE" in
    lens) "$BINARY" -g "$GEOMETRY" -b -l --config "$CONFIG_FILE" $SSD_ARGS "${EXTRA_ARGS[@]}" ;;
    opt)  "$BINARY" -g "$GEOMETRY" -b -o --config "$CONFIG_FILE" $SSD_ARGS "${EXTRA_ARGS[@]}" ;;
  esac
  echo "[DONE]  Simulazione completata."
  exit 0
fi

# run parallelo (N_JOBS > 1)
echo "[INFO]  Avvio $MODE ($TARGET) con $N_JOBS processi paralleli"

# Genera tutte le coppie (x1, x2) valide e le distribuisce in N chunk bilanciati
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
TMPDIR_CONFIGS="/tmp/riptide_chunks_$TIMESTAMP"
mkdir -p "$TMPDIR_CONFIGS"

# --- Gestione Segnali e Pulizia ---
PIDS=()
CHUNK_OUTPUTS=()
MONITOR_PID=""

cleanup() {
  echo -e "\n\n[WARN]  Interruzione rilevata. Terminazione processi in corso..."
  
  # Ferma il monitor/dashboard
  if [[ -n "$MONITOR_PID" ]]; then
    kill "$MONITOR_PID" 2>/dev/null || true
  fi

  # Ferma tutti i chunk di simulazione
  for pid in "${PIDS[@]}"; do
    kill "$pid" 2>/dev/null || true
  done

  # Rimuovi file temporanei e chunk parziali
  echo "[INFO]  Pulizia file temporanei..."
  rm -rf "$TMPDIR_CONFIGS"
  for f in "${CHUNK_OUTPUTS[@]}"; do
    rm -f "$f" 2>/dev/null || true
  done

  echo "[DONE]  Pulizia completata. Uscita."
  exit 1
}

# Registra il trap per SIGINT (Ctrl+C) e SIGTERM
trap cleanup SIGINT SIGTERM
# ----------------------------------

# Python: calcola tutte le coppie, le divide in N chunk con lo stesso numero di coppie
python3 - <<EOF
import json, numpy as np, os

with open('$CONFIG_FILE') as f:
    c = json.load(f)

x_min = c['x_min']; x_max = c['x_max']; dx = c['dx']
n_jobs = $N_JOBS
mode = '$MODE'
lens75_id = '$LENS75_ID'
lens60_id = '$LENS60_ID'
all_lenses = $ALL_LENSES

# Load Thorlabs data
thorlabs_data = {}
n_models = 1

def load_tsv(filename, has_rotation=False):
    path = os.path.join('$PROJECT_ROOT', 'lens_cutter/lens_data', filename)
    if os.path.exists(path):
        with open(path) as f:
            lines = f.readlines()[1:] # skip header
            for line in lines:
                parts = line.split('\t')
                if len(parts) >= 6:
                    lid = parts[0]
                    tc = float(parts[4])
                    te = float(parts[5])
                    rot = float(parts[7]) if has_rotation and len(parts) >= 8 else 0.0
                    # Offset center: (tc - te)/2 * sign
                    # rotation 0 -> Z->X, sign +1
                    # rotation 180 -> Z->-X, sign -1
                    sign = -1.0 if abs(rot - 180.0) < 1e-6 else 1.0
                    offset = (tc - te) / 2.0 * sign if has_rotation else 0.0
                    thorlabs_data[lid] = {'h': tc, 'offset': offset}

load_tsv('thorlabs_biconvex.tsv', has_rotation=False)
load_tsv('thorlabs_planoconvex.tsv', has_rotation=True)

if mode == 'opt' and all_lenses:
    n_models = len(thorlabs_data) * len(thorlabs_data)

if mode == 'lens':
    # Parametri sorgente GPS (mappa 3D PSF)
    s_x_min = c.get('source_x_min', -30.0)
    s_x_max = c.get('source_x_max', 30.0)
    s_dx    = c.get('source_dx', 5.0)
    
    s_y_min = c.get('source_y_min', 0.0)
    s_y_max = c.get('source_y_max', 10.0 * np.sqrt(2.0))
    s_dy    = c.get('source_dy', 1.0)

    nx = len(np.arange(s_x_min, s_x_max + 1e-9, s_dx))
    ny = len(np.arange(s_y_min, s_y_max + 1e-9, s_dy))
    n_runs_per_pair = nx * ny
else:
    n_runs_per_pair = 1    # una sola esecuzione per configurazione

# Genera tutte le coppie valide (Adaptive Loop)
pairs = []
# Default GDML lens parameters (valori di backup se non trovati nel TSV)
h1 = 12.5
h2 = 16.3
margin = 3.0 if mode == 'lens' else 1.0

# Carica spessori reali se disponibili
if lens75_id in thorlabs_data:
    h1 = thorlabs_data[lens75_id]['h']
if lens60_id in thorlabs_data:
    h2 = thorlabs_data[lens60_id]['h']

for x1 in np.arange(x_min, x_max + 1e-9, dx):
    # x2_min per evitare collisioni: x1 + h1/2 + margin < x2 - h2/2
    # Poiché x1 e x2 sono i centri geometrici, la collisione è solo sulla somma dei semi-spessori.
    x2_min_collision = x1 + (h1 + h2) / 2.0 + margin
    # Align x2_start to the dx grid relative to x_min
    x2_start_raw = max(x1 + dx, x2_min_collision)
    x2_start = x_min + np.ceil(round((x2_start_raw - x_min) / dx, 8)) * dx
    
    for x2 in np.arange(x2_start, x_max + 1e-9, dx):
        pairs.append((round(float(x1), 6), round(float(x2), 6)))

total_pairs = len(pairs)
chunk_size  = (total_pairs + n_jobs - 1) // n_jobs

print(f"[INFO]  Coppie (x1,x2) totali: {total_pairs}, ~{chunk_size} per chunk")

config_offset = 0
run_offset    = 0

for chunk in range(n_jobs):
    start = chunk * chunk_size
    end   = min(start + chunk_size, total_pairs)
    if start >= total_pairs:
        break

    chunk_pairs = pairs[start:end]
    chunk_runs  = len(chunk_pairs) * n_runs_per_pair * n_models

    chunk_config = dict(c)
    chunk_config['pairs']            = chunk_pairs
    chunk_config['config_id_offset'] = config_offset
    chunk_config['run_id_offset']    = run_offset

    path = f"$TMPDIR_CONFIGS/config_chunk_{chunk}.json"
    with open(path, 'w') as f:
        json.dump(chunk_config, f, indent=2)

    # Scrivi info chunk su file per bash
    with open(f"$TMPDIR_CONFIGS/chunk_{chunk}_info.txt", 'w') as f:
        f.write(f"{chunk_runs}\n")

    print(f"[INFO]  Chunk {chunk}: {len(chunk_pairs)} coppie, {chunk_runs} run, config_id_offset={config_offset}")

    config_offset += len(chunk_pairs) * n_models
    run_offset    += chunk_runs
EOF

# Conta i chunk effettivamente generati da Python
ACTUAL_JOBS=$(ls "$TMPDIR_CONFIGS"/config_chunk_*.json 2>/dev/null | wc -l)

CHUNK_TOTALS=()

for (( chunk=0; chunk<ACTUAL_JOBS; chunk++ )); do
  CHUNK_CONFIG="$TMPDIR_CONFIGS/config_chunk_${chunk}.json"

  if [[ "$MODE" == "opt" ]]; then
    SUBDIR="optimization"
    EXT="events"
  else
    SUBDIR="lens_simulation"
    EXT="lens"
  fi

  if [[ "$TARGET" == "ssd" ]]; then
    CHUNK_OUTPUT="$SSD_MOUNT/riptide/runs/run_${TIMESTAMP}/chunk_${chunk}.root"
  else
    CHUNK_OUTPUT="$PROJECT_ROOT/output/${SUBDIR}/chunk_${TIMESTAMP}_${chunk}.root"
  fi
  CHUNK_OUTPUTS+=("$CHUNK_OUTPUT")

  CHUNK_TOTAL=$(cat "$TMPDIR_CONFIGS/chunk_${chunk}_info.txt" 2>/dev/null | tr -d '[:space:]')
  [[ -z "$CHUNK_TOTAL" || ! "$CHUNK_TOTAL" =~ ^[0-9]+$ ]] && CHUNK_TOTAL=1
  CHUNK_TOTALS+=("$CHUNK_TOTAL")

  # Esegue il binario corretto con i flag corretti
  LOG="$TMPDIR_CONFIGS/chunk_${chunk}.log"
  # In parallelo non passiamo mai --ssd al binario perché run.sh gestisce già il path di output.
  # Passiamo solo --ssd-mount se necessario.
  PARALLEL_SSD_ARGS="${SSD_ARGS/--ssd/}"

  if [[ "$MODE" == "opt" ]]; then
    # In optimization mode, we use -o
    "$BINARY" -g "$GEOMETRY" -b -o --config "$CHUNK_CONFIG" --output "$CHUNK_OUTPUT" $PARALLEL_SSD_ARGS "${EXTRA_ARGS[@]}" > "$LOG" 2>&1 &
  else
    # In lens simulation mode, we use -l
    "$BINARY" -g "$GEOMETRY" -b -l --config "$CHUNK_CONFIG" --output "$CHUNK_OUTPUT" $PARALLEL_SSD_ARGS "${EXTRA_ARGS[@]}" > "$LOG" 2>&1 &
  fi
  PIDS+=($!)
done

# Avvia il monitor in foreground su /dev/tty (terminale diretto)
# I chunk girano in background, il monitor prende il controllo del display
if [[ -f "$(dirname "$0")/dashboard.py" ]]; then
  chmod +x "$(dirname "$0")/dashboard.py"
  "$(dirname "$0")/dashboard.py" "$TMPDIR_CONFIGS" "$ACTUAL_JOBS" "${CHUNK_TOTALS[@]}" </dev/null &
  MONITOR_PID=$!
elif [[ -f "$(dirname "$0")/monitor.sh" ]]; then
  chmod +x "$(dirname "$0")/monitor.sh"
  "$(dirname "$0")/monitor.sh" "$TMPDIR_CONFIGS" "$ACTUAL_JOBS" "${CHUNK_TOTALS[@]}" </dev/null &
  MONITOR_PID=$!
fi

# Attendi tutti i processi e controlla gli exit code
echo "[INFO]  Attesa completamento di ${#PIDS[@]} processi..."
FAILED=0
for i in "${!PIDS[@]}"; do
  wait "${PIDS[$i]}" || { echo "[ERROR] Chunk $i fallito (vedi $TMPDIR_CONFIGS/chunk_${i}.log)"; FAILED=1; }
  echo "[INFO]  Chunk $i completato"
done

# Ferma il monitor
[[ -n "$MONITOR_PID" ]] && wait "$MONITOR_PID" 2>/dev/null || true

[[ $FAILED -eq 1 ]] && { echo "[ERROR] Uno o più chunk sono falliti. Merge annullato."; exit 1; }

# Merge con hadd
if [[ "$MODE" == "opt" ]]; then
  OUT_NAME="events"
  OUT_DIR="optimization"
else
  OUT_NAME="lens"
  OUT_DIR="lens_simulation"
fi

if [[ "$TARGET" == "ssd" ]]; then
  MERGED_OUTPUT="$SSD_MOUNT/riptide/runs/run_${TIMESTAMP}/${OUT_NAME}.root"
else
  MERGED_OUTPUT="output/${OUT_DIR}/${OUT_NAME}_${TIMESTAMP}.root"
fi

echo "[INFO]  Merge dei chunk in $MERGED_OUTPUT"
hadd -fk "$MERGED_OUTPUT" "${CHUNK_OUTPUTS[@]}"

# Rimuovi i chunk intermedi
echo "[INFO]  Rimozione chunk intermedi..."
rm -f "${CHUNK_OUTPUTS[@]}"
rm -rf "$TMPDIR_CONFIGS"

echo "[DONE]  Output finale: $MERGED_OUTPUT"
