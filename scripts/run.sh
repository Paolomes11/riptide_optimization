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
    --config)
      CONFIG_FILE="$2"; EXTRA_ARGS+=("--config" "$2"); shift 2 ;;
    *)
      EXTRA_ARGS+=("$1"); shift ;;
  esac
done

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
    lens) "$BINARY" -g "$GEOMETRY" -b -l $SSD_ARGS "${EXTRA_ARGS[@]}" ;;
    opt)  "$BINARY" -g "$GEOMETRY" -b -o $SSD_ARGS "${EXTRA_ARGS[@]}" ;;
  esac
  echo "[DONE]  Simulazione completata."
  exit 0
fi

# run parallelo (N_JOBS > 1) — solo lens_simulation
if [[ "$MODE" != "lens" ]]; then
  echo "[ERROR] La parallelizzazione --jobs è supportata solo per 'lens'." >&2
  exit 1
fi

echo "[INFO]  Avvio lens_simulation ($TARGET) con $N_JOBS processi paralleli"

# Genera tutte le coppie (x1, x2) valide e le distribuisce in N chunk bilanciati
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
TMPDIR_CONFIGS="/tmp/riptide_chunks_$TIMESTAMP"
mkdir -p "$TMPDIR_CONFIGS"

# Python: calcola tutte le coppie, le divide in N chunk con lo stesso numero di coppie
python3 - <<EOF
import json, numpy as np, os

with open('$CONFIG_FILE') as f:
    c = json.load(f)

x_min = c['x_min']; x_max = c['x_max']; dx = c['dx']
r1 = c['r1']; h1 = c['h1']; r2 = c['r2']; h2 = c['h2']
n_jobs = $N_JOBS
n_runs_per_pair = 101  # y_source 0..10 step 0.1

# Genera tutte le coppie valide
pairs = []
for x1 in np.arange(x_min - r1 + h1, x_max - h2 - 1 - r1 + 1e-9, dx):
    for x2 in np.arange(x1 + r1 + r2 + 1, x_max + r2 - h2 + 1e-9, dx):
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
    chunk_runs  = len(chunk_pairs) * n_runs_per_pair

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

    config_offset += len(chunk_pairs)
    run_offset    += chunk_runs
EOF

# Conta i chunk effettivamente generati da Python
ACTUAL_JOBS=$(ls "$TMPDIR_CONFIGS"/config_chunk_*.json 2>/dev/null | wc -l)

PIDS=()
CHUNK_OUTPUTS=()
CHUNK_TOTALS=()

for (( chunk=0; chunk<ACTUAL_JOBS; chunk++ )); do
  CHUNK_CONFIG="$TMPDIR_CONFIGS/config_chunk_${chunk}.json"

  if [[ "$TARGET" == "ssd" ]]; then
    CHUNK_OUTPUT="$SSD_MOUNT/riptide/runs/run_${TIMESTAMP}/chunk_${chunk}.root"
  else
    CHUNK_OUTPUT="$PROJECT_ROOT/output/lens_simulation/chunk_${TIMESTAMP}_${chunk}.root"
  fi
  CHUNK_OUTPUTS+=("$CHUNK_OUTPUT")

  CHUNK_TOTAL=$(cat "$TMPDIR_CONFIGS/chunk_${chunk}_info.txt" 2>/dev/null | tr -d '[:space:]')
  [[ -z "$CHUNK_TOTAL" || ! "$CHUNK_TOTAL" =~ ^[0-9]+$ ]] && CHUNK_TOTAL=1
  CHUNK_TOTALS+=("$CHUNK_TOTAL")

  "$BINARY_LENS" -g "$GEOMETRY" -b -l \
    --config "$CHUNK_CONFIG" \
    --output "$CHUNK_OUTPUT" \
    "${EXTRA_ARGS[@]}" \
    > "$TMPDIR_CONFIGS/chunk_${chunk}.log" 2>&1 &
  PIDS+=($!)
done

# Avvia il monitor in foreground su /dev/tty (terminale diretto)
# I chunk girano in background, il monitor prende il controllo del display
MONITOR_PID=""
if [[ -f "$(dirname "$0")/monitor.sh" ]]; then
  chmod +x "$(dirname "$0")/monitor.sh"
  "$(dirname "$0")/monitor.sh" "$TMPDIR_CONFIGS" "$N_JOBS" "${CHUNK_TOTALS[@]}" </dev/null &
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
if [[ "$TARGET" == "ssd" ]]; then
  MERGED_OUTPUT="$SSD_MOUNT/riptide/runs/run_${TIMESTAMP}/lens.root"
else
  MERGED_OUTPUT="output/lens_simulation/lens_${TIMESTAMP}.root"
fi

echo "[INFO]  Merge dei chunk in $MERGED_OUTPUT"
hadd -fk "$MERGED_OUTPUT" "${CHUNK_OUTPUTS[@]}"

# Rimuovi i chunk intermedi
echo "[INFO]  Rimozione chunk intermedi..."
rm -f "${CHUNK_OUTPUTS[@]}"
rm -rf "$TMPDIR_CONFIGS"

echo "[DONE]  Output finale: $MERGED_OUTPUT"
