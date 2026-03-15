#!/bin/bash
# =============================================================================
# run.sh — Lancia lens_simulation o optimization, locale o su SSD esterna
#
# Uso:
#   ./run.sh [lens|opt] [local|ssd] [opzioni extra passate al binario]
#
# Esempi:
#   ./run.sh lens local
#   ./run.sh lens ssd
#   ./run.sh opt  ssd  --ssd-mount /mnt/myusb
#   ./run.sh lens local --output /tmp/test.root --config config/config_alt.json
# =============================================================================

set -e

# ---------- configurazione di default ----------
BUILD_DIR="build"
GEOMETRY="geometry/main.gdml"
SSD_MOUNT="/mnt/external_ssd"
BINARY_LENS="$BUILD_DIR/Release/lens_simulation_main"
BINARY_OPT="$BUILD_DIR/Release/optimization_main"

# ---------- parsing argomenti posizionali ----------
MODE="${1:-lens}"   # lens | opt
TARGET="${2:-local}" # local | ssd
shift 2 || true      # il resto va direttamente ai binari

# ---------- verifica build ----------
if [[ "$MODE" == "lens" && ! -f "$BINARY_LENS" ]]; then
  echo "[ERROR] Binario non trovato: $BINARY_LENS — esegui cmake + make prima." >&2
  exit 1
fi
if [[ "$MODE" == "opt" && ! -f "$BINARY_OPT" ]]; then
  echo "[ERROR] Binario non trovato: $BINARY_OPT — esegui cmake + make prima." >&2
  exit 1
fi

# ---------- modalità SSD ----------
SSD_ARGS=""
if [[ "$TARGET" == "ssd" ]]; then
  # Verifica che l'SSD sia montata
  if ! mountpoint -q "$SSD_MOUNT"; then
    echo "[ERROR] $SSD_MOUNT non è un mountpoint attivo." >&2
    echo "        Monta il device prima di procedere:" >&2
    echo "        sudo mount -o noatime,nodiratime,discard /dev/nvme1n1p1 $SSD_MOUNT" >&2
    exit 1
  fi
  echo "[INFO]  SSD rilevata su $SSD_MOUNT"
  SSD_ARGS="--ssd --ssd-mount $SSD_MOUNT"
fi

# ---------- lancio ----------
case "$MODE" in
  lens)
    echo "[INFO]  Avvio lens_simulation ($TARGET)"
    "$BINARY_LENS" -g "$GEOMETRY" -b -l $SSD_ARGS "$@"
    ;;
  opt)
    echo "[INFO]  Avvio optimization ($TARGET)"
    "$BINARY_OPT" -g "$GEOMETRY" -b -o $SSD_ARGS "$@"
    ;;
  *)
    echo "Uso: $0 [lens|opt] [local|ssd] [opzioni extra]" >&2
    exit 1
    ;;
esac

echo "[DONE]  Simulazione completata."