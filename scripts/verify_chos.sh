#!/bin/bash
set -e

L1=${1:-CHOS60R}
L2=${2:-CHOS25R}
BASE="output/lens_simulations/${L1}_${L2}"
DATA="$BASE/data"
PLOTS="$BASE/plots"

if [ ! -f "$BASE/simulation_summary.json" ]; then
  echo "Errore: $BASE/simulation_summary.json non trovato. Eseguire prima lens_runner.py." >&2
  exit 1
fi

X1=$(python3 -c "import json; d=json.load(open('$BASE/simulation_summary.json')); print(d['best_x1'])")
X2=$(python3 -c "import json; d=json.load(open('$BASE/simulation_summary.json')); print(d['best_x2'])")

echo "=== Verifica $L1 + $L2 (x1=$X1 mm, x2=$X2 mm) ==="

echo "→ Efficiency map"
./build/analysis/Release/plot2D \
  --input "$DATA/events.root" \
  --output "$PLOTS/verify_eff2D.png"

echo "→ PSF trace @ best focus"
./build/analysis/psf_analysis/Release/trace_viewer \
  --psf "$DATA/psf_data.root" \
  --x1 "$X1" --x2 "$X2" --y0 5.0 --fit \
  --output "$PLOTS/verify_trace.png"

echo "→ Q-map (log scale)"
./build/analysis/psf_analysis/Release/q_map \
  --psf "$DATA/psf_data.root" --log \
  --output "$PLOTS/verify_qmap.png"

echo "Done — plots in $PLOTS/"
