#!/bin/bash
set -e

L1=${1:-CHOS60R}
L2=${2:-CHOS25R}
BASE="output/lens_simulations/${L1}_${L2}"
DATA="$BASE/data"
PLOTS="$BASE/plots"
PARETO_TSV="$BASE/pareto/balanced/pareto_results.tsv"

if [ ! -f "$PARETO_TSV" ]; then
  echo "Errore: $PARETO_TSV non trovato. Eseguire prima pareto_runner.py." >&2
  exit 1
fi

read -r X1 X2 <<< "$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next} $c["pareto_rank"]==1{print $c["x1"], $c["x2"]}' "$PARETO_TSV")"

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
