#!/bin/bash
# =============================================================================
# monitor.sh — Monitoraggio in tempo reale dei chunk paralleli
# Uso: ./monitor.sh <tmpdir> <n_jobs> <total_chunk_0> <total_chunk_1> ...
# =============================================================================

TMPDIR="$1"
N_JOBS="$2"
shift 2
TOTALS=("$@")

REFRESH=2
START_TIME=$(date +%s)

format_time() {
  local secs=$1
  printf "%02d:%02d:%02d" $(( secs/3600 )) $(( (secs%3600)/60 )) $(( secs%60 ))
}

progress_bar() {
  local pct=$1
  local width=30
  local filled=$(( pct * width / 100 ))
  local empty=$(( width - filled ))
  local bar="[" i
  for (( i=0; i<filled; i++ )); do bar="${bar}#"; done
  for (( i=0; i<empty;  i++ )); do bar="${bar}-"; done
  echo -n "${bar}]"
}

tput civis 2>/dev/null
trap 'tput cnorm 2>/dev/null; echo' EXIT

PREV_LINES=0

while true; do
  NOW=$(date +%s)
  ELAPSED=$(( NOW - START_TIME ))

  COMPLETED=()
  CHUNK_DONE_FLAGS=()
  ALL_DONE=true
  TOTAL_COMPLETED=0
  TOTAL_EXPECTED=0

  for (( i=0; i<N_JOBS; i++ )); do
    LOG="$TMPDIR/chunk_${i}.log"

    # Totale atteso: usa il valore scritto da lens_scan nel log se disponibile,
    # altrimenti usa la stima passata come argomento
    EXPECTED="${TOTALS[$i]}"
    EXPECTED=${EXPECTED//[$'\t\r\n ']/}
    [[ -z "$EXPECTED" || ! "$EXPECTED" =~ ^[0-9]+$ ]] && EXPECTED=0

    if [[ -f "$LOG" ]]; then
      REAL_TOTAL=$(grep "Total runs completed:" "$LOG" 2>/dev/null | tail -1 | grep -oP '\d+$')
      REAL_TOTAL=${REAL_TOTAL//[$'\t\r\n ']/}
      if [[ -n "$REAL_TOTAL" && "$REAL_TOTAL" =~ ^[0-9]+$ ]]; then
        EXPECTED=$REAL_TOTAL
      fi
    fi

    TOTAL_EXPECTED=$(( TOTAL_EXPECTED + EXPECTED ))

    if [[ ! -f "$LOG" ]]; then
      COMPLETED+=( 0 )
      CHUNK_DONE_FLAGS+=( false )
      ALL_DONE=false
      continue
    fi

    DONE=$(grep -c "Run done:" "$LOG" 2>/dev/null)
    DONE=${DONE//[$'\t\r\n ']/}
    [[ -z "$DONE" || ! "$DONE" =~ ^[0-9]+$ ]] && DONE=0
    COMPLETED+=( "$DONE" )
    TOTAL_COMPLETED=$(( TOTAL_COMPLETED + DONE ))

    # Chunk finito = "Optimization completed" nel log
    if grep -q "Optimization completed" "$LOG" 2>/dev/null; then
      CHUNK_DONE_FLAGS+=( true )
    else
      CHUNK_DONE_FLAGS+=( false )
      ALL_DONE=false
    fi
  done

  # ETA
  GLOBAL_PCT=0
  ETA_STR="--:--:--"
  if [[ $TOTAL_EXPECTED -gt 0 ]]; then
    GLOBAL_PCT=$(( TOTAL_COMPLETED * 100 / TOTAL_EXPECTED ))
    if [[ $TOTAL_COMPLETED -gt 0 && $ELAPSED -gt 0 ]]; then
      REMAINING=$(( TOTAL_EXPECTED - TOTAL_COMPLETED ))
      ETA_SECS=$(( REMAINING * ELAPSED / TOTAL_COMPLETED ))
      ETA_STR=$(format_time $ETA_SECS)
    fi
  fi

  # Ridisegna
  if [[ $PREV_LINES -gt 0 ]]; then
    for (( l=0; l<PREV_LINES; l++ )); do
      tput cuu1 2>/dev/null
      tput el   2>/dev/null
    done
  fi

  LINE_COUNT=0

  echo "=== riptide lens_simulation - progresso parallelo ==="
  (( LINE_COUNT++ ))
  echo "  Tempo trascorso : $(format_time $ELAPSED)   ETA: $ETA_STR"
  (( LINE_COUNT++ ))
  echo "  Run completati  : $TOTAL_COMPLETED / $TOTAL_EXPECTED"
  (( LINE_COUNT++ ))
  echo -n "  Progresso globale: $(progress_bar $GLOBAL_PCT) ${GLOBAL_PCT}%"
  echo ""
  (( LINE_COUNT++ ))
  echo ""
  (( LINE_COUNT++ ))

  for (( i=0; i<N_JOBS; i++ )); do
    DONE="${COMPLETED[$i]}"
    EXP="${TOTALS[$i]}"
    EXP=${EXP//[$'\t\r\n ']/}
    [[ -z "$EXP" || ! "$EXP" =~ ^[0-9]+$ ]] && EXP=1

    # Aggiorna EXP con il valore reale dal log se disponibile
    LOG="$TMPDIR/chunk_${i}.log"
    if [[ -f "$LOG" ]]; then
      REAL_TOTAL=$(grep "Total runs completed:" "$LOG" 2>/dev/null | tail -1 | grep -oP '\d+$')
      REAL_TOTAL=${REAL_TOTAL//[$'\t\r\n ']/}
      if [[ -n "$REAL_TOTAL" && "$REAL_TOTAL" =~ ^[0-9]+$ && "$REAL_TOTAL" -gt 0 ]]; then
        EXP=$REAL_TOTAL
      fi
    fi

    PCT=0
    [[ $EXP -gt 0 ]] && PCT=$(( DONE * 100 / EXP ))
    [[ $PCT -gt 100 ]] && PCT=100

    if ${CHUNK_DONE_FLAGS[$i]}; then
      STATUS="[DONE]"
      PCT=100
    elif [[ ! -f "$LOG" ]]; then
      STATUS="[wait]"
    else
      STATUS="[ >> ]"
    fi

    printf "  Chunk %d %s $(progress_bar $PCT) %3d%%  (%d/%d)\n" \
      "$i" "$STATUS" "$PCT" "$DONE" "$EXP"
    (( LINE_COUNT++ ))
  done

  PREV_LINES=$LINE_COUNT

  if $ALL_DONE && [[ $N_JOBS -gt 0 ]]; then
    echo ""
    echo "COMPLETATO in $(format_time $ELAPSED)"
    break
  fi

  sleep "$REFRESH"
done

tput cnorm 2>/dev/null
