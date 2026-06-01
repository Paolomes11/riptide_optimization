#!/usr/bin/env python3
"""
pareto_runner.py — Post-processing Pareto per una coppia di lenti già simulata.

Legge config/pareto_runs.json con multipli set di pesi Pareto ed esegue
pareto_selector per ognuno. I risultati vanno in:
  output/lens_simulations/{l1_id}_{l2_id}/pareto/{label}/

I file di input (events.root, q_map.tsv, chi2_map.tsv, dof_map.tsv,
resolution_map.tsv) vengono letti da:
  output/lens_simulations/{l1_id}_{l2_id}/data/
(oppure da data_dir nel config, se specificato).

Richiede che lens_runner.py sia stato eseguito in precedenza per la coppia.
"""
from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path

# ─── Costanti ────────────────────────────────────────────────────────────────

PROJECT_ROOT = Path(__file__).resolve().parent.parent
BUILD        = PROJECT_ROOT / "build"

_PARETO_BIN  = BUILD / "analysis" / "pareto_analysis" / "Release" / "pareto_selector"
TIMEOUT_SEC  = 10 * 60


def find_pareto_binary() -> Path:
    if _PARETO_BIN.exists():
        return _PARETO_BIN
    alt = Path(str(_PARETO_BIN).replace("/Release/", "/"))
    if alt.exists():
        return alt
    raise FileNotFoundError(
        f"Binario 'pareto_selector' non trovato: {_PARETO_BIN}\n"
        f"  Alternativa cercata: {alt}\n"
        f"  Esegui il build su questa macchina."
    )


def run_cmd(cmd: list[str], log_path: Path, dry_run: bool = False) -> bool:
    cmd_str = " ".join(str(c) for c in cmd)
    logging.info(f"CMD: {cmd_str}")
    if dry_run:
        logging.info("[DRY-RUN] (non eseguito)")
        return True
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "a") as lf:
        lf.write(f"\n{'='*60}\nCMD: {cmd_str}\n{'='*60}\n")
        try:
            result = subprocess.run(
                [str(c) for c in cmd],
                stdout=lf, stderr=subprocess.STDOUT,
                timeout=TIMEOUT_SEC,
            )
            ok = result.returncode == 0
            if not ok:
                logging.error(f"FALLITO (rc={result.returncode}): {cmd_str}")
            return ok
        except subprocess.TimeoutExpired:
            logging.error(f"TIMEOUT ({TIMEOUT_SEC}s): {cmd_str}")
            return False
        except Exception as e:
            logging.error(f"ECCEZIONE: {e}")
            return False


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Pareto post-processing per una coppia di lenti RIPTIDE",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Esempi:
  # Usa config/pareto_runs.json (default)
  python3 scripts/pareto_runner.py

  # Override della coppia di lenti
  python3 scripts/pareto_runner.py --l1-id LA4078 --l2-id LA4464R

  # Config custom + dry-run
  python3 scripts/pareto_runner.py --config config/pareto_runs_custom.json --dry-run
""",
    )
    parser.add_argument(
        "--config", type=str,
        default=str(PROJECT_ROOT / "config" / "pareto_runs.json"),
        help="Config JSON con coppia di lenti e set di pesi (default: config/pareto_runs.json)",
    )
    parser.add_argument("--l1-id",  type=str, default=None,
                        help="Override l1_id dal config")
    parser.add_argument("--l2-id",  type=str, default=None,
                        help="Override l2_id dal config")
    parser.add_argument("--dry-run", action="store_true",
                        help="Stampa comandi senza eseguire")
    args = parser.parse_args()

    cfg_path = Path(args.config)
    if not cfg_path.exists():
        print(f"Errore: config non trovato: {cfg_path}", file=sys.stderr)
        print("Crea config/pareto_runs.json o specifica --config PATH", file=sys.stderr)
        return 1

    cfg   = json.loads(cfg_path.read_text())
    l1_id = args.l1_id or cfg.get("l1_id", "")
    l2_id = args.l2_id or cfg.get("l2_id", "")
    runs  = cfg.get("runs", [])

    if not l1_id or not l2_id:
        print("Errore: specificare l1_id e l2_id nel config o via --l1-id/--l2-id", file=sys.stderr)
        return 1
    if not runs:
        print("Errore: 'runs' è vuoto nel config", file=sys.stderr)
        return 1

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )

    pair_tag = f"{l1_id}_{l2_id}"
    data_dir_cfg = cfg.get("data_dir")
    data_dir = (
        Path(data_dir_cfg)
        if data_dir_cfg
        else PROJECT_ROOT / "output" / "lens_simulations" / pair_tag / "data"
    )
    pareto_base = PROJECT_ROOT / "output" / "lens_simulations" / pair_tag / "pareto"

    logging.info(f"=== Pareto Runner — {pair_tag} — {len(runs)} run ===")
    logging.info(f"Data dir : {data_dir}")
    logging.info(f"Output   : {pareto_base}")

    required = {
        "events.root":        data_dir / "events.root",
        "q_map.tsv":          data_dir / "q_map.tsv",
        "chi2_map.tsv":       data_dir / "chi2_map.tsv",
        "dof_map.tsv":        data_dir / "dof_map.tsv",
        "resolution_map.tsv": data_dir / "resolution_map.tsv",
    }
    if not args.dry_run:
        missing = [name for name, p in required.items() if not p.exists()]
        if missing:
            logging.error(f"File mancanti in {data_dir}: {missing}")
            logging.error("Esegui prima: python3 scripts/lens_runner.py --l1-id %s --l2-id %s", l1_id, l2_id)
            return 1

    pareto_bin = find_pareto_binary()
    n_ok = 0

    for run in runs:
        label = run.get("label", f"run{runs.index(run)}")
        logging.info(f"--- {label} ---")

        out_dir  = pareto_base / label
        out_dir.mkdir(parents=True, exist_ok=True)
        tsv_out  = out_dir / "pareto_results.tsv"
        png_out  = out_dir / "pareto_plot.png"
        log_path = out_dir / "pareto.log"

        cmd = [
            pareto_bin,
            "--events",     str(required["events.root"]),
            "--qmap",       str(required["q_map.tsv"]),
            "--chi2map",    str(required["chi2_map.tsv"]),
            "--dofmap",     str(required["dof_map.tsv"]),
            "--resolution", str(required["resolution_map.tsv"]),
            "--tsv",        str(tsv_out),
            "--ee80-max",   str(run.get("ee80_max", 10.0)),
            "--w-eta",      str(run.get("w_eta",  0.1)),
            "--w-Q",        str(run.get("w_Q",    0.1)),
            "--w-dof",      str(run.get("w_dof",  0.3)),
            "--w-M",        str(run.get("w_M",    0.5)),
            "--output",     str(png_out),
            "--l1-id",      l1_id,
            "--l2-id",      l2_id,
        ]

        ok = run_cmd(cmd, log_path, args.dry_run)
        if ok:
            logging.info(f"  TSV : {tsv_out}")
            logging.info(f"  PNG : {png_out}")
            n_ok += 1
        else:
            logging.error(f"  Run '{label}' FALLITA")

    status = "COMPLETATI" if n_ok == len(runs) else "PARZIALI"
    logging.info(f"{status}: {n_ok}/{len(runs)} run Pareto per {pair_tag}")
    return 0 if n_ok == len(runs) else 1


if __name__ == "__main__":
    sys.exit(main())
