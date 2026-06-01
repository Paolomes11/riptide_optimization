#!/usr/bin/env python3
"""
autonomous_optimizer.py — Driver autonomo per sweep parametrico RIPTIDE.

Esegue pipeline completa (opt → lens → psf_extractor → q_map → chi2_map →
dof → dof_map → psf-dof → resolution_map → pareto_selector) con strategia
a due fasi, registry di resume, validazione output e fix automatici.
"""
from __future__ import annotations

import argparse
import datetime
import json
import logging
import os
import shutil
import signal
import subprocess
import sys
import textwrap
import time
from pathlib import Path

# ─── Costanti ────────────────────────────────────────────────────────────────

PROJECT_ROOT = Path(__file__).resolve().parent.parent
BUILD        = PROJECT_ROOT / "build"

BINARIES = {
    "run_sh":             PROJECT_ROOT / "scripts" / "run.sh",
    "plot2D":             BUILD / "analysis" / "Release" / "plot2D",
    "psf_extractor":      BUILD / "analysis" / "Release" / "psf_extractor",
    "psf_dof_extractor":  BUILD / "analysis" / "dof_analysis" / "Release" / "psf_dof_extractor",
    "q_map":              BUILD / "analysis" / "psf_analysis" / "Release" / "q_map",
    "chi2_map":           BUILD / "analysis" / "psf_analysis" / "Release" / "chi2_map",
    "dof_map":            BUILD / "analysis" / "dof_analysis" / "Release" / "dof_map",
    "resolution_map":     BUILD / "analysis" / "resolution_analysis" / "Release" / "resolution_map",
    "pareto_selector":    BUILD / "analysis" / "pareto_analysis" / "Release" / "pareto_selector",
}

# Timeout in secondi per ogni step
TIMEOUTS_SEC = {
    "opt":               60  * 60,
    "lens":              210 * 60,
    "dof":               90  * 60,
    "psf-dof":           90  * 60,
    "resolution_map":    30  * 60,
    "psf_dof_extractor": 30  * 60,
    "analysis":          10  * 60,   # psf_extractor, q_map, chi2_map, dof_map, pareto_selector
}

SWEEP_GEOM = [
    {"name": "nominal",  "source_dx": 0.5, "source_dy": 0.5,
     "source_x_min": -30.0, "source_x_max": 30.0, "source_y_max": 14.14},
    {"name": "coarse",   "source_dx": 1.0, "source_dy": 1.0,
     "source_x_min": -30.0, "source_x_max": 30.0, "source_y_max": 14.14},
    {"name": "extended", "source_dx": 1.0, "source_dy": 1.0,
     "source_x_min": -45.0, "source_x_max": 45.0, "source_y_max": 20.0},
]
MARGIN_VALUES = [0.5, 1.0, 2.0]

# ─── Stato globale per cleanup su SIGINT ────────────────────────────────────

_active_proc: subprocess.Popen | None = None
_registry_path: Path | None           = None
_registry: dict                       = {}


def _save_registry() -> None:
    if _registry_path is None:
        return
    _registry_path.parent.mkdir(parents=True, exist_ok=True)
    _registry_path.write_text(json.dumps(_registry, indent=2))


def _sigint_handler(sig, frame):
    logging.warning("SIGINT ricevuto — terminazione processi e salvataggio registry...")
    if _active_proc is not None:
        try:
            _active_proc.terminate()
        except Exception:
            pass
    _save_registry()
    sys.exit(130)


signal.signal(signal.SIGINT,  _sigint_handler)
signal.signal(signal.SIGTERM, _sigint_handler)

# ─── Logging ─────────────────────────────────────────────────────────────────

def setup_logging(sweep_dir: Path, timestamp: str) -> logging.Logger:
    sweep_dir.mkdir(parents=True, exist_ok=True)
    log_path = sweep_dir / f"optimizer_{timestamp}.log"
    fmt = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(level=logging.INFO, format=fmt,
                        handlers=[logging.FileHandler(log_path),
                                  logging.StreamHandler(sys.stdout)])
    return logging.getLogger("optimizer")

# ─── Binari & config ─────────────────────────────────────────────────────────

def find_binary(name: str) -> Path:
    p = BINARIES[name]
    if p.exists():
        return p
    # Fallback: single-config build (Unix Makefiles, nessun subdir Release/)
    alt = Path(str(p).replace("/Release/", "/"))
    if alt.exists():
        return alt
    raise FileNotFoundError(
        f"Binario '{name}' non trovato: {p}\n"
        f"  Alternativa cercata: {alt}\n"
        f"  Esegui il build su questa macchina."
    )


def load_base_config() -> dict:
    path = PROJECT_ROOT / "config" / "config.json"
    return json.loads(path.read_text())


def load_analysis_params(override_path: str | None) -> dict:
    path = Path(override_path) if override_path else PROJECT_ROOT / "config" / "analysis_params.json"
    return json.loads(path.read_text())


def generate_sweep_config(base: dict, geom: dict, margin: float,
                          fast: bool, tag: str, mobile_focus: bool = False) -> Path:
    cfg = dict(base)
    cfg.update({
        "source_dx":       geom["source_dx"],
        "source_dy":       geom["source_dy"],
        "source_x_min":    geom["source_x_min"],
        "source_x_max":    geom["source_x_max"],
        "source_y_max":    geom["source_y_max"],
        "lens_gap_margin": margin,
    })
    if fast:
        cfg["x_max"] = 150.0
    elif not mobile_focus and "x_det" in cfg:
        limit = cfg["x_det"] - cfg.get("lens_det_gap", 0.0)
        if cfg.get("x_max", limit) > limit:
            cfg["x_max"] = limit
    path = PROJECT_ROOT / "config" / f"config_sweep_{tag}.json"
    path.write_text(json.dumps(cfg, indent=2))
    return path

# ─── SSD discovery ───────────────────────────────────────────────────────────

def find_ssd_root_from_stdout(stdout: str, prefix: str) -> Path | None:
    """Estrae il path dal messaggio '[DONE]  Output finale: ...' di run.sh."""
    for line in stdout.splitlines():
        if "[DONE]" in line and "Output finale:" in line:
            parts = line.split("Output finale:")
            if len(parts) == 2:
                candidate = Path(parts[1].strip())
                if candidate.exists():
                    return candidate
    return None


def find_latest_ssd_root(ssd_mount: Path, prefix: str, after: float = 0.0) -> Path | None:
    runs_dir = ssd_mount / "riptide" / "runs"
    if not runs_dir.exists():
        return None
    candidates = sorted(
        (p for p in runs_dir.glob(f"run_*/{prefix}.root") if p.stat().st_mtime > after),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    return candidates[0] if candidates else None

# ─── Esecuzione comandi ───────────────────────────────────────────────────────

def run_cmd(cmd: list[str], log_path: Path, timeout_sec: int,
            dry_run: bool = False, env: dict | None = None) -> tuple[bool, str, str]:
    """Esegue un comando, tee su log_path. Restituisce (ok, stdout, stderr)."""
    global _active_proc
    cmd_str = " ".join(str(c) for c in cmd)
    logging.info(f"CMD: {cmd_str}")
    if dry_run:
        logging.info("[DRY-RUN] (non eseguito)")
        return True, "", ""

    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "a") as lf:
        lf.write(f"\n{'='*60}\nCMD: {cmd_str}\n{'='*60}\n")
        try:
            proc = subprocess.Popen(
                [str(c) for c in cmd],
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, env=env
            )
            _active_proc = proc
            stdout_lines = []
            for line in proc.stdout:
                lf.write(line)
                stdout_lines.append(line)
            proc.wait(timeout=timeout_sec)
            _active_proc = None
            stdout = "".join(stdout_lines)
            ok = proc.returncode == 0
            if not ok:
                logging.error(f"FALLITO (rc={proc.returncode}): {cmd_str}")
                tail = stdout_lines[-30:] if len(stdout_lines) > 30 else stdout_lines
                if tail:
                    logging.error("--- Output (ultime righe) ---")
                    for line in tail:
                        logging.error("  | %s", line.rstrip())
                    logging.error("--- Fine output ---")
            return ok, stdout, ""
        except subprocess.TimeoutExpired:
            proc.kill()
            _active_proc = None
            logging.error(f"TIMEOUT ({timeout_sec}s): {cmd_str}")
            return False, "", "timeout"
        except Exception as e:
            _active_proc = None
            logging.error(f"ECCEZIONE: {e}")
            return False, "", str(e)

# ─── Pre-flight checks ───────────────────────────────────────────────────────

MIN_FREE_LENS_GB = 1300   # picco atteso durante hadd rolling lens (2×lens_merged + batch×chunk)

def check_disk_space() -> tuple[bool, str, str]:
    out_dir = PROJECT_ROOT / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    stat = shutil.disk_usage(out_dir)
    free_gb = stat.free / (1024 ** 3)
    if free_gb < MIN_FREE_LENS_GB:
        return False, "disk_space", (
            f"Spazio libero insufficiente su {out_dir.resolve()}: "
            f"{free_gb:.1f} GB liberi, richiesti {MIN_FREE_LENS_GB} GB per hadd rolling lens."
        )
    logging.info(f"Pre-flight disk: {free_gb:.1f} GB liberi su {out_dir.resolve()}")
    return True, "", ""


def check_margin_consistency() -> tuple[bool, str, str]:
    key = 'config.value("lens_gap_margin"'
    for rel in ["src/lens_simulation/lens_scan.cpp",
                "src/optimization/optimizer.cpp"]:
        path = PROJECT_ROOT / rel
        if not path.exists():
            continue
        if key not in path.read_text():
            return False, "margin_mismatch", f"{rel}: lens_gap_margin non letto da config"
    return True, "", ""


def check_nhits_branch_type() -> tuple[bool, str, str]:
    bad_key = 'CreateNtupleIColumn("n_hits"'
    for rel in ["src/lens_simulation/lens_scan.cpp",
                "src/optimization/optimizer.cpp",
                "src/psf_dof_scan/psf_dof_scan.cpp"]:
        path = PROJECT_ROOT / rel
        if not path.exists():
            continue
        if bad_key in path.read_text():
            return False, "nhits_type_mismatch", f"{rel}: n_hits è IColumn (deve essere DColumn)"
    return True, "", ""


def preflight_checks() -> list[str]:
    failures = []
    for checker in [check_margin_consistency, check_nhits_branch_type, check_disk_space]:
        ok, fid, msg = checker()
        if not ok:
            logging.error(f"Pre-flight FAIL [{fid}]: {msg}")
            failures.append(fid)
        else:
            logging.info(f"Pre-flight OK: {checker.__name__}")
    return failures

# ─── Validazione output ───────────────────────────────────────────────────────

def _read_tsv_header(tsv: Path) -> list[str]:
    if not tsv.exists():
        return []
    with open(tsv) as f:
        header = f.readline().strip()
    return [h.strip() for h in header.split("\t")]


def check_ee80_mean_column(tsv: Path) -> tuple[bool, str, str]:
    if not tsv.exists():
        return False, "ee80_missing", f"{tsv}: file non trovato"
    cols = _read_tsv_header(tsv)
    if "EE80_mean" not in cols:
        return False, "ee80_missing", f"{tsv}: colonna EE80_mean assente (trovate: {cols})"
    return True, "", ""


def check_delta_m_column(tsv: Path) -> tuple[bool, str, str]:
    if not tsv.exists():
        return False, "delta_m_missing", f"{tsv}: file non trovato"
    cols = _read_tsv_header(tsv)
    if "M_abs_err" not in cols:
        return False, "delta_m_missing", f"{tsv}: colonna M_abs_err assente (trovate: {cols})"
    return True, "", ""


def check_mtot_normalization(tsv: Path) -> tuple[bool, str, str]:
    if not tsv.exists():
        return False, "mtot_unnormalized", f"{tsv}: file non trovato"
    with open(tsv) as f:
        lines = f.readlines()
    if len(lines) < 2:
        return True, "", ""   # nessuna riga dati, non validabile
    header = [h.strip() for h in lines[0].split("\t")]
    if "Mtot" not in header:
        return False, "mtot_unnormalized", f"{tsv}: colonna Mtot assente"
    idx = header.index("Mtot")
    for line in lines[1:]:
        parts = line.split("\t")
        if idx >= len(parts):
            continue
        try:
            val = float(parts[idx].strip())
        except ValueError:
            continue
        if not (0.0 <= val <= 1.0 + 1e-6):
            return False, "mtot_unnormalized", f"{tsv}: Mtot={val} fuori [0,1]"
    return True, "", ""


def validate_outputs(sweep_dir: Path) -> list[str]:
    failures = []
    checks = [
        (check_ee80_mean_column,    sweep_dir / "resolution_map.tsv"),
        (check_delta_m_column,      sweep_dir / "dof_map.tsv"),
        (check_mtot_normalization,  sweep_dir / "pareto_results.tsv"),
    ]
    for checker, tsv in checks:
        ok, fid, msg = checker(tsv)
        if not ok:
            logging.error(f"Validazione FAIL [{fid}]: {msg}")
            failures.append(fid)
        else:
            logging.info(f"Validazione OK: {tsv.name}")
    return failures

# ─── Fix dispatch ─────────────────────────────────────────────────────────────

FIX_TARGETS = {
    "ee80_missing":        ["resolution_map"],
    "delta_m_missing":     ["dof_map"],
    "mtot_unnormalized":   ["pareto_selector"],
    "margin_mismatch":     ["lens_simulation_main", "optimization_main"],
    "nhits_type_mismatch": ["lens_simulation_main", "optimization_main", "psf_dof_scan_main"],
}


_MANUAL_FIX_MSG = {
    "margin_mismatch": (
        "Correzione manuale richiesta: lens_gap_margin deve essere letto da config "
        "(config.value(\"lens_gap_margin\", ...)) in:\n"
        "  src/lens_simulation/lens_scan.cpp\n"
        "  src/optimization/optimizer.cpp"
    ),
    "nhits_type_mismatch": (
        "Correzione manuale richiesta: sostituire CreateNtupleIColumn(\"n_hits\" con "
        "CreateNtupleDColumn(\"n_hits\" in:\n"
        "  src/lens_simulation/lens_scan.cpp\n"
        "  src/optimization/optimizer.cpp\n"
        "  src/psf_dof_scan/psf_dof_scan.cpp"
    ),
}


def apply_fix(failure_id: str, dry_run: bool = False) -> bool:
    if failure_id in _MANUAL_FIX_MSG:
        logging.error(f"[{failure_id}] Fix automatico non disponibile.\n{_MANUAL_FIX_MSG[failure_id]}")
        return False

    targets = FIX_TARGETS.get(failure_id, [])
    for target in targets:
        cmd = ["cmake", "--build", str(BUILD), "--config", "Release", "--target", target]
        logging.info(f"Ricompilo target: {target}")
        if dry_run:
            logging.info(f"[DRY-RUN] {' '.join(cmd)}")
            continue
        r = subprocess.run(cmd, cwd=PROJECT_ROOT, capture_output=True, text=True)
        if r.returncode != 0:
            logging.error(f"Build fallita per {target}:\n{r.stdout}\n{r.stderr}")
            return False
        logging.info(f"Build OK: {target}")
    return True

# ─── Pipeline per variante ───────────────────────────────────────────────────

def run_simulation_step(mode: str, sweep_cfg_path: Path, tag: str,
                        sweep_dir: Path, ssd_mount: Path, jobs: int,
                        l1_id: str, l2_id: str,
                        dry_run: bool,
                        all_lenses: bool = False,
                        lens_subset: str = "",
                        use_local: bool = False,
                        focus_tsv: Path | None = None) -> Path | None:
    """Esegue run.sh per la modalità indicata e restituisce il path del file ROOT prodotto."""
    out_names = {"opt": "events", "lens": "lens", "dof": "focal", "psf-dof": "psf_dof"}
    prefix = out_names[mode]

    reg_key = tag
    if reg_key in _registry and prefix in _registry[reg_key]:
        candidate = Path(_registry[reg_key][prefix])
        if candidate.exists():
            logging.info(f"Registry hit [{mode}]: {candidate}")
            return candidate
        else:
            logging.warning(f"Registry entry obsoleta per {mode} ({candidate}), rieseguo.")

    run_target = "local" if use_local else "ssd"
    cmd = [
        find_binary("run_sh"), mode, run_target,
        "--jobs", str(jobs),
        "--config", str(sweep_cfg_path),
    ]
    if not use_local:
        cmd += ["--ssd-mount", str(ssd_mount)]
    if l1_id:
        cmd += ["--l1-id", l1_id]
    if l2_id:
        cmd += ["--l2-id", l2_id]
    if all_lenses:
        cmd.append("--all-lenses")
    if lens_subset:
        cmd += ["--lens-subset", lens_subset]
    if focus_tsv and mode in ("opt", "lens"):
        cmd += ["--focus-tsv", str(focus_tsv)]

    log_path = sweep_dir / "pipeline.log"
    step_start = time.time()
    ok, stdout, _ = run_cmd(cmd, log_path, TIMEOUTS_SEC[mode], dry_run)
    if not ok:
        return None

    # Trova il path dall'output di run.sh
    if dry_run:
        return Path(f"/dev/null/{prefix}.root")   # sentinel per dry-run

    root_path = find_ssd_root_from_stdout(stdout, prefix)
    if root_path is None and not use_local:
        root_path = find_latest_ssd_root(ssd_mount, prefix, after=step_start)
    if root_path is None:
        logging.error(f"File ROOT '{prefix}.root' non trovato dopo {mode}")
        return None

    # Aggiorna registry
    if reg_key not in _registry:
        _registry[reg_key] = {}
    _registry[reg_key][prefix] = str(root_path)
    _save_registry()

    # Copia locale (riferimento, non duplicato pesante)
    local_link = sweep_dir / f"{prefix}.root"
    if not dry_run and not local_link.exists():
        try:
            local_link.symlink_to(root_path)
        except Exception:
            pass  # se symlink fallisce, ok — usiamo il path SSD direttamente

    return root_path


def run_pipeline(tag: str, sweep_cfg: Path, sweep_dir: Path,
                 ssd_mount: Path, jobs: int, l1_id: str, l2_id: str,
                 ap: dict, dry_run: bool,
                 prebuilt_events: Path | None = None,
                 use_local: bool = False,
                 focus_tsv: Path | None = None,
                 prebuilt_dof_tsv: Path | None = None,
                 keep_lens: bool = False,
                 keep_psf_dof: bool = False) -> dict | None:
    """
    Esegue i 10 step della pipeline per una variante.
    Restituisce dict con i path dei TSV/ROOT prodotti, o None se step critico fallisce.
    """
    out = {}
    log = sweep_dir / "pipeline.log"
    sweep_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: opt
    if prebuilt_events is not None:
        events_root = prebuilt_events
        logging.info(f"[opt] Uso events pre-costruito: {events_root}")
    else:
        events_root = run_simulation_step("opt", sweep_cfg, tag, sweep_dir,
                                          ssd_mount, jobs, l1_id, l2_id, dry_run,
                                          use_local=use_local)
        if events_root is None:
            return None
    out["events"] = events_root

    # Salva mappa efficienza (non-ranking)
    plots_dir = sweep_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    eff_png = plots_dir / "efficiency2D.png"
    cmd_p2d = [find_binary("plot2D"),
               "--input",  str(events_root),
               "--config", str(sweep_cfg),
               "--output", str(eff_png)]
    ok_p2d, _, _ = run_cmd(cmd_p2d, log, TIMEOUTS_SEC["analysis"], dry_run)
    if ok_p2d:
        out["efficiency_png"] = eff_png
    else:
        logging.warning("plot2D per mappa efficienza fallito — continuo")

    # Step 2: lens
    lens_root = run_simulation_step("lens", sweep_cfg, tag, sweep_dir,
                                    ssd_mount, jobs, l1_id, l2_id, dry_run,
                                    use_local=use_local, focus_tsv=focus_tsv)
    if lens_root is None:
        return None
    out["lens"] = lens_root

    # Step 3: psf_extractor
    psf_data = sweep_dir / "psf_data.root"
    min_hits = ap["psf_extractor"]["min_hits"]
    cmd = [find_binary("psf_extractor"), str(lens_root), str(psf_data), str(min_hits)]
    ok, _, _ = run_cmd(cmd, log, TIMEOUTS_SEC["analysis"], dry_run)
    if not ok:
        return None
    out["psf_data"] = psf_data

    if not keep_lens and not dry_run and lens_root and lens_root.exists():
        lens_size_gb = lens_root.stat().st_size / (1024 ** 3)
        logging.info(f"Eliminazione lens.root ({lens_size_gb:.1f} GB) — psf_data.root prodotto.")
        try:
            lens_root.unlink()
        except Exception as e:
            logging.warning(f"Impossibile eliminare {lens_root}: {e}")

    # Step 4: q_map
    q_tsv = sweep_dir / "q_map.tsv"
    qp = ap["q_map"]
    cmd = [
        find_binary("q_map"),
        "--psf", str(psf_data),
        "--config", str(sweep_cfg),
        "--n-tracks", str(qp["n_tracks"]),
        "--dt", str(qp["dt"]),
        "--min-hits", str(qp["min_hits"]),
        "--trace-frac", str(qp["trace_frac"]),
        "--tsv", str(q_tsv),
        "--jobs", str(jobs),
        "--output", str(sweep_dir / "plots" / "q_map.png"),
    ]
    if qp.get("dist_to_target"):
        cmd.append("--dist-to-target")
    ok, _, _ = run_cmd(cmd, log, TIMEOUTS_SEC["analysis"], dry_run)
    if not ok:
        return None
    out["q_tsv"] = q_tsv

    # Step 5: chi2_map
    chi2_tsv = sweep_dir / "chi2_map.tsv"
    cp = ap["chi2_map"]
    cmd = [
        find_binary("chi2_map"),
        "--psf", str(psf_data),
        "--config", str(sweep_cfg),
        "--min-hits", str(cp["min_hits"]),
        "--p-low", str(cp["p_low"]),
        "--p-high", str(cp["p_high"]),
        "--tsv", str(chi2_tsv),
        "--jobs", str(jobs),
        "--output", str(sweep_dir / "plots" / "chi2_map.png"),
    ]
    if cp.get("adaptive_target"):
        cmd.append("--adaptive-target")
    ok, _, _ = run_cmd(cmd, log, TIMEOUTS_SEC["analysis"], dry_run)
    if not ok:
        return None
    out["chi2_tsv"] = chi2_tsv

    # Step 6: dof  /  Step 7: dof_map
    if prebuilt_dof_tsv is not None:
        dof_tsv = prebuilt_dof_tsv
        logging.info(f"[dof] Uso dof_tsv pre-costruito: {dof_tsv}")
        local_dof = sweep_dir / "dof_map.tsv"
        if not local_dof.exists():
            local_dof.symlink_to(prebuilt_dof_tsv)
    else:
        focal_root = run_simulation_step("dof", sweep_cfg, tag, sweep_dir,
                                         ssd_mount, jobs, l1_id, l2_id, dry_run,
                                         use_local=use_local)
        if focal_root is None:
            return None
        out["focal"] = focal_root

        dof_tsv = sweep_dir / "dof_map.tsv"
        dp = ap["dof_map"]
        cmd = [
            find_binary("dof_map"),
            "--input", str(focal_root),
            "--config", str(sweep_cfg),
            "--core-fraction", str(dp["core_fraction"]),
            "--m-target", str(dp["m_target"]),
            "--tsv", str(dof_tsv),
            "--jobs", str(jobs),
            "--output", str(sweep_dir / "plots"),
        ]
        ok, _, _ = run_cmd(cmd, log, TIMEOUTS_SEC["analysis"], dry_run)
        if not ok:
            return None
    out["dof_tsv"] = dof_tsv

    # Step 8: psf-dof
    psf_dof_root = run_simulation_step("psf-dof", sweep_cfg, tag, sweep_dir,
                                       ssd_mount, jobs, l1_id, l2_id, dry_run,
                                       use_local=use_local)
    if psf_dof_root is None:
        return None
    out["psf_dof"] = psf_dof_root

    # Step 9: psf_dof_extractor → resolution_map.tsv + psf_dof_data.root
    psf_dof_data = sweep_dir / "psf_dof_data.root"
    res_tsv      = sweep_dir / "resolution_map.tsv"
    cmd = [
        find_binary("psf_dof_extractor"),
        str(psf_dof_root),
        str(psf_dof_data),
        str(res_tsv),
        str(sweep_cfg),
    ]
    ok, _, _ = run_cmd(cmd, log, TIMEOUTS_SEC["psf_dof_extractor"], dry_run)
    if not ok:
        return None
    out["psf_dof_data"] = psf_dof_data
    out["res_tsv"]      = res_tsv

    if not keep_psf_dof and not dry_run and psf_dof_root and psf_dof_root.exists():
        psf_dof_size_gb = psf_dof_root.stat().st_size / (1024 ** 3)
        logging.info(f"Eliminazione psf_dof.root ({psf_dof_size_gb:.1f} GB) — psf_dof_data.root prodotto.")
        try:
            psf_dof_root.unlink()
        except Exception as e:
            logging.warning(f"Impossibile eliminare {psf_dof_root}: {e}")

    # Step 10: pareto_selector
    pareto_tsv = sweep_dir / "pareto_results.tsv"
    pp = ap["pareto_selector"]
    cmd = [
        find_binary("pareto_selector"),
        "--events", str(events_root),
        "--qmap", str(q_tsv),
        "--chi2map", str(chi2_tsv),
        "--dofmap", str(dof_tsv),
        "--resolution", str(res_tsv),
        "--tsv", str(pareto_tsv),
        "--ee80-max", str(pp["ee80_max"]),
        "--w-eta", str(pp["w_eta"]),
        "--w-Q", str(pp["w_Q"]),
        "--w-dof", str(pp["w_dof"]),
        "--w-M", str(pp["w_M"]),
        "--output", str(sweep_dir / "plots" / "pareto_plot.png"),
    ]
    if l1_id:
        cmd += ["--l1-id", l1_id]
    if l2_id:
        cmd += ["--l2-id", l2_id]
    ok, _, _ = run_cmd(cmd, log, TIMEOUTS_SEC["analysis"], dry_run)
    if not ok:
        return None
    out["pareto_tsv"] = pareto_tsv

    return out

# ─── Ranking Pareto ───────────────────────────────────────────────────────────

def _mean_mtot(pareto_tsv: Path) -> float:
    if not pareto_tsv.exists():
        return -1.0
    with open(pareto_tsv) as f:
        lines = f.readlines()
    if len(lines) < 2:
        return -1.0
    header = [h.strip() for h in lines[0].split("\t")]
    if "Mtot" not in header:
        return -1.0
    idx = header.index("Mtot")
    vals = []
    for line in lines[1:]:
        parts = line.split("\t")
        if idx < len(parts):
            try:
                vals.append(float(parts[idx].strip()))
            except ValueError:
                pass
    return sum(vals) / len(vals) if vals else -1.0


def top_k_variants(results: dict, k: int) -> list[str]:
    """Ordina i tag per Mtot medio decrescente e restituisce i top-k."""
    ranked = sorted(results.keys(),
                    key=lambda t: _mean_mtot(results[t].get("pareto_tsv", Path("/dev/null"))),
                    reverse=True)
    return ranked[:k]

# ─── Report Markdown ──────────────────────────────────────────────────────────

def _load_pareto_rows(tsv: Path) -> list[dict]:
    if not tsv or not tsv.exists():
        return []
    with open(tsv) as f:
        lines = f.readlines()
    if len(lines) < 2:
        return []
    header = [h.strip() for h in lines[0].split("\t")]
    rows = []
    for line in lines[1:]:
        parts = [p.strip() for p in line.split("\t")]
        if len(parts) == len(header):
            rows.append(dict(zip(header, parts)))
    return rows


def generate_report(fast_results: dict, full_results: dict,
                    sweep_dir: Path, timestamp: str,
                    l1_id: str, l2_id: str,
                    t_start: float) -> Path:
    elapsed = time.time() - t_start
    h, rem = divmod(int(elapsed), 3600)
    m = rem // 60
    duration = f"{h}h {m}m"

    n_done = len(fast_results) + len(full_results)
    report_path = sweep_dir / f"report_{timestamp}.md"

    # Raccoglie tutte le righe Pareto da full (poi fast) con tag variante
    all_rows: list[tuple[str, dict]] = []
    for tag, out in full_results.items():
        tsv = out.get("pareto_tsv")
        for row in _load_pareto_rows(tsv):
            all_rows.append((tag, row))
    for tag, out in fast_results.items():
        if tag in full_results:
            continue
        tsv = out.get("pareto_tsv")
        for row in _load_pareto_rows(tsv):
            all_rows.append((tag, row))

    # Ordina per Mtot
    def _mtot(pair):
        try:
            return float(pair[1].get("Mtot", -1))
        except ValueError:
            return -1.0

    all_rows.sort(key=_mtot, reverse=True)
    top10 = all_rows[:10]

    lines = [
        "# RIPTIDE Autonomous Optimization Report",
        "",
        f"Generated: {datetime.datetime.now().isoformat(timespec='seconds')}  "
        f"|  Duration: {duration}  |  Variants: {n_done}/9",
        f"Lenti: L1={l1_id or 'N/A'}  L2={l2_id or 'N/A'}",
        "",
        "## Sweep Matrix",
        "",
        "| Variante | source_dx | margin | Tipo | Status |",
        "|----------|-----------|--------|------|--------|",
    ]
    for tag, out in {**fast_results, **full_results}.items():
        phase = "full" if tag in full_results else "fast"
        parts = tag.split("_")
        # tag formato: {geom_name}_m{margin}_{fast|full}
        status = "OK" if out else "FAIL"
        margin_str = tag.split("_m")[1].split("_")[0] if "_m" in tag else "?"
        geom_name  = tag.split("_m")[0]       if "_m" in tag else tag
        try:
            geom = next(g for g in SWEEP_GEOM if g["name"] == geom_name)
            dx = geom["source_dx"]
        except StopIteration:
            dx = "?"
        lines.append(f"| {tag} | {dx} | {margin_str} | {phase} | {status} |")

    lines += [
        "",
        "## Top-10 Pareto Configurations (ranked by Mtot)",
        "",
        "| Rank | x1 | x2 | η | Q | M_abs_err | EE80 | DoF | Mtot | Variante |",
        "|------|----|----|---|---|-----------|-----------|-----|------|----------|",
    ]
    for rank, (tag, row) in enumerate(top10, 1):
        x1   = row.get("x1",       "—")
        x2   = row.get("x2",       "—")
        eta  = row.get("eta",      "—")
        Q    = row.get("Q",        "—")
        merr = row.get("M_abs_err","—")
        ee80 = row.get("EE80", "—")
        dof  = row.get("DoF",  "—")
        mtot = row.get("Mtot",     "—")
        lines.append(f"| {rank} | {x1} | {x2} | {eta} | {Q} | {merr} | {ee80} | {dof} | {mtot} | {tag} |")

    lines += [
        "",
        "## Interpretazione fisica",
        "",
        "M ≈ (x_det − x2) / x1 (regime telescopio).  ",
        "η = efficienza di raccolta fotoni.  ",
        "M_abs_err = errore assoluto di magnificazione rispetto al target.  ",
        "EE80 = raggio del cerchio che contiene l'80% dell'energia (mm).  ",
        "DoF = profondità di campo (mm).  ",
        "",
        "Le configurazioni Pareto-ottimali massimizzano Mtot = w_eta·η̃ + w_Q·Q̃ + w_dof·DoF̃ + w_M·M̃.",
        "",
    ]

    report_path.write_text("\n".join(lines))
    logging.info(f"Report scritto: {report_path}")
    return report_path

# ─── Screening lenti ─────────────────────────────────────────────────────────

def _escape_root_path(p: Path) -> str:
    """Escape path for inline use in a ROOT C-macro string literal."""
    return str(p).replace("\\", "\\\\").replace('"', '\\"')


def generate_screening_config(base_cfg: dict, n_photons: int) -> Path:
    import copy
    cfg = copy.deepcopy(base_cfg)
    cfg["n_photons"] = n_photons
    path = PROJECT_ROOT / "config" / "config_screening.json"
    path.write_text(json.dumps(cfg, indent=2))
    return path


def read_efl(item_id: str) -> tuple[float, float]:
    """Restituisce (focal_length_mm, center_thickness_mm) da planoconvex o biconvex."""
    import csv as _csv
    catalogs = [
        PROJECT_ROOT / "lens_cutter" / "lens_data" / "thorlabs_planoconvex.tsv",
        PROJECT_ROOT / "lens_cutter" / "lens_data" / "thorlabs_biconvex.tsv",
    ]
    for tsv in catalogs:
        with open(tsv) as f:
            for row in _csv.DictReader(f, delimiter='\t'):
                if row['Item #'].strip() == item_id.strip():
                    return float(row['Focal Length']), float(row['Center Thickness'])
    raise KeyError(f"Lente {item_id!r} non trovata nei cataloghi")


def compute_focus_tsv(cfg: dict, lens_subset: str, out_path: Path) -> Path:
    """Stima thin-lens del piano focale per ogni (x1, x2) della griglia.

    Per ogni coppia ordinata (lens_a, lens_b) nel subset, calcola x_focus con la
    formula thin-lens. Per ogni punto (x1, x2) prende la mediana dei fuochi validi
    su tutte le coppie, così la stima è rappresentativa dell'intero subset.
    """
    import math as _math, statistics as _stats
    from collections import defaultdict as _dd

    subset_ids = [s.strip() for s in lens_subset.split(",") if s.strip()]
    if len(subset_ids) < 2:
        raise ValueError(f"lens_subset deve contenere ≥2 lenti, ricevuto: {lens_subset!r}")

    # Leggi EFL e spessore per tutte le lenti del subset (entrambi i cataloghi)
    lens_data: dict[str, tuple[float, float]] = {}
    for lid in subset_ids:
        try:
            lens_data[lid] = read_efl(lid)
        except KeyError as e:
            logging.warning(f"compute_focus_tsv: {e} — lente ignorata")

    if len(lens_data) < 2:
        raise ValueError("Meno di 2 lenti trovate nei cataloghi per il subset indicato")

    found = list(lens_data.keys())
    # Tutte le coppie ordinate (a come lens1, b come lens2 e viceversa)
    pairs = [(a, b) for i, a in enumerate(found) for b in found if a != b]

    x_min        = cfg['x_min']
    x_max        = cfg['x_max']
    dx           = cfg['dx']
    lens_det_gap = cfg.get('lens_det_gap', 10.0)
    margin_col   = 1.0  # margine collisione identico a run.sh

    def arange_grid(start, stop, step):
        vals, v = [], start
        while v < stop - 1e-9:
            vals.append(round(v, 6))
            v = round(v + step, 9)
        return vals

    x1_vals = arange_grid(x_min, x_max, dx)
    x2_vals = arange_grid(x_min, x_max, dx)

    focus_map: dict[tuple[float, float], list[float]] = _dd(list)

    for lid1, lid2 in pairs:
        efl1, h1 = lens_data[lid1]
        efl2, h2 = lens_data[lid2]
        for x1 in x1_vals:
            for x2 in x2_vals:
                if x2 <= x1 + (h1 + h2) / 2.0 + margin_col:
                    continue
                try:
                    if x1 < 1e-6 or abs(x1 - efl1) < 1e-6:
                        continue
                    v1 = 1.0 / (1.0/efl1 - 1.0/x1)
                    u2 = v1 - (x2 - x1)
                    if abs(u2) < 1e-6 or abs(u2 - efl2) < 1e-6:
                        continue
                    v2 = 1.0 / (1.0/efl2 - 1.0/u2)
                    x_focus = x2 + v2
                    if x_focus < x2 + lens_det_gap:
                        continue
                    focus_map[(x1, x2)].append(x_focus)
                except (ZeroDivisionError, OverflowError):
                    continue

    rows = []
    for (x1, x2), focuses in sorted(focus_map.items()):
        rows.append((x1, x2, _stats.median(focuses)))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write("x1\tx2\tx_focus\n")
        for r in rows:
            f.write(f"{r[0]:.3f}\t{r[1]:.3f}\t{r[2]:.3f}\n")
    logging.info(
        f"focus TSV (thin-lens, {len(pairs)} coppie, {len(lens_data)} lenti): "
        f"{len(rows)} punti griglia → {out_path}"
    )
    return out_path


def extract_focus_tsv(dof_tsv: Path, out_path: Path, dry_run: bool = False) -> Path:
    """Estrae (x1, x2, x_focus) da dof_map.tsv, filtrando focus_before_lens2=1."""
    if dry_run:
        logging.info(f"[DRY-RUN] extract_focus_tsv {dof_tsv} → {out_path}")
        return out_path
    import csv as _csv
    rows = []
    with open(dof_tsv) as f:
        for row in _csv.DictReader(f, delimiter='\t'):
            if row.get('focus_before_lens2', '0').strip() == '1':
                continue
            rows.append((float(row['x1']), float(row['x2']), float(row['x_focus'])))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        f.write("x1\tx2\tx_focus\n")
        for r in rows:
            f.write(f"{r[0]:.3f}\t{r[1]:.3f}\t{r[2]:.3f}\n")
    logging.info(f"focus_accurate.tsv: {len(rows)} righe → {out_path}")
    return out_path


def run_plot2d(events_root: Path, cfg_path: Path, out_dir: Path,
               dry_run: bool, ranking_only: bool = False) -> Path | None:
    out_png  = out_dir / "screening_efficiency2D.png"
    out_rank = out_dir / "screening_efficiency2D_ranking.csv"
    cmd = [find_binary("plot2D"),
           "--input",  str(events_root),
           "--config", str(cfg_path),
           "--output", str(out_png)]
    if ranking_only:
        cmd.append("--ranking-only")
    log_path = out_dir / "plot2d_screening.log"
    ok, _, _ = run_cmd(cmd, log_path, TIMEOUTS_SEC["analysis"], dry_run)
    if not ok:
        return None
    return out_rank if (dry_run or out_rank.exists()) else None


def parse_ranking_csv(ranking_csv: Path, top_n: int) -> list[tuple[str, str]]:
    import csv as csv_mod
    pairs: list[tuple[str, str]] = []
    with open(ranking_csv, newline="") as f:
        for row in csv_mod.DictReader(f):
            pairs.append((row["lens1_id"], row["lens2_id"]))
            if len(pairs) == top_n:
                break
    return pairs


def filter_events_root(events_root: Path, selected_pairs: list[tuple[str, str]],
                       out_path: Path, dry_run: bool) -> bool:
    if dry_run:
        logging.info(f"[DRY-RUN] filter {events_root} → {out_path}")
        return True

    conds = " || ".join(
        f'(strcmp(l1_id,"{l75}")==0 && strcmp(l2_id,"{l60}")==0)'
        for l75, l60 in selected_pairs
    )

    events_root_s = _escape_root_path(events_root)
    out_path_s    = _escape_root_path(out_path)
    macro = f"""
void do_filter() {{
  TFile *fin  = TFile::Open("{events_root_s}");
  TTree *cfg  = (TTree*)fin->Get("Configurations");
  TTree *eff  = (TTree*)fin->Get("Efficiency");

  TFile *fout = new TFile("{out_path_s}", "RECREATE");
  fout->SetCompressionAlgorithm(ROOT::kLZ4);
  fout->SetCompressionLevel(4);

  TTree *cfg2 = cfg->CopyTree("{conds}");
  std::set<int> keep_ids;
  int cid;
  cfg2->SetBranchAddress("config_id", &cid);
  for (long long i = 0; i < cfg2->GetEntries(); ++i) {{
    cfg2->GetEntry(i);
    keep_ids.insert(cid);
  }}
  cfg2->Write("Configurations");

  int eid;
  eff->SetBranchAddress("config_id", &eid);
  TTree *eff2 = eff->CloneTree(0);
  for (long long i = 0; i < eff->GetEntries(); ++i) {{
    eff->GetEntry(i);
    if (keep_ids.count(eid)) eff2->Fill();
  }}
  eff2->Write("Efficiency");
  fout->Close();
  fin->Close();
}}
"""
    macro_path = out_path.parent / "_filter_lens.C"
    macro_path.write_text(macro)
    r = subprocess.run(
        ["root", "-b", "-q", str(macro_path)],
        cwd=PROJECT_ROOT, capture_output=True, text=True, timeout=120
    )
    macro_path.unlink(missing_ok=True)
    if r.returncode != 0:
        logging.error(f"Filtraggio ROOT fallito:\n{r.stderr}")
        return False
    return True


# ─── Loop principale ──────────────────────────────────────────────────────────

def build_variants(geom_filter: str | None, margin_filter: float | None) -> list[tuple[dict, float]]:
    variants = []
    for geom in SWEEP_GEOM:
        if geom_filter and geom["name"] != geom_filter:
            continue
        for margin in MARGIN_VALUES:
            if margin_filter is not None and abs(margin - margin_filter) > 1e-9:
                continue
            variants.append((geom, margin))
    return variants


def make_tag(geom: dict, margin: float, phase: str) -> str:
    return f"{geom['name']}_m{margin}_{phase}"


def make_tag_lens(geom: dict, margin: float, l1: str, l2: str, phase: str) -> str:
    return f"{geom['name']}_m{margin}_{l1}_{l2}_{phase}"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Autonomous optimization driver per RIPTIDE",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Esempi:
              # Dry-run singola variante
              python3 scripts/autonomous_optimizer.py --fast --geom nominal --margin 1.0 \\
                  --l1-id LA4464 --l2-id LA4464R --dry-run

              # Run completo
              python3 scripts/autonomous_optimizer.py --l1-id LA4464 --l2-id LA4464R
        """),
    )
    parser.add_argument("--fast",         action="store_true",
                        help="Solo fase 1 (x_max=150), senza fase 2")
    parser.add_argument("--no-two-phase", dest="two_phase", action="store_false", default=True)
    parser.add_argument("--top-k",        type=int,   default=3,
                        help="Varianti da promuovere in fase 2 (default: 3)")
    parser.add_argument("--geom",         type=str,   default=None,
                        help="Esegui solo la geometria specificata (nominal|coarse|extended)")
    parser.add_argument("--margin",       type=float, default=None,
                        help="Esegui solo il margin specificato (0.5|1.0|2.0)")
    parser.add_argument("--l1-id",    type=str,   default="",
                        help="ID lente L1 Thorlabs (es. LA4464)")
    parser.add_argument("--l2-id",    type=str,   default="",
                        help="ID lente L2 Thorlabs (es. LA4464R)")
    parser.add_argument("--local",         action="store_true",
                        help="Usa modalità local di run.sh (output in output/sweep/, nessuna SSD)")
    parser.add_argument("--ssd-mount",    type=str,   default="/mnt/external_ssd",
                        help="Override SSD mount (default: /mnt/external_ssd)")
    parser.add_argument("--analysis-params", type=str, default=None,
                        help="Override config/analysis_params.json")
    parser.add_argument("--max-hours",    type=float, default=48.0,
                        help="Budget tempo totale in ore (default: 48)")
    parser.add_argument("--jobs",         type=int,   default=os.cpu_count() or 1,
                        help="Numero di job paralleli per run.sh")
    parser.add_argument("--dry-run",      action="store_true",
                        help="Stampa comandi senza eseguire")
    parser.add_argument("--no-commit",    action="store_true",
                        help="Disabilita git commit automatici")
    parser.add_argument("--all-lenses",        action="store_true",
                        help="Fase 0: screening su tutte le coppie di lenti")
    parser.add_argument("--lens-subset",       type=str, default="",
                        help="IDs lenti comma-separated per lo screening (es. LB4592,LB4553)")
    parser.add_argument("--top-n-lenses",      type=int, default=3,
                        help="Top-N coppie da selezionare dallo screening (default: 3)")
    parser.add_argument("--screening-photons", type=int, default=1000,
                        help="n_photons per la fase 0 di screening (default: 1000)")
    parser.add_argument("--mobile-focus",      action="store_true",
                        help="Fuoco mobile: thin-lens screening + DOF accurato per posizione detector")
    parser.add_argument("--refine-photons",    type=int, default=0,
                        help="n_photons per il re-run opt con fuoco accurato (0=usa config.json)")
    parser.add_argument("--dx",               type=float, default=None,
                        help="Passo griglia x1/x2 in mm (default: valore da config.json)")
    parser.add_argument("--keep-lens",       action="store_true", default=False,
                        help="Non eliminare lens.root dopo psf_extractor (debug)")
    parser.add_argument("--keep-psf-dof",   action="store_true", default=False,
                        help="Non eliminare psf_dof.root dopo psf_dof_extractor (debug)")
    args = parser.parse_args()

    timestamp  = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    sweep_dir  = PROJECT_ROOT / "output" / "sweep"
    t_start    = time.time()
    max_sec    = args.max_hours * 3600
    use_local  = args.local
    ssd_mount  = Path(args.ssd_mount)

    global _registry_path, _registry
    _registry_path = sweep_dir / "registry.json"
    if _registry_path.exists():
        _registry = json.loads(_registry_path.read_text())

    setup_logging(sweep_dir, timestamp)
    logging.info(f"=== RIPTIDE Autonomous Optimizer — {timestamp} ===")
    logging.info(f"dry_run={args.dry_run}  two_phase={args.two_phase and not args.fast}  "
                 f"jobs={args.jobs}  max_hours={args.max_hours}")

    base_cfg = load_base_config()
    if args.dx is not None:
        base_cfg["dx"] = args.dx
    ap       = load_analysis_params(args.analysis_params)

    # Pre-flight
    pf_failures = preflight_checks()
    for fid in pf_failures:
        if not apply_fix(fid, args.dry_run):
            logging.error(f"Fix pre-flight [{fid}] fallito — interrotto.")
            return 1

    variants = build_variants(args.geom, args.margin)
    if not variants:
        logging.error("Nessuna variante selezionata con i filtri forniti.")
        return 1
    logging.info(f"Varianti selezionate: {len(variants)}")

    fast_results: dict[str, dict | None] = {}
    full_results: dict[str, dict | None] = {}

    # ── FASE 0: screening lenti (opzionale) ──────────────────────────────────
    prebuilt_events: Path | None = None
    prebuilt_dof_tsv: Path | None = None
    focus_accurate_tsv: Path | None = None
    lens_pairs: list[tuple[str, str]]

    if args.all_lenses:
        screen_dir = sweep_dir / "lens_screening"
        screen_dir.mkdir(parents=True, exist_ok=True)
        screen_cfg = generate_screening_config(base_cfg, args.screening_photons)

        # Thin-lens focus TSV iniziale (solo modalità mobile)
        thin_focus_tsv: Path | None = None
        if args.mobile_focus and args.lens_subset:
            try:
                thin_focus_tsv = compute_focus_tsv(
                    base_cfg, args.lens_subset,
                    screen_dir / "thin_focus.tsv")
            except (KeyError, ValueError) as e:
                logging.warning(f"compute_focus_tsv: {e} — screening senza focus-tsv")

        events_all = run_simulation_step(
            "opt", screen_cfg, "screening", screen_dir, ssd_mount,
            args.jobs, "", "", args.dry_run,
            all_lenses=True, lens_subset=args.lens_subset, use_local=use_local,
            focus_tsv=thin_focus_tsv)
        if events_all is None:
            logging.error("Screening lenti fallito.")
            return 1

        # Ranking-only: solo CSV, nessuna immagine
        ranking_csv = run_plot2d(events_all, screen_cfg, screen_dir, args.dry_run,
                                 ranking_only=True)
        if ranking_csv is None:
            logging.error("plot2D fallito nel calcolo del ranking.")
            return 1

        if args.dry_run:
            top_pairs = [("DRY_L75", "DRY_L60")] * args.top_n_lenses
            logging.info(f"[DRY-RUN] top_pairs sentinel: {top_pairs}")
        else:
            top_pairs = parse_ranking_csv(ranking_csv, args.top_n_lenses)
        logging.info(f"Top-{args.top_n_lenses} coppie: {top_pairs}")

        if args.mobile_focus:
            # DOF per le top coppie → fuochi accurati
            subset_ids = ",".join({l for pair in top_pairs for l in pair})
            dof_screen_dir = screen_dir / "dof_screening"
            dof_screen_dir.mkdir(parents=True, exist_ok=True)
            focal_screen = run_simulation_step(
                "dof", screen_cfg, "screening_dof", dof_screen_dir, ssd_mount,
                args.jobs, "", "", args.dry_run,
                all_lenses=True, lens_subset=subset_ids, use_local=use_local)
            if focal_screen is None:
                logging.error("DOF screening fallito.")
                return 1

            dof_screen_tsv = dof_screen_dir / "dof_map.tsv"
            dp = ap["dof_map"]
            cmd_dm = [
                find_binary("dof_map"),
                "--input", str(focal_screen),
                "--config", str(screen_cfg),
                "--core-fraction", str(dp["core_fraction"]),
                "--m-target", str(dp["m_target"]),
                "--tsv", str(dof_screen_tsv),
            ]
            ok_dm, _, _ = run_cmd(cmd_dm, dof_screen_dir / "dof_map.log",
                                  TIMEOUTS_SEC["analysis"], args.dry_run)
            if not ok_dm:
                logging.error("dof_map screening fallito.")
                return 1

            focus_accurate_tsv = extract_focus_tsv(
                dof_screen_tsv, screen_dir / "focus_accurate.tsv", args.dry_run)
            prebuilt_dof_tsv = dof_screen_tsv

            # Re-run opt con fuoco accurato per efficienza precisa
            n_refine = args.refine_photons if args.refine_photons > 0 else base_cfg.get("n_photons", 10000)
            screen_cfg_refined = generate_screening_config(base_cfg, n_refine)
            events_refined = run_simulation_step(
                "opt", screen_cfg_refined, "screening_refined", screen_dir, ssd_mount,
                args.jobs, "", "", args.dry_run,
                all_lenses=True, lens_subset=subset_ids, use_local=use_local,
                focus_tsv=focus_accurate_tsv)
            if events_refined is None:
                logging.warning("Re-run opt con fuoco accurato fallito — uso events_all")
                events_refined = events_all

            # Plot con immagini (salva in plots/)
            plots_screen_dir = screen_dir / "plots"
            plots_screen_dir.mkdir(parents=True, exist_ok=True)
            run_plot2d(events_refined, screen_cfg_refined, plots_screen_dir,
                       args.dry_run, ranking_only=False)

            prebuilt_events = events_refined
            lens_pairs = top_pairs
        else:
            # Fuoco fisso: filtra events per top coppie
            events_top = screen_dir / "events_top.root"
            if not filter_events_root(events_all, top_pairs, events_top, args.dry_run):
                logging.warning("Filtraggio ROOT fallito — uso file completo")
                events_top = events_all
            prebuilt_events = events_top
            lens_pairs = top_pairs
    else:
        lens_pairs = [(args.l1_id, args.l2_id)]

    # ── FASE 1: fast su tutte le varianti, per ogni coppia ───────────────────
    fast_tag_lens: dict[str, tuple[str, str]] = {}
    for l1, l2 in lens_pairs:
        for geom, margin in variants:
            if time.time() - t_start > max_sec:
                logging.warning(f"Budget tempo ({args.max_hours}h) esaurito — interruzione fase 1.")
                break

            tag = (make_tag_lens(geom, margin, l1, l2, "fast")
                   if args.all_lenses else make_tag(geom, margin, "fast"))
            fast_tag_lens[tag] = (l1, l2)
            sweep_sub = sweep_dir / tag
            logging.info(f"--- Fase 1: {tag} ---")

            sweep_cfg = generate_sweep_config(base_cfg, geom, margin, fast=True, tag=tag,
                                               mobile_focus=args.mobile_focus)
            out = run_pipeline(tag, sweep_cfg, sweep_sub, ssd_mount,
                               args.jobs, l1, l2, ap, args.dry_run,
                               prebuilt_events=prebuilt_events, use_local=use_local,
                               focus_tsv=focus_accurate_tsv,
                               prebuilt_dof_tsv=prebuilt_dof_tsv,
                               keep_lens=args.keep_lens,
                               keep_psf_dof=args.keep_psf_dof)
            fast_results[tag] = out or {}

    # ── FASE 2: full sui top-k (se abilitata) ────────────────────────────────
    if not args.fast and args.two_phase and fast_results:
        top_tags = top_k_variants(fast_results, args.top_k)
        logging.info(f"Top-{args.top_k} da fase 1: {top_tags}")

        for fast_tag in top_tags:
            if time.time() - t_start > max_sec:
                logging.warning("Budget tempo esaurito — interruzione fase 2.")
                break

            # Ricostruisce geom e margin dal tag
            geom_name  = fast_tag.split("_m")[0]
            margin_str = fast_tag.split("_m")[1].split("_")[0]
            try:
                geom   = next(g for g in SWEEP_GEOM if g["name"] == geom_name)
                margin = float(margin_str)
            except (StopIteration, ValueError):
                logging.error(f"Impossibile decodificare tag: {fast_tag}")
                continue

            tag_l1, tag_l2 = fast_tag_lens.get(
                fast_tag, (args.l1_id, args.l2_id))
            full_tag = (make_tag_lens(geom, margin, tag_l1, tag_l2, "full")
                        if args.all_lenses else make_tag(geom, margin, "full"))

            sweep_sub = sweep_dir / full_tag
            logging.info(f"--- Fase 2: {full_tag} ---")

            sweep_cfg = generate_sweep_config(base_cfg, geom, margin, fast=False, tag=full_tag,
                                               mobile_focus=args.mobile_focus)

            retries = 0
            success = False
            while retries < 3 and not success:
                out = run_pipeline(full_tag, sweep_cfg, sweep_sub, ssd_mount,
                                   args.jobs, tag_l1, tag_l2, ap, args.dry_run,
                                   prebuilt_events=prebuilt_events, use_local=use_local,
                                   focus_tsv=focus_accurate_tsv,
                                   prebuilt_dof_tsv=prebuilt_dof_tsv,
                                   keep_lens=args.keep_lens,
                                   keep_psf_dof=args.keep_psf_dof)
                if out is None:
                    retries += 1
                    logging.warning(f"Pipeline fallita per {full_tag} (tentativo {retries}/3)")
                    continue

                failures = validate_outputs(sweep_sub)
                if not failures:
                    success = True
                    full_results[full_tag] = out
                    logging.info(f"Variante {full_tag} completata con successo.")

                    if not args.no_commit and not args.dry_run and (PROJECT_ROOT / ".git").exists():
                        subprocess.run(
                            ["git", "add", str(sweep_sub / "pareto_results.tsv")],
                            cwd=PROJECT_ROOT, check=False
                        )
                        subprocess.run(
                            ["git", "commit", "--allow-empty", "-m",
                             f"sweep: {full_tag} completato"],
                            cwd=PROJECT_ROOT, check=False
                        )
                else:
                    for fid in failures:
                        apply_fix(fid, args.dry_run)
                    retries += 1
                    logging.warning(f"Validazione fallita per {full_tag}, retry {retries}/3")

            if not success:
                logging.error(f"Variante {full_tag} fallita dopo 3 tentativi.")
                full_results[full_tag] = {}

    # ── Report finale ─────────────────────────────────────────────────────────
    if len(lens_pairs) == 1:
        rep_l75, rep_l60 = lens_pairs[0]
    elif lens_pairs:
        rep_l75 = "multi_lens"
        rep_l60 = f"{len(lens_pairs)}_pairs"
    else:
        rep_l75, rep_l60 = args.l1_id, args.l2_id
    report = generate_report(fast_results, full_results, sweep_dir, timestamp,
                             rep_l75, rep_l60, t_start)

    if not args.no_commit and not args.dry_run and (PROJECT_ROOT / ".git").exists():
        subprocess.run(
            ["git", "add", str(report)],
            cwd=PROJECT_ROOT, check=False
        )
        subprocess.run(
            ["git", "commit", "--allow-empty", "-m",
             f"sweep: report finale {timestamp}"],
            cwd=PROJECT_ROOT, check=False
        )

    _save_registry()
    elapsed = time.time() - t_start
    h, rem  = divmod(int(elapsed), 3600)
    m       = rem // 60
    logging.info(f"=== Completato in {h}h {m}m. Report: {report} ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
