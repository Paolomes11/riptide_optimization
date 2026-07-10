#!/usr/bin/env python3
"""
lens_runner.py — Pipeline completa per una singola coppia di lenti RIPTIDE.

Esegue la catena completa (opt → plot2D → lens → psf_extractor → q_map → coverage_map →
chi2_map → cor_map → dof → dof_map → psf-dof → psf_dof_extractor → resolution_map) per una specifica
coppia L1/L2 senza screening iniziale né Pareto selector finale.

File ROOT pesanti → storage esterno (SSD) o locale.
Grafici e TSV        → output/lens_simulations/{l1_id}_{l2_id}/
"""
from __future__ import annotations

import argparse
import csv
import datetime
import json
import logging
import os
import signal
import subprocess
import sys
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
}

TIMEOUTS_SEC = {
    "opt":               60  * 60,
    "lens":              210 * 60,
    "dof":               90  * 60,
    "psf_dof_extractor": 30  * 60,
    "resolution_map":    30  * 60,
    "analysis":          10  * 60,
}

LOG_CSV    = PROJECT_ROOT / "output" / "simulation_log.csv"
LOG_FIELDS = [
    "date", "l1_id", "l2_id", "detector_mode",
    "x_min", "x_max", "dx",
    "n_photons", "lens_n_photons", "dof_n_photons", "psf_dof_n_photons",
    "lens_gap_margin", "storage_type", "storage_path", "output_dir", "notes",
]

_active_proc: subprocess.Popen | None = None


def _sigint_handler(sig, frame):
    logging.warning("SIGINT ricevuto — terminazione processo in corso...")
    if _active_proc is not None:
        try:
            _active_proc.terminate()
        except Exception:
            pass
    sys.exit(130)


signal.signal(signal.SIGINT,  _sigint_handler)
signal.signal(signal.SIGTERM, _sigint_handler)

# ─── Logging ─────────────────────────────────────────────────────────────────

def setup_logging(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "pipeline.log"
    fmt = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(
        level=logging.INFO, format=fmt,
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )

# ─── Binari & config ─────────────────────────────────────────────────────────

def find_binary(name: str) -> Path:
    p = BINARIES[name]
    if p.exists():
        return p
    alt = Path(str(p).replace("/Release/", "/"))
    if alt.exists():
        return alt
    raise FileNotFoundError(
        f"Binario '{name}' non trovato: {p}\n"
        f"  Alternativa cercata: {alt}\n"
        f"  Esegui il build su questa macchina."
    )


def load_base_config(path: str | None = None) -> dict:
    p = Path(path) if path else PROJECT_ROOT / "config" / "config.json"
    return json.loads(p.read_text())


def load_analysis_params(path: str | None = None) -> dict:
    p = Path(path) if path else PROJECT_ROOT / "config" / "analysis_params.json"
    return json.loads(p.read_text())

# ─── SSD discovery ───────────────────────────────────────────────────────────

def find_ssd_root_from_stdout(stdout: str, prefix: str) -> Path | None:
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
            dry_run: bool = False) -> tuple[bool, str, str]:
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
                text=True,
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


def check_ssd_mount(ssd_mount: Path) -> tuple[bool, str, str]:
    if not ssd_mount.exists():
        return False, "ssd_not_found", (
            f"'{ssd_mount}' non esiste o il dispositivo non è montato. "
            "Monta il dispositivo o usa --local."
        )
    if not os.path.ismount(str(ssd_mount)):
        return False, "ssd_not_mountpoint", (
            f"'{ssd_mount}' non è un mountpoint. "
            "Passa il mountpoint (es. /mnt/external_ssd) o usa --local."
        )
    return True, "", ""


def preflight_checks() -> list[str]:
    failures = []
    for checker in [check_margin_consistency, check_nhits_branch_type]:
        ok, fid, msg = checker()
        if not ok:
            logging.error(f"Pre-flight FAIL [{fid}]: {msg}")
            failures.append(fid)
        else:
            logging.info(f"Pre-flight OK: {checker.__name__}")
    return failures

# ─── Fuoco mobile (DOF-based) ────────────────────────────────────────────────

def _remove_focus_outliers(rows: list[tuple[str, str, str]],
                            dx: float = 3.0,
                            thresh: float = 15.0) -> list[tuple[str, str, str]]:
    """Rimuove iterativamente outlier isolati in x_focus all'interno di ogni x1-slice.

    Un punto (x1, x2, xf) è un outlier se differisce di più di `thresh` mm
    dalla mediana dei vicini entro ±2 passi dx sulla stessa riga x1.
    L'approccio iterativo svela outlier "mascherati" da altri outlier vicini
    (tipico con coppie asimmetriche dove lo stray-light crea cascate di
    minimi spuri nel DOF scan).
    """
    import statistics
    from collections import defaultdict

    by_x1: dict[float, list[tuple[float, float]]] = defaultdict(list)
    for x1s, x2s, xfs in rows:
        by_x1[float(x1s)].append((float(x2s), float(xfs)))

    cur: dict[float, list[tuple[float, float]]] = dict(by_x1)
    total_removed = 0

    for _ in range(20):
        next_cur: dict[float, list[tuple[float, float]]] = {}
        removed_this = 0
        for x1, pts in cur.items():
            pts_s = sorted(pts, key=lambda p: p[0])
            x2s_list = [p[0] for p in pts_s]
            xfs_list = [p[1] for p in pts_s]
            kept: list[tuple[float, float]] = []
            for i, (x2, xf) in enumerate(pts_s):
                neighbors = [xfs_list[j] for j in range(len(pts_s))
                             if j != i and abs(x2s_list[j] - x2) <= 2 * dx + 1e-6]
                if len(neighbors) >= 2:
                    med = statistics.median(neighbors)
                    if abs(xf - med) > thresh:
                        logging.warning(
                            f"focus outlier rimosso: x1={x1:.0f} x2={x2:.0f} "
                            f"x_focus={xf:.1f} mm (mediana vicini={med:.1f} mm, "
                            f"delta={abs(xf-med):.1f} mm > soglia={thresh:.0f} mm)"
                        )
                        removed_this += 1
                        total_removed += 1
                        continue
                kept.append((x2, xf))
            next_cur[x1] = kept
        cur = next_cur
        if removed_this == 0:
            break

    if total_removed:
        logging.info(f"focus outlier: rimossi {total_removed} punti su {len(rows)}")

    return [(str(x1), str(x2), str(xf))
            for x1, pts in sorted(cur.items())
            for x2, xf in pts]


def extract_focus_tsv(dof_tsv: Path, out_path: Path, dry_run: bool,
                      dx: float = 3.0) -> Path | None:
    """Estrae (x1, x2, x_focus) da dof_map.tsv, incluse le configurazioni col fuoco
    prima di L2 (già clampate al bordo fisico da dof_map), rimuovendo gli outlier
    isolati causati da minimi spuri nel DOF scan."""
    if dry_run:
        return out_path
    if not dof_tsv.exists():
        logging.error(f"extract_focus_tsv: {dof_tsv} non trovato")
        return None
    rows = []
    with open(dof_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            try:
                rows.append((row["x1"], row["x2"], row["x_focus"]))
            except KeyError:
                continue
    if not rows:
        logging.warning("extract_focus_tsv: nessun punto focale valido")
        return None
    rows = _remove_focus_outliers(rows, dx=dx)
    if not rows:
        logging.warning("extract_focus_tsv: nessun punto rimasto dopo rimozione outlier")
        return None
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write("x1\tx2\tx_focus\n")
        for r in rows:
            f.write(f"{r[0]}\t{r[1]}\t{r[2]}\n")
    logging.info(f"Focus TSV accurato: {len(rows)} punti → {out_path}")
    return out_path

# ─── Disk-guard utilities ─────────────────────────────────────────────────────

def _arange(start: float, stop: float, step: float) -> list[float]:
    result = []
    i = 0
    while True:
        val = start + i * step
        if val > stop + 1e-9:
            break
        result.append(val)
        i += 1
    return result


def _load_thorlabs_data() -> dict[str, dict]:
    """Carica spessori e offset da TSV Thorlabs (stessa logica di run.sh)."""
    data: dict[str, dict] = {}

    def _parse(filename: str, has_rotation: bool) -> None:
        path = PROJECT_ROOT / "lens_cutter" / "lens_data" / filename
        if not path.exists():
            return
        with open(path) as f:
            for line in f.readlines()[1:]:
                parts = line.split('\t')
                if len(parts) < 6:
                    continue
                lid  = parts[0].strip()
                tc   = float(parts[4])
                rot  = float(parts[7]) if has_rotation and len(parts) >= 8 else 0.0
                sign = -1.0 if abs(rot - 180.0) < 1e-6 else 1.0
                te   = float(parts[5])
                data[lid] = {'h': tc, 'offset': (tc - te) / 2.0 * sign if has_rotation else 0.0}

    _parse('thorlabs_biconvex.tsv', has_rotation=False)
    _parse('thorlabs_planoconvex.tsv', has_rotation=True)
    return data


def generate_pairs(cfg: dict, l1_id: str, l2_id: str,
                   mode: str = "lens") -> list[tuple[float, float]]:
    """Genera la lista completa di coppie (x1, x2) — stessa logica di run.sh."""
    import math
    thorlabs = _load_thorlabs_data()
    h1     = thorlabs.get(l1_id, {}).get('h', cfg.get('h1', 12.5))
    h2     = thorlabs.get(l2_id, {}).get('h', cfg.get('h2', 16.3))
    margin = 3.0 if mode == 'lens' else 1.0

    x_min = cfg['x_min']
    x_max = cfg['x_max']
    dx    = cfg['dx']

    pairs: list[tuple[float, float]] = []
    for x1 in _arange(x_min, x_max, dx):
        x2_min_collision = x1 + (h1 + h2) / 2.0 + margin
        x2_start_raw     = max(x1 + dx, x2_min_collision)
        x2_start = x_min + math.ceil(round((x2_start_raw - x_min) / dx, 8)) * dx
        for x2 in _arange(x2_start, x_max, dx):
            pairs.append((round(float(x1), 6), round(float(x2), 6)))
    return pairs


def estimate_lens_gb_per_pair(cfg: dict) -> float:
    """Stima GB per coppia per lens.root (conservativa: 20% hit rate, 75 B/hit)."""
    nx = len(_arange(cfg.get('source_x_min', -30.0), cfg.get('source_x_max', 30.0),
                     cfg.get('source_dx', 0.5)))
    ny = len(_arange(cfg.get('source_y_min', 0.0), cfg.get('source_y_max', 14.14),
                     cfg.get('source_dy', 0.5)))
    n_photons = cfg.get('lens_n_photons', 1000)
    return nx * ny * n_photons * 0.20 * 75 / 1e9


def run_lens_with_disk_guard(
    cfg: dict, cfg_path: Path, out_dir: Path, ssd_mount: Path,
    jobs: int, l1_id: str, l2_id: str, use_local: bool, dry_run: bool,
    focus_tsv: Path | None, disk_guard_gb: float, keep_lens: bool,
    min_hits: int = 5,
) -> Path | None:
    """Esegue lens + psf_extractor in super-batch reattivi al disco.

    Ogni super-batch simula un sottoinsieme di coppie, estrae subito psf_data
    parziale e cancella lens.root, tenendo il picco su disco sotto il budget.
    """
    import shutil as _shutil

    log   = out_dir / "pipeline.log"
    data  = out_dir / "data"

    all_pairs   = generate_pairs(cfg, l1_id, l2_id, mode="lens")
    total       = len(all_pairs)
    gb_per_pair = estimate_lens_gb_per_pair(cfg)

    avail_gb = _shutil.disk_usage(out_dir).free / 1e9
    safe_gb  = max(avail_gb - disk_guard_gb, 0.0)
    batch_size = max(1, int(safe_gb / max(gb_per_pair, 1e-9)))
    batch_size = min(batch_size, total)

    logging.info(
        f"[disk-guard] {total} coppie totali, ~{gb_per_pair:.3f} GB/coppia, "
        f"{avail_gb:.0f} GB liberi → batch_size iniziale={batch_size}"
    )

    nx = len(_arange(cfg.get('source_x_min', -30.0), cfg.get('source_x_max', 30.0),
                     cfg.get('source_dx', 0.5)))
    ny = len(_arange(cfg.get('source_y_min', 0.0), cfg.get('source_y_max', 14.14),
                     cfg.get('source_dy', 0.5)))
    n_runs_per_pair = nx * ny

    psf_partials:     list[Path] = []
    psf_dof_partials: list[Path] = []
    offset_cfg = 0
    offset_run = 0
    batch_idx  = 0
    pos        = 0

    while pos < total:
        # Ricalcola batch_size in base allo spazio corrente
        avail_gb = _shutil.disk_usage(out_dir).free / 1e9
        safe_gb  = max(avail_gb - disk_guard_gb, 0.0)
        if safe_gb < gb_per_pair and not dry_run:
            logging.error(
                f"[disk-guard] Spazio insufficiente ({avail_gb:.0f} GB liberi, "
                f"guardia={disk_guard_gb:.0f} GB). Interruzione."
            )
            return None
        batch_size = max(1, int(safe_gb / max(gb_per_pair, 1e-9)))
        batch_size = min(batch_size, total - pos)

        batch_pairs = all_pairs[pos : pos + batch_size]

        # Config temporaneo per il super-batch
        tmp_cfg = dict(cfg)
        tmp_cfg['pairs']                 = [list(p) for p in batch_pairs]
        tmp_cfg['config_id_offset'] = offset_cfg
        tmp_cfg['run_id_offset']    = offset_run
        tmp_cfg_path = data / f"superbatch_{batch_idx}.json"
        if not dry_run:
            tmp_cfg_path.write_text(json.dumps(tmp_cfg, indent=2))

        partial_lens = data / f"partial_lens_{batch_idx}.root"
        partial_psf  = data / f"partial_psf_data_{batch_idx}.root"

        logging.info(
            f"[disk-guard] Super-batch {batch_idx}: coppie {pos}–{pos+batch_size-1} "
            f"({batch_size} coppie, ~{batch_size*gb_per_pair:.1f} GB attesi)"
        )

        # Simula questo batch
        lens_root = run_simulation_step(
            "lens", tmp_cfg_path, out_dir, ssd_mount, jobs,
            l1_id, l2_id, use_local, dry_run, focus_tsv,
        )
        if lens_root is None:
            logging.error(f"[disk-guard] Simulazione lens fallita al super-batch {batch_idx}")
            return None

        # Sposta lens.root → partial_lens_K.root (shutil.move gestisce cross-device)
        if not dry_run and lens_root.exists() and lens_root != partial_lens:
            import shutil as _shutil_mv
            _shutil_mv.move(str(lens_root), str(partial_lens))
            # Il symlink data/lens.root punta al file appena rinominato: rimuovilo
            stale_link = data / "lens.root"
            if stale_link.is_symlink():
                stale_link.unlink()
        elif dry_run:
            partial_lens = Path(f"/dev/null/partial_lens_{batch_idx}.root")

        # Estrai psf_data parziale
        cmd = [find_binary("psf_extractor"), str(partial_lens), str(partial_psf), str(min_hits)]
        ok, _, _ = run_cmd(cmd, log, TIMEOUTS_SEC["analysis"], dry_run)
        if not ok:
            logging.error(f"[disk-guard] psf_extractor fallito al super-batch {batch_idx}")
            return None

        psf_partials.append(partial_psf)

        # Estrai psf_dof_data parziale dal partial_lens (Tecnica E embedded)
        partial_psf_dof  = data / f"partial_psf_dof_data_{batch_idx}.root"
        partial_res_tmp  = data / f"partial_res_tsv_{batch_idx}.tsv"
        cmd_dof = [
            find_binary("psf_dof_extractor"),
            str(partial_lens),
            str(partial_psf_dof),
            str(partial_res_tmp),
            str(cfg_path),
        ]
        ok_dof, _, _ = run_cmd(cmd_dof, log, TIMEOUTS_SEC["psf_dof_extractor"], dry_run)
        if not ok_dof:
            logging.warning(
                f"[disk-guard] psf_dof_extractor fallito al super-batch {batch_idx} — continuo"
            )
        else:
            psf_dof_partials.append(partial_psf_dof)

        # Libera il lens parziale
        if not keep_lens and not dry_run and partial_lens.exists():
            try:
                sz = partial_lens.stat().st_size / (1024 ** 3)
                logging.info(f"[disk-guard] Eliminato partial_lens_{batch_idx}.root ({sz:.1f} GB)")
                partial_lens.unlink()
            except Exception as e:
                logging.warning(f"[disk-guard] Impossibile eliminare {partial_lens}: {e}")

        offset_cfg += len(batch_pairs)
        offset_run += len(batch_pairs) * n_runs_per_pair
        pos        += batch_size
        batch_idx  += 1

    # Combina tutti i psf_data parziali
    psf_data = data / "psf_data.root"
    if dry_run:
        logging.info(f"[disk-guard] [DRY] hadd {psf_data} ← {batch_idx} parziali")
        logging.info(f"[disk-guard] [DRY] merge TSV parziali → {data / 'resolution_map.tsv'}")
        return psf_data

    if len(psf_partials) == 1:
        psf_partials[0].rename(psf_data)
    else:
        hadd_cmd = ["hadd", "-fk", str(psf_data)] + [str(p) for p in psf_partials]
        res = subprocess.run(hadd_cmd, capture_output=True, text=True)
        if res.returncode != 0:
            logging.error(f"[disk-guard] hadd fallito:\n{res.stderr}")
            return None
        for p in psf_partials:
            try:
                p.unlink()
            except Exception:
                pass

    # Combina psf_dof parziali (Tecnica E)
    psf_dof_data_final = data / "psf_dof_data.root"
    if psf_dof_partials and not dry_run:
        if len(psf_dof_partials) == 1:
            psf_dof_partials[0].rename(psf_dof_data_final)
        else:
            hadd_dof = ["hadd", "-fk", str(psf_dof_data_final)] + [str(p) for p in psf_dof_partials]
            res_dof  = subprocess.run(hadd_dof, capture_output=True, text=True)
            if res_dof.returncode != 0:
                logging.error(f"[disk-guard] hadd psf_dof fallito:\n{res_dof.stderr}")
            else:
                for p in psf_dof_partials:
                    try: p.unlink()
                    except Exception: pass
        logging.info(f"[disk-guard] psf_dof_data.root finale: {psf_dof_data_final}")
    elif dry_run and psf_dof_partials:
        logging.info(f"[disk-guard] [DRY] hadd {psf_dof_data_final} ← {len(psf_dof_partials)} parziali psf_dof")

    # Merge TSV parziali di psf_dof in resolution_map.tsv finale
    res_tsv_final = data / "resolution_map.tsv"
    partial_tsvs = [data / f"partial_res_tsv_{i}.tsv"
                    for i in range(batch_idx)
                    if (data / f"partial_res_tsv_{i}.tsv").exists()]
    if partial_tsvs and not dry_run:
        with open(res_tsv_final, "w") as fout:
            header_written = False
            for p in partial_tsvs:
                with open(p) as fin:
                    lines = fin.readlines()
                if not header_written and lines:
                    fout.writelines(lines[:1])
                    header_written = True
                fout.writelines(lines[1:])
        logging.info(f"[disk-guard] resolution_map.tsv merged: {res_tsv_final}")

    # Pulizia TSV temporanei di psf_dof
    for i in range(batch_idx):
        tmp_tsv = data / f"partial_res_tsv_{i}.tsv"
        if tmp_tsv.exists():
            tmp_tsv.unlink()

    # Pulizia config temporanei
    for i in range(batch_idx):
        tmp = data / f"superbatch_{i}.json"
        if tmp.exists():
            tmp.unlink()

    logging.info(f"[disk-guard] psf_data.root finale: {psf_data}")
    return psf_data


def _delete_ssd_file(path: Path | None, data_dir: Path | None,
                     label: str, ssd_mount: Path, dry_run: bool) -> None:
    """Cancella un file ROOT dall'SSD e tenta di rimuovere la dir run ora vuota."""
    if dry_run or path is None:
        return
    try:
        if not str(path.resolve()).startswith(str(ssd_mount)):
            return
    except Exception:
        return
    if path.exists():
        try:
            sz = path.stat().st_size / (1024 ** 3)
            logging.info(f"[ssd-cleanup] Eliminato {label} ({sz:.1f} GB): {path.name}")
            path.unlink()
        except Exception as e:
            logging.warning(f"[ssd-cleanup] Impossibile eliminare {path}: {e}")
    if data_dir is not None:
        link = data_dir / path.name
        if link.is_symlink() and not link.exists():
            try:
                link.unlink()
            except Exception:
                pass
    run_dir = path.parent
    try:
        run_dir.rmdir()
        logging.info(f"[ssd-cleanup] Rimossa dir run vuota: {run_dir.name}")
    except OSError:
        pass


# ─── Simulation step launcher ─────────────────────────────────────────────────

def run_simulation_step(mode: str, cfg_path: Path, out_dir: Path,
                        ssd_mount: Path, jobs: int, l1_id: str, l2_id: str,
                        use_local: bool, dry_run: bool,
                        focus_tsv: Path | None = None) -> Path | None:
    out_names = {"opt": "events", "lens": "lens", "dof": "focal"}
    prefix = out_names[mode]

    run_target = "local" if use_local else "ssd"
    cmd: list = [
        find_binary("run_sh"), mode, run_target,
        "--jobs",   str(jobs),
        "--config", str(cfg_path),
    ]
    if not use_local:
        cmd += ["--ssd-mount", str(ssd_mount)]
    if l1_id:
        cmd += ["--l1-id", l1_id]
    if l2_id:
        cmd += ["--l2-id", l2_id]
    if focus_tsv and mode in ("opt", "lens"):
        cmd += ["--focus-tsv", str(focus_tsv)]

    log_path   = out_dir / "pipeline.log"
    step_start = time.time()
    ok, stdout, _ = run_cmd(cmd, log_path, TIMEOUTS_SEC[mode], dry_run)
    if not ok:
        return None

    if dry_run:
        return Path(f"/dev/null/{prefix}.root")

    root_path = find_ssd_root_from_stdout(stdout, prefix)
    if root_path is None and not use_local:
        root_path = find_latest_ssd_root(ssd_mount, prefix, after=step_start)
    if root_path is None:
        logging.error(f"File ROOT '{prefix}.root' non trovato dopo step '{mode}'")
        return None

    # Symlink locale in data/ per riferimento senza duplicare
    link = out_dir / "data" / f"{prefix}.root"
    link.parent.mkdir(parents=True, exist_ok=True)
    if not link.exists():
        try:
            link.symlink_to(root_path)
        except Exception:
            pass

    return root_path

# ─── Pipeline ────────────────────────────────────────────────────────────────

def run_pipeline(cfg_path: Path, out_dir: Path, ssd_mount: Path, jobs: int,
                 l1_id: str, l2_id: str, ap: dict,
                 dry_run: bool, use_local: bool,
                 mobile_focus: bool,
                 keep_lens: bool,
                 skip_opt: bool, skip_lens: bool,
                 skip_dof: bool, skip_psf_dof: bool,
                 disk_guard_gb: float = 200.0) -> dict | None:

    out   = {}
    log   = out_dir / "pipeline.log"
    plots = out_dir / "plots"
    data  = out_dir / "data"
    plots.mkdir(parents=True, exist_ok=True)
    data.mkdir(parents=True, exist_ok=True)

    try:
        _run_cfg = json.loads(cfg_path.read_text())
    except Exception:
        _run_cfg = {}

    if not use_local:
        ok, fid, msg = check_ssd_mount(ssd_mount)
        if not ok:
            logging.error(f"Pre-flight FAIL [{fid}]: {msg}")
            return None

    # Fuoco mobile: esegue dof+dof_map come pre-pipeline per ottenere i fuochi reali
    focus_tsv: Path | None = None
    if mobile_focus:
        dof_tsv_pre = data / "dof_map.tsv"
        if skip_dof and (dry_run or dof_tsv_pre.exists()):
            logging.info("[mobile-focus] Uso dof_map.tsv esistente (--skip-dof)")
        else:
            logging.info("[mobile-focus] Pre-pipeline: lancio dof per stimare il fuoco reale...")
            focal_root_pre = run_simulation_step(
                "dof", cfg_path, out_dir, ssd_mount, jobs,
                l1_id, l2_id, use_local, dry_run,
            )
            if focal_root_pre is None:
                return None
            dp  = ap["dof_map"]
            cmd = [
                find_binary("dof_map"),
                "--input",         str(focal_root_pre),
                "--config",        str(cfg_path),
                "--core-fraction", str(dp["core_fraction"]),
                "--m-target",      str(dp["m_target"]),
                "--tsv",           str(dof_tsv_pre),
                "--jobs",          str(jobs),
                "--output",        str(plots),
            ]
            ok, _, _ = run_cmd(cmd, log, TIMEOUTS_SEC["analysis"], dry_run)
            if not ok:
                return None
            skip_dof = True  # dof già eseguito, non rieseguire negli step 7-8
        focus_tsv = extract_focus_tsv(dof_tsv_pre, out_dir / "focus_accurate.tsv", dry_run,
                                       dx=_run_cfg.get("dx", 3.0))
        if focus_tsv is None:
            logging.warning("[mobile-focus] Nessun fuoco estratto — continuo con fuoco fisso")

    # ── Step 1: opt ───────────────────────────────────────────────────────────
    if skip_opt:
        events_root = data / "events.root"
        if not dry_run and not events_root.exists():
            logging.error("--skip-opt richiede data/events.root esistente")
            return None
        logging.info(f"[opt] SKIP — uso {events_root}")
    else:
        events_root = run_simulation_step(
            "opt", cfg_path, out_dir, ssd_mount, jobs,
            l1_id, l2_id, use_local, dry_run, focus_tsv,
        )
        if events_root is None:
            return None
    out["events"] = events_root

    # ── Step 2: plot2D → efficiency2D.png ────────────────────────────────────
    eff_png = plots / "efficiency2D.png"
    cmd = [find_binary("plot2D"),
           "--input",  str(events_root),
           "--config", str(cfg_path),
           "--output", str(eff_png)]
    ok, _, _ = run_cmd(cmd, log, TIMEOUTS_SEC["analysis"], dry_run)
    if ok:
        out["efficiency_png"] = eff_png
    else:
        logging.warning("plot2D fallito — continuo")

    # events.root sull'SSD non è più necessario dopo plot2D
    if not skip_opt:
        _delete_ssd_file(events_root, data, "events.root", ssd_mount, dry_run)

    # ── Steps 3–4: lens + psf_extractor ──────────────────────────────────────
    psf_data = data / "psf_data.root"
    # lens_root_for_psfdof: lens.root tenuto vivo per psf_dof_extractor (Tecnica E)
    lens_root_for_psfdof: Path | None = None

    if skip_lens:
        if not dry_run and not psf_data.exists():
            logging.error("--skip-lens richiede data/psf_data.root esistente")
            return None
        logging.info(f"[lens+psf_extractor] SKIP — uso {psf_data}")
    elif disk_guard_gb > 0:
        # Super-batch reattivo al disco (psf_dof estratto per-batch internamente)
        if run_lens_with_disk_guard(
            cfg=_run_cfg, cfg_path=cfg_path, out_dir=out_dir,
            ssd_mount=ssd_mount, jobs=jobs, l1_id=l1_id, l2_id=l2_id,
            use_local=use_local, dry_run=dry_run, focus_tsv=focus_tsv,
            disk_guard_gb=disk_guard_gb, keep_lens=keep_lens,
            min_hits=ap["psf_extractor"]["min_hits"],
        ) is None:
            return None
    else:
        lens_root = run_simulation_step(
            "lens", cfg_path, out_dir, ssd_mount, jobs,
            l1_id, l2_id, use_local, dry_run, focus_tsv,
        )
        if lens_root is None:
            return None
        out["lens"] = lens_root

        min_hits = ap["psf_extractor"]["min_hits"]
        cmd = [find_binary("psf_extractor"), str(lens_root), str(psf_data), str(min_hits)]
        ok, _, _ = run_cmd(cmd, log, TIMEOUTS_SEC["analysis"], dry_run)
        if not ok:
            return None

        # Mantieni lens.root vivo per psf_dof_extractor (Tecnica E); cancellazione dopo step 11
        lens_root_for_psfdof = lens_root

    out["psf_data"] = psf_data

    # ── Step 5: q_map ─────────────────────────────────────────────────────────
    q_tsv = data / "q_map.tsv"
    qp    = ap["q_map"]
    cmd   = [
        find_binary("q_map"),
        "--psf",        str(psf_data),
        "--config",     str(cfg_path),
        "--n-tracks",   str(qp["n_tracks"]),
        "--dt",         str(qp["dt"]),
        "--min-hits",   str(qp["min_hits"]),
        "--trace-frac", str(qp["trace_frac"]),
        "--tsv",        str(q_tsv),
        "--jobs",       str(jobs),
        "--output",     str(plots / "q_map.png"),
    ]
    if qp.get("dist_to_target"):
        cmd.append("--dist-to-target")
    ok, _, _ = run_cmd(cmd, log, TIMEOUTS_SEC["analysis"], dry_run)
    if not ok:
        return None
    out["q_tsv"] = q_tsv

    # ── Step 5b: coverage_map ─────────────────────────────────────────────────
    cov_tsv = data / "coverage_map.tsv"
    cmd_cov = [
        find_binary("q_map"),
        "--psf",        str(psf_data),
        "--config",     str(cfg_path),
        "--n-tracks",   str(qp["n_tracks"]),
        "--dt",         str(qp["dt"]),
        "--min-hits",   str(qp["min_hits"]),
        "--trace-frac", str(qp["trace_frac"]),
        "--coverage",
        "--tsv",        str(cov_tsv),
        "--jobs",       str(jobs),
        "--output",     str(plots / "coverage_map.png"),
    ]
    ok_cov, _, _ = run_cmd(cmd_cov, log, TIMEOUTS_SEC["analysis"], dry_run)
    if ok_cov:
        out["coverage_tsv"] = cov_tsv
    else:
        logging.warning("coverage_map fallito — proseguo comunque")

    # ── Step 6: chi2_map ──────────────────────────────────────────────────────
    chi2_tsv = data / "chi2_map.tsv"
    cp       = ap["chi2_map"]
    cmd      = [
        find_binary("chi2_map"),
        "--psf",      str(psf_data),
        "--config",   str(cfg_path),
        "--min-hits", str(cp["min_hits"]),
        "--p-low",    str(cp["p_low"]),
        "--p-high",   str(cp["p_high"]),
        "--tsv",      str(chi2_tsv),
        "--output",   str(plots / "chi2_map.png"),
    ]
    if cp.get("adaptive_target"):
        cmd.append("--adaptive-target")
    ok, _, _ = run_cmd(cmd, log, TIMEOUTS_SEC["analysis"], dry_run)
    if not ok:
        return None
    out["chi2_tsv"] = chi2_tsv

    # ── Step 6b: cor_map ──────────────────────────────────────────────────────
    corr_tsv = data / "corr_map.tsv"
    cmd_corr = [
        find_binary("chi2_map"),
        "--psf",      str(psf_data),
        "--config",   str(cfg_path),
        "--min-hits", str(cp["min_hits"]),
        "--p-low",    str(cp["p_low"]),
        "--p-high",   str(cp["p_high"]),
        "--corr-map",
        "--tsv",      str(corr_tsv),
        "--jobs",     str(jobs),
        "--output",   str(plots / "cor_map.png"),
    ]
    ok_corr, _, _ = run_cmd(cmd_corr, log, TIMEOUTS_SEC["analysis"], dry_run)
    if not ok_corr:
        logging.warning("cor_map fallito — proseguo comunque")

    # ── Steps 7–8: dof + dof_map ──────────────────────────────────────────────
    dof_tsv = data / "dof_map.tsv"
    if skip_dof:
        if not dry_run and not dof_tsv.exists():
            logging.error("--skip-dof richiede data/dof_map.tsv esistente")
            return None
        logging.info(f"[dof+dof_map] SKIP — uso {dof_tsv}")
    else:
        focal_root = run_simulation_step(
            "dof", cfg_path, out_dir, ssd_mount, jobs,
            l1_id, l2_id, use_local, dry_run,
        )
        if focal_root is None:
            return None
        out["focal"] = focal_root

        dp  = ap["dof_map"]
        cmd = [
            find_binary("dof_map"),
            "--input",         str(focal_root),
            "--config",        str(cfg_path),
            "--core-fraction", str(dp["core_fraction"]),
            "--m-target",      str(dp["m_target"]),
            "--tsv",           str(dof_tsv),
            "--jobs",          str(jobs),
            "--output",        str(plots),
        ]
        ok, _, _ = run_cmd(cmd, log, TIMEOUTS_SEC["analysis"], dry_run)
        if not ok:
            return None
        # focal.root sull'SSD non è più necessario dopo dof_map
        _delete_ssd_file(focal_root, data, "focal.root", ssd_mount, dry_run)

    out["dof_tsv"] = dof_tsv

    # ── Steps 9–10: psf_dof_extractor (legge PsfDofRuns embedded in lens.root) ──
    # La simulazione psf-dof separata è eliminata (Tecnica E): i dati sono già in lens.root.
    psf_dof_data = data / "psf_dof_data.root"
    res_tsv      = data / "resolution_map.tsv"
    psf_dof_root_raw: Path | None = None  # file grezzo per resolution_map

    if skip_psf_dof:
        if not dry_run and not psf_dof_data.exists():
            logging.error("--skip-psf-dof richiede data/psf_dof_data.root esistente")
            return None
        logging.info(f"[psf_dof_extractor] SKIP — uso {psf_dof_data}")
    elif lens_root_for_psfdof is not None and (dry_run or lens_root_for_psfdof.exists()):
        # Caso non-disk-guard: estrai direttamente da lens.root
        psf_dof_root_raw = lens_root_for_psfdof
        cmd = [
            find_binary("psf_dof_extractor"),
            str(lens_root_for_psfdof),
            str(psf_dof_data),
            str(res_tsv),
            str(cfg_path),
        ]
        ok, _, _ = run_cmd(cmd, log, TIMEOUTS_SEC["psf_dof_extractor"], dry_run)
        if not ok:
            logging.warning("[psf_dof_extractor] fallito — proseguo senza dati DoF")
            psf_dof_root_raw = None
    elif not dry_run and psf_dof_data.exists():
        # Caso disk-guard: psf_dof_data già costruito dalle super-batch
        logging.info(f"[psf_dof_extractor] psf_dof_data.root dal disk-guard: {psf_dof_data}")
    else:
        logging.info("[psf_dof_extractor] SKIP — lens.root non disponibile e nessun psf_dof_data.root")

    out["psf_dof_data"] = psf_dof_data
    out["res_tsv"]      = res_tsv

    # ── Step 11: resolution_map (usa lens.root con PsfDofConfigs+PsfDofRuns embedded) ──
    _focus_warn_args: list[str] = []
    _dof_tsv_path = dof_tsv if dof_tsv is not None else (data / "dof_map.tsv")
    if _dof_tsv_path.exists():
        _focus_warn_args = ["--focus-warn-tsv", str(_dof_tsv_path)]

    if psf_dof_root_raw is not None and (dry_run or psf_dof_root_raw.exists()):
        cmd = [
            find_binary("resolution_map"),
            "--input",  str(psf_dof_root_raw),
            "--config", str(cfg_path),
            "--tsv",    str(res_tsv),
            "--jobs",   str(jobs),
            "--output", str(plots),
        ] + _focus_warn_args
        ok, _, _ = run_cmd(cmd, log, TIMEOUTS_SEC["resolution_map"], dry_run)
        if not ok:
            logging.warning("resolution_map fallito — i TSV sono comunque disponibili in data/")
    elif dry_run or res_tsv.exists():
        cmd = [
            find_binary("resolution_map"),
            "--from-tsv", str(res_tsv),
            "--config",   str(cfg_path),
            "--output",   str(plots),
        ] + _focus_warn_args
        ok, _, _ = run_cmd(cmd, log, TIMEOUTS_SEC["resolution_map"], dry_run)
        if not ok:
            logging.warning("resolution_map (--from-tsv) fallito — TSV disponibile in data/")
    else:
        logging.info("[resolution_map] SKIP — psf_dof_root_raw e resolution_map.tsv non disponibili")

    # Elimina lens.root dopo resolution_map (Tecnica E: sostituisce eliminazione psf_dof.root)
    if not keep_lens:
        _delete_ssd_file(lens_root_for_psfdof, data, "lens.root", ssd_mount, dry_run)

    return out

# ─── Simulation log CSV ───────────────────────────────────────────────────────

def append_to_log(cfg: dict, l1_id: str, l2_id: str,
                  detector_mode: str,
                  storage_type: str, storage_path: str,
                  out_dir: Path, notes: str) -> None:
    LOG_CSV.parent.mkdir(parents=True, exist_ok=True)
    write_header = not LOG_CSV.exists()
    with open(LOG_CSV, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=LOG_FIELDS)
        if write_header:
            writer.writeheader()
        writer.writerow({
            "date":              datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "l1_id":             l1_id,
            "l2_id":             l2_id,
            "detector_mode":     detector_mode,
            "x_min":             cfg.get("x_min", "unknown"),
            "x_max":             cfg.get("x_max", "unknown"),
            "dx":                cfg.get("dx", "unknown"),
            "n_photons":         cfg.get("n_photons", "unknown"),
            "lens_n_photons":    cfg.get("lens_n_photons", "unknown"),
            "dof_n_photons":     cfg.get("dof_n_photons", "unknown"),
            "psf_dof_n_photons": cfg.get("psf_dof_n_photons", "unknown"),
            "lens_gap_margin":   cfg.get("lens_gap_margin", "unknown"),
            "storage_type":      storage_type,
            "storage_path":      storage_path,
            "output_dir":        str(out_dir),
            "notes":             notes,
        })
    logging.info(f"Log aggiornato: {LOG_CSV}")

# ─── Main ─────────────────────────────────────────────────────────────────────

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Pipeline RIPTIDE completa per una singola coppia di lenti",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Esempi:
  # Run completo su SSD (ROOT pesanti su /mnt/external_ssd)
  python3 scripts/lens_runner.py --l1-id LA4464 --l2-id LA4464R

  # Run locale, detector mobile sul piano focale
  python3 scripts/lens_runner.py --l1-id LA4464 --l2-id LA4464R --local --mobile-focus

  # Dry-run: mostra comandi senza eseguire
  python3 scripts/lens_runner.py --l1-id LA4464 --l2-id LA4464R --local --dry-run

  # Resume: salta opt e lens già completati
  python3 scripts/lens_runner.py --l1-id LA4464 --l2-id LA4464R --local --skip-opt --skip-lens

  # Storage su SSD esterno custom
  python3 scripts/lens_runner.py --l1-id LA4464 --l2-id LA4464R --ssd-mount /mnt/runs
""",
    )
    req = parser.add_argument_group("obbligatori")
    req.add_argument("--l1-id", type=str, required=True, help="ID lente L1 Thorlabs (es. LA4464)")
    req.add_argument("--l2-id", type=str, required=True, help="ID lente L2 Thorlabs (es. LA4464R)")

    det = parser.add_argument_group("modalità detector")
    det.add_argument("--mobile-focus", action="store_true",
                     help="Fuoco mobile: esegue dof+dof_map prima della pipeline per ottenere il piano focale reale")

    sto = parser.add_argument_group("storage")
    sto.add_argument("--ssd-mount", type=str, default="/mnt/external_ssd",
                     help="Path SSD esterno per ROOT files pesanti (default: /mnt/external_ssd)")
    sto.add_argument("--local", action="store_true",
                     help="Salva ROOT files in locale (no SSD mount)")

    cfg = parser.add_argument_group("configurazione")
    cfg.add_argument("--config",          type=str, default=None,
                     help="Override config.json (default: config/config.json)")
    cfg.add_argument("--analysis-params", type=str, default=None,
                     help="Override analysis_params.json (default: config/analysis_params.json)")
    cfg.add_argument("--jobs",            type=int, default=os.cpu_count() or 1,
                     help="Job paralleli per run.sh (default: nproc)")

    opt = parser.add_argument_group("opzioni avanzate")
    opt.add_argument("--dry-run",       action="store_true", help="Stampa comandi senza eseguire")
    opt.add_argument("--keep-lens",     action="store_true", help="Non eliminare lens.root dopo psf_extractor")
    opt.add_argument("--skip-opt",      action="store_true", help="Salta opt (usa data/events.root esistente)")
    opt.add_argument("--skip-lens",     action="store_true", help="Salta lens+psf_extractor (usa data/psf_data.root)")
    opt.add_argument("--skip-dof",      action="store_true", help="Salta dof+dof_map (usa data/dof_map.tsv)")
    opt.add_argument("--skip-psf-dof",  action="store_true", help="Salta psf-dof+psf_dof_extractor (usa data/psf_dof_data.root)")
    opt.add_argument("--disk-guard-gb", type=float, default=200.0,
                     help="Super-batch reattivo: estr. psf_data e cancella lens.root mantenendo "
                          "almeno N GB liberi su disco (default: 200, minimo forzato: 200)")
    opt.add_argument("--notes",         type=str, default="",  help="Note aggiuntive nel CSV di log")

    args = parser.parse_args()
    if args.disk_guard_gb > 0:
        args.disk_guard_gb = max(args.disk_guard_gb, 1.0)

    pair_tag  = f"{args.l1_id}_{args.l2_id}"
    out_dir   = PROJECT_ROOT / "output" / "lens_simulations" / pair_tag
    ssd_mount = Path(args.ssd_mount)

    setup_logging(out_dir)
    t_start = time.time()
    logging.info(f"=== RIPTIDE Lens Runner — {pair_tag} ===")
    logging.info(
        f"dry_run={args.dry_run}  mobile_focus={args.mobile_focus}  "
        f"storage={'local' if args.local else ssd_mount}  jobs={args.jobs}"
    )

    cfg_path = Path(args.config) if args.config else PROJECT_ROOT / "config" / "config.json"
    base_cfg = load_base_config(args.config)
    ap       = load_analysis_params(args.analysis_params)

    pf_failures = preflight_checks()
    if pf_failures:
        logging.warning(
            f"Pre-flight WARNINGS: {pf_failures} — "
            "usa autonomous_optimizer.py per auto-fix automatico"
        )

    if not args.local and not args.dry_run and not ssd_mount.exists():
        logging.error(f"SSD mount non trovato: {ssd_mount}. Usa --local o --ssd-mount PATH corretto.")
        return 1

    detector_mode = "mobile" if args.mobile_focus else "fixed"
    storage_type  = "local" if args.local else "ssd"
    storage_path  = (
        str(PROJECT_ROOT / "output")
        if args.local
        else str(ssd_mount / "riptide" / "runs")
    )

    result = run_pipeline(
        cfg_path         = cfg_path,
        out_dir          = out_dir,
        ssd_mount        = ssd_mount,
        jobs             = args.jobs,
        l1_id            = args.l1_id,
        l2_id            = args.l2_id,
        ap               = ap,
        dry_run          = args.dry_run,
        use_local        = args.local,
        mobile_focus     = args.mobile_focus,
        keep_lens        = args.keep_lens,
        skip_opt         = args.skip_opt,
        skip_lens        = args.skip_lens,
        skip_dof         = args.skip_dof,
        skip_psf_dof     = args.skip_psf_dof,
        disk_guard_gb    = args.disk_guard_gb,
    )

    elapsed = time.time() - t_start
    if result is None:
        logging.error(f"Pipeline FALLITA dopo {elapsed/60:.1f} min")
        return 1

    logging.info(f"Pipeline COMPLETATA in {elapsed/60:.1f} min")
    logging.info(f"Output plots : {out_dir / 'plots'}")
    logging.info(f"Output data  : {out_dir / 'data'}")

    # Salva simulation_summary.json
    summary = {
        "date":              datetime.datetime.now().isoformat(timespec="seconds"),
        "l1_id":             args.l1_id,
        "l2_id":             args.l2_id,
        "detector_mode":     detector_mode,
        "storage_type":      storage_type,
        "storage_path":      storage_path,
        "output_dir":        str(out_dir),
        "config":            base_cfg,
        "analysis_params":   ap,
        "duration_min":      round(elapsed / 60, 2),
        "notes":             args.notes,
    }
    (out_dir / "simulation_summary.json").write_text(json.dumps(summary, indent=2))

    if not args.dry_run:
        append_to_log(
            cfg           = base_cfg,
            l1_id         = args.l1_id,
            l2_id         = args.l2_id,
            detector_mode = detector_mode,
            storage_type  = storage_type,
            storage_path  = storage_path,
            out_dir       = out_dir,
            notes         = args.notes,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
