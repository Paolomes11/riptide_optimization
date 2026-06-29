#!/usr/bin/env python3
"""
Validazione Tecnica E: confronta PsfDofRuns embedded (lens_simulation_main)
vs standalone (psf_dof_scan_main) spot-per-spot.

Uso:
  python3 scripts/validate_tecnica_e.py embedded_lens.root standalone_psf_dof.root
  python3 scripts/validate_tecnica_e.py embedded.root standalone.root --threshold 3 --csv report.csv
"""
import argparse
import sys
import numpy as np

try:
    import uproot
except ImportError:
    sys.exit("uproot non trovato — pip install uproot")

FIELDS = ["mu_y", "mu_z", "sigma_y", "sigma_z", "cov_y_dy", "cov_z_dz"]


def load_runs(path: str) -> dict:
    """
    Restituisce dict (x1, x2, x_source, y_source) → {field: value, n_hits: value}.
    Usa (x1,x2,x_source,y_source) arrotondati a 3 cifre come chiave per match cross-run.
    """
    with uproot.open(path) as f:
        configs = f["PsfDofConfigs"].arrays(["config_id", "x1", "x2"], library="np")
        runs = f["PsfDofRuns"].arrays(
            ["config_id", "x_source", "y_source", "n_hits"] + FIELDS, library="np"
        )

    cid_to_x = {
        int(configs["config_id"][i]): (round(float(configs["x1"][i]), 3),
                                        round(float(configs["x2"][i]), 3))
        for i in range(len(configs["config_id"]))
    }

    result = {}
    for i in range(len(runs["config_id"])):
        x1, x2 = cid_to_x.get(int(runs["config_id"][i]), (None, None))
        if x1 is None:
            continue
        key = (x1, x2,
               round(float(runs["x_source"][i]), 3),
               round(float(runs["y_source"][i]), 3))
        result[key] = {f: float(runs[f][i]) for f in FIELDS}
        result[key]["n_hits"] = float(runs["n_hits"][i])
    return result


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Confronta PsfDofRuns embedded (Tecnica E) vs standalone"
    )
    ap.add_argument("embedded",   help="lens.root prodotto da lens_simulation_main")
    ap.add_argument("standalone", help="psf_dof.root prodotto da psf_dof_scan_main")
    ap.add_argument("--threshold", type=float, default=3.0,
                    help="Soglia outlier in unità shot-noise (default 3.0)")
    ap.add_argument("--outlier-frac", type=float, default=0.05,
                    help="Frazione massima outlier ammessa (default 0.05 = 5%%)")
    ap.add_argument("--csv", help="Salva report dettagliato per-spot in CSV")
    args = ap.parse_args()

    print(f"Carico embedded:   {args.embedded}")
    emb = load_runs(args.embedded)
    print(f"Carico standalone: {args.standalone}")
    std = load_runs(args.standalone)

    common = set(emb) & set(std)
    n_common = len(common)
    print(f"\nSpot confrontabili: {n_common} / emb={len(emb)} std={len(std)}")

    if n_common < 10:
        print("ATTENZIONE: meno di 10 spot in comune — confronto non statisticamente robusto.")

    if n_common == 0:
        print("ERRORE: nessun match — verificare che i run usino lo stesso focus TSV e config.")
        return 1

    rows: list[dict] = []
    for key in sorted(common):
        e = emb[key]
        s = std[key]
        n = (e["n_hits"] + s["n_hits"]) / 2.0
        row: dict = {"x1": key[0], "x2": key[1], "x_source": key[2], "y_source": key[3]}
        for f in FIELDS:
            delta = abs(e[f] - s[f])
            ref   = (abs(e[f]) + abs(s[f])) / 2.0
            # shot noise relativo per sigma-like quantities: sigma/sqrt(N)
            # usiamo sqrt(2/N) per avere margine conservativo
            shot  = max(ref, 1e-9) * np.sqrt(2.0 / max(n, 1.0))
            row[f"{f}_emb"]   = e[f]
            row[f"{f}_std"]   = s[f]
            row[f"{f}_delta"] = e[f] - s[f]
            row[f"{f}_norm"]  = delta / shot if shot > 1e-15 else 0.0
        rows.append(row)

    mean_n = float(np.mean([r.get("n_hits_emb", emb[k]["n_hits"]) for k in common
                             for r in [{"n_hits_emb": emb[k]["n_hits"]}]]))
    mean_n = float(np.mean([emb[k]["n_hits"] for k in common]))
    print(f"N_hits medio per spot: {mean_n:.0f}")
    if mean_n < 1000:
        print(f"ATTENZIONE: N={mean_n:.0f} < 1000 — rumore MC elevato, confrontare con N≥1000.")

    # Report per campo
    print(f"\n{'Campo':12s}  {'mean_norm':>10s}  {'max_norm':>9s}  "
          f"{'outliers>thr':>14s}  {'esito':>6s}")
    print("-" * 65)

    all_ok = True
    for f in FIELDS:
        norms = np.array([r[f"{f}_norm"] for r in rows])
        frac  = float(np.mean(norms > args.threshold))
        ok    = frac <= args.outlier_frac
        all_ok &= ok
        status = "OK  " if ok else "FAIL"
        print(f"{f:12s}  {np.mean(norms):10.2f}  {np.max(norms):9.2f}  "
              f"{100*frac:13.1f}%  [{status}]")

    print()
    if all_ok:
        print("RISULTATO: Tecnica E VALIDATA — momenti embedded coerenti con standalone")
    else:
        print("RISULTATO: VALIDAZIONE FALLITA — verificare implementazione Tecnica E")

    if args.csv:
        import csv
        fieldnames = (["x1", "x2", "x_source", "y_source"]
                      + [f"{f}_{s}" for f in FIELDS for s in ("emb", "std", "delta", "norm")])
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)
        print(f"Report CSV salvato: {args.csv}")

    return 0 if all_ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
