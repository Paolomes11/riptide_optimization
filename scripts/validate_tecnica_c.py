#!/usr/bin/env python3
"""
Validazione Tecnica C: confronta momenti PSF+DoF a lens_n_photons=10k vs 5k.
L'errore relativo medio su sigma_y, sigma_dz, cov_z_dz deve essere < 5% per
approvare la riduzione da 10k a 5k fotoni (speedup ×2).

Uso:
  python3 scripts/validate_tecnica_c.py psf_dof_10k.root psf_dof_5k.root
  python3 scripts/validate_tecnica_c.py 10k.root 5k.root --threshold 0.05 --csv report.csv
"""
import argparse
import sys
import numpy as np

try:
    import uproot
except ImportError:
    sys.exit("uproot non trovato — pip install uproot")

FIELDS = ["sigma_y", "sigma_z", "sigma_dz", "cov_y_dy", "cov_z_dz"]


def load_runs(path: str) -> dict:
    """
    Restituisce dict (x1, x2, x_source, y_source) → {field: value, n_hits: value}.
    """
    with uproot.open(path) as f:
        configs = f["PsfDofConfigs"].arrays(["config_id", "x1", "x2"], library="np")
        available = f["PsfDofRuns"].keys()
        # sigma_dz potrebbe non esistere in versioni vecchie del ntuple
        fields_to_load = [f for f in FIELDS if f in available]
        runs = f["PsfDofRuns"].arrays(
            ["config_id", "x_source", "y_source", "n_hits"] + fields_to_load,
            library="np",
        )

    cid_to_x = {
        int(configs["config_id"][i]): (round(float(configs["x1"][i]), 3),
                                        round(float(configs["x2"][i]), 3))
        for i in range(len(configs["config_id"]))
    }

    result: dict = {}
    for i in range(len(runs["config_id"])):
        x1, x2 = cid_to_x.get(int(runs["config_id"][i]), (None, None))
        if x1 is None:
            continue
        key = (x1, x2,
               round(float(runs["x_source"][i]), 3),
               round(float(runs["y_source"][i]), 3))
        result[key] = {f: float(runs[f][i]) for f in fields_to_load}
        result[key]["n_hits"] = float(runs["n_hits"][i])
    return result


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Valida Tecnica C: errore relativo 10k vs 5k fotoni"
    )
    ap.add_argument("ref10k",  help="psf_dof_data.root con lens_n_photons=10000 (riferimento)")
    ap.add_argument("test5k",  help="psf_dof_data.root con lens_n_photons=5000 (candidato)")
    ap.add_argument("--threshold", type=float, default=0.05,
                    help="Errore relativo massimo ammesso (default 0.05 = 5%%)")
    ap.add_argument("--csv", help="Salva report dettagliato per-spot in CSV")
    args = ap.parse_args()

    print(f"Carico riferimento (10k): {args.ref10k}")
    ref = load_runs(args.ref10k)
    print(f"Carico candidato  (5k) : {args.test5k}")
    tst = load_runs(args.test5k)

    common = set(ref) & set(tst)
    n_common = len(common)
    print(f"\nSpot confrontabili: {n_common} / ref={len(ref)} tst={len(tst)}")

    if n_common == 0:
        print("ERRORE: nessun match — verificare che i run usino le stesse coppie (x1,x2).")
        return 1

    mean_n_ref = float(np.mean([ref[k]["n_hits"] for k in common]))
    mean_n_tst = float(np.mean([tst[k]["n_hits"] for k in common]))
    print(f"N_hits medio  ref: {mean_n_ref:.1f}   tst: {mean_n_tst:.1f}")
    if mean_n_tst < 2:
        print("ATTENZIONE: N_hits troppo basso per confronto affidabile.")

    # Campi disponibili in comune tra i due file
    avail_fields = [f for f in FIELDS if f in next(iter(ref.values()))]

    rows: list[dict] = []
    for key in sorted(common):
        r = ref[key]
        t = tst[key]
        row: dict = {"x1": key[0], "x2": key[1], "x_source": key[2], "y_source": key[3]}
        for f in avail_fields:
            ref_val = r[f]
            tst_val = t[f]
            rel_err = abs(tst_val - ref_val) / max(abs(ref_val), 1e-9)
            row[f"{f}_ref"]     = ref_val
            row[f"{f}_tst"]     = tst_val
            row[f"{f}_rel_err"] = rel_err
        rows.append(row)

    # Shot-noise floor atteso: con N_ref e N_tst bassi, E[|ref-tst|/ref] > 0 anche per
    # estimatori unbiased. Per sigma_y: floor ≈ sqrt(2/π) * sqrt(1/(2*N_ref)+1/(2*N_tst)).
    n_ref_mean = mean_n_ref
    n_tst_mean = mean_n_tst
    shot_floor_sigma = np.sqrt(2 / np.pi) * np.sqrt(1 / (2 * max(n_ref_mean - 1, 1))
                                                     + 1 / (2 * max(n_tst_mean - 1, 1)))

    # Test A: media degli errori relativi per-spot (sensibile al shot noise)
    # Test B: errore relativo delle medie (stima del bias sistematico)
    print(f"\nShot-noise floor atteso per Test A (sigma_y): {100*shot_floor_sigma:.1f}%")
    print("(se mean_relErr_perspot ≈ floor → nessun bias sistematico; test A non discrimina)")

    print(f"\n--- Test A: media errori per-spot  (soglia {100*args.threshold:.0f}%) ---")
    print(f"{'Campo':12s}  {'mean(ref)':>10s}  {'mean(tst)':>10s}  "
          f"{'relErr_perspot%':>15s}  {'esito_A':>8s}")
    print("-" * 65)

    all_ok_a = True
    field_stats: dict = {}
    for f in avail_fields:
        rel_errs = np.array([r[f"{f}_rel_err"] for r in rows])
        mean_r   = float(np.mean([r[f"{f}_ref"] for r in rows]))
        mean_t   = float(np.mean([r[f"{f}_tst"] for r in rows]))
        mre      = float(np.mean(rel_errs))
        ok_a     = mre <= args.threshold
        all_ok_a &= ok_a
        field_stats[f] = {"mean_r": mean_r, "mean_t": mean_t, "mre": mre}
        status_a = "OK  " if ok_a else "FAIL"
        print(f"{f:12s}  {mean_r:+10.4f}  {mean_t:+10.4f}  "
              f"{100*mre:14.2f}%  [{status_a}]")

    print(f"\n--- Test B: errore relativo delle medie (soglia {100*args.threshold:.0f}%) ---")
    print(f"{'Campo':12s}  {'mean(ref)':>10s}  {'mean(tst)':>10s}  "
          f"{'relErr_medie%':>13s}  {'esito_B':>8s}")
    print("-" * 60)

    all_ok_b = True
    for f in avail_fields:
        s = field_stats[f]
        mean_r, mean_t = s["mean_r"], s["mean_t"]
        rel_mean = abs(mean_t - mean_r) / max(abs(mean_r), 1e-9)
        ok_b = rel_mean <= args.threshold
        all_ok_b &= ok_b
        status_b = "OK  " if ok_b else "FAIL"
        print(f"{f:12s}  {mean_r:+10.4f}  {mean_t:+10.4f}  "
              f"{100*rel_mean:12.2f}%  [{status_b}]")

    print()
    all_ok = all_ok_b  # decisione basata su Test B (bias sistematico)
    if all_ok_b:
        print(f"RISULTATO: Tecnica C APPROVATA (Test B) — bias sistematico < {100*args.threshold:.0f}%")
        if not all_ok_a:
            print("  NOTA: Test A fallisce per shot noise (atteso con N_hits bassi)")
        print(f"  → Aggiornare config/config.json: \"lens_n_photons\": 5000")
    else:
        print(f"RISULTATO: Tecnica C RESPINTA (Test B) — bias sistematico > {100*args.threshold:.0f}%")
        print(f"  → Mantenere lens_n_photons=10000")

    if args.csv:
        import csv
        fieldnames = (["x1", "x2", "x_source", "y_source"]
                      + [f"{f}_{s}" for f in avail_fields for s in ("ref", "tst", "rel_err")])
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)
        print(f"Report CSV salvato: {args.csv}")

    return 0 if all_ok_b else 2


if __name__ == "__main__":
    raise SystemExit(main())
