#!/usr/bin/env python3
"""Decompose an eDNA detection deficit into sequential in silico filters.

Input: a pipe-separated table with one row per species recorded by a reference
method, with columns
  s       species name as reported by the study
  g       group (K = cryptobenthic, C = conspicuous)
  ed      1 if detected by eDNA, 0 otherwise
  lista   accepted name in the regional checklist (empty = outside checklist)
  filtro  original filter label (0_fuera_de_lista, 1_sin_referencia, 2_cebador,
          3_no_resoluble, 4_disponibilidad_u_otro)
  dmin    edit distance of the amplicon to the nearest non-conspecific
          reference (">6" allowed; empty = not computed)

Filters 0 to 2 (checklist, reference, primer) are taken from `filtro`.
Filters 3 and 4 are reassigned from `dmin` with the chosen threshold:
a species is resolvable when dmin >= threshold.

Usage:
  python decompose_deficit.py table.psv --thresholds 1 2 3 4 --out counts.csv \
      [--species-out per_species.csv --species-threshold 3]
"""
import argparse, csv, collections

ORDER = ["0_fuera_de_lista", "1_sin_referencia", "2_cebador", "3_no_resoluble", "4_sin_explicacion"]


def dval(x):
    x = (x or "").strip()
    if not x:
        return None
    return 99 if x.startswith(">") else float(x)


def classify(row, t):
    f = row["filtro"]
    if f[0] in "012":
        return f
    d = dval(row.get("dmin"))
    return "3_no_resoluble" if (d is not None and d < t) else "4_sin_explicacion"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--thresholds", type=float, nargs="+", default=[3])
    ap.add_argument("--out", required=True)
    ap.add_argument("--species-out")
    ap.add_argument("--species-threshold", type=float, default=3)
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.table, encoding="utf-8"), delimiter="|"))
    missed = [r for r in rows if r["ed"] == "0"]
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["umbral_dmin", "grupo", "n_no_detectadas"] + ORDER)
        for t in a.thresholds:
            for g in sorted({r["g"] for r in missed}):
                c = collections.Counter(classify(r, t) for r in missed if r["g"] == g)
                w.writerow([t, g, sum(c.values())] + [c.get(k, 0) for k in ORDER])
    if a.species_out:
        with open(a.species_out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["especie_estudio", "nombre_lista", "grupo", "detectada_eDNA", "dmin",
                        "vecino", "filtro_umbral_%g" % a.species_threshold])
            for r in rows:
                w.writerow([r["s"], r.get("lista", ""), r["g"], r["ed"], r.get("dmin", ""),
                            r.get("vec", ""), classify(r, a.species_threshold)])


if __name__ == "__main__":
    main()
