#!/usr/bin/env python3
"""Weekly monitor of public reference sequences for a species list (ENA/INSDC).

For each species: resolve the NCBI/ENA taxon id, count sequence records by marker
(12S, 16S, COI, complete mitogenome) and SRA/ENA read runs, and list accessions.
Compare with a previous snapshot and write a Markdown report of new records.

Usage:
  python ena_reference_monitor.py --species species.csv --out snapshot_YYYY-MM-DD.csv \
      [--previous snapshot_prev.csv] [--report report.md] [--name-col especie]
"""
import argparse, csv, json, os, re, ssl, sys, time, urllib.parse, urllib.request

ENA_TAX = "https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/"
ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"


def _ctx():
    ca = os.environ.get("SSL_CERT_FILE") or os.environ.get("REQUESTS_CA_BUNDLE")
    if not ca:
        cand = os.path.expanduser("~/.ccr/ca-bundle.crt")
        ca = cand if os.path.exists(cand) else None
    return ssl.create_default_context(cafile=ca) if ca else ssl.create_default_context()


CTX = _ctx()


def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, context=CTX, timeout=90) as r:
                return r.read().decode()
        except Exception as e:  # network hiccup or 5xx
            if i == tries - 1:
                raise
            time.sleep(3 * (i + 1))


def taxid(name):
    try:
        j = json.loads(get(ENA_TAX + urllib.parse.quote(name)) or "[]")
    except Exception:
        return ""
    for t in j:
        if t.get("scientificName", "").lower() == name.lower():
            return t["taxId"]
    return j[0]["taxId"] if j else ""


def search(result, query, fields):
    q = urllib.parse.urlencode({"result": result, "query": query, "fields": fields,
                                "limit": 0, "format": "tsv"})
    txt = get(f"{ENA_SEARCH}?{q}")
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    return rows


def classify(desc, length):
    d = desc.lower()
    out = set()
    if "complete genome" in d and "mitochondri" in d:
        out.add("mito")
    if re.search(r"12s|small subunit ribosomal rna", d) and "mitochondri" in d:
        out.add("12S")
    if re.search(r"16s|large subunit ribosomal rna", d) and "mitochondri" in d:
        out.add("16S")
    if re.search(r"\bcoi\b|\bcox1\b|\bco1\b|cytochrome (c )?oxidase subunit (1|i)\b", d):
        out.add("COI")
    if "mito" in out:
        out |= {"12S", "16S", "COI"}
    return out


def species_record(name):
    tid = taxid(name)
    rec = {"species": name, "taxid": tid, "n12S": 0, "n16S": 0, "nCOI": 0, "nmito": 0,
           "n_runs": 0, "acc12S": "", "acc16S": "", "accCOI": "", "accmito": "", "runs": ""}
    if not tid:
        return rec
    seqs = search("sequence", f"tax_eq({tid})", "accession,description,base_count,first_public")
    acc = {"12S": [], "16S": [], "COI": [], "mito": []}
    for s in seqs:
        for m in classify(s.get("description", ""), s.get("base_count", "0")):
            acc[m].append(s["accession"])
    for m in acc:
        rec["n" + m] = len(acc[m])
        rec["acc" + m] = ";".join(sorted(set(acc[m])))
    runs = search("read_run", f"tax_eq({tid})", "run_accession,library_strategy,first_public")
    rec["n_runs"] = len(runs)
    rec["runs"] = ";".join(sorted(f"{r['run_accession']}:{r.get('library_strategy','')}" for r in runs))
    return rec


def load(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--name-col", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--previous")
    ap.add_argument("--report")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    rows = load(a.species)
    col = a.name_col or next(c for c in rows[0] if c.lower() in ("especie", "species", "species2", "nombre"))
    names = [r[col].strip() for r in rows if r[col].strip()]
    if a.limit:
        names = names[: a.limit]

    fields = ["species", "taxid", "n12S", "n16S", "nCOI", "nmito", "n_runs",
              "acc12S", "acc16S", "accCOI", "accmito", "runs"]
    out = []
    for i, n in enumerate(names, 1):
        try:
            out.append(species_record(n))
        except Exception as e:
            print(f"[warn] {n}: {e}", file=sys.stderr)
            out.append({"species": n, "taxid": "ERROR"})
        if i % 50 == 0:
            print(f"{i}/{len(names)}", file=sys.stderr)
        time.sleep(0.2)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(out)

    if a.previous and a.report and os.path.exists(a.previous):
        prev = {r["species"]: r for r in load(a.previous)}
        lines = ["# New public records since previous snapshot", ""]
        n_new = 0
        for r in out:
            p = prev.get(r["species"])
            if not p or r.get("taxid") in ("", "ERROR"):
                continue
            for m in ("12S", "16S", "COI", "mito"):
                new = set(filter(None, r.get("acc" + m, "").split(";"))) - set(filter(None, p.get("acc" + m, "").split(";")))
                if new:
                    first = " (FIRST reference for this marker)" if not p.get("acc" + m) else ""
                    lines.append(f"- {r['species']}: {m} +{len(new)}{first}: {', '.join(sorted(new))}")
                    n_new += 1
            newruns = set(filter(None, r.get("runs", "").split(";"))) - set(filter(None, p.get("runs", "").split(";")))
            if newruns:
                lines.append(f"- {r['species']}: new read runs: {', '.join(sorted(newruns))}")
                n_new += 1
        if n_new == 0:
            lines.append("No new records.")
        with open(a.report, "w") as f:
            f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
