#!/usr/bin/env python3
"""
Analisis in silico de puntos ciegos del eDNA para peces criptobentonicos
del Golfo de California (proyecto CBF-2026-4124).

Pasos:
  1. Lista regional de especies (OBIS, poligono del Golfo de California).
  2. Cobertura de referencias en GenBank: 12S por region de marcador
     (MiFish-U, teleo, Riaz 12S-V5), COI y mitogenomas completos.
  3. Desajustes de cebadores por especie (secuencias que abarcan el sitio)
     y por familia (mitogenomas completos de GenBank, cualquier region).
  4. Resolucion: distancias entre amplicones de especies del Golfo.

Requisitos: Python >= 3.9, requests.  (pip install requests)
Uso:        python cripticos_insilico.py --email tu@correo --outdir resultados
            (opcional --api-key de NCBI para 10 consultas/s en lugar de 3)

El script reproduce el analisis preliminar del 21 de noviembre de 2026,
que se ejecuto en el navegador con la misma logica.
"""
import argparse, csv, json, os, sys, time
from collections import defaultdict
import requests

# ---------------------------------------------------------------- configuracion
WKT = ("POLYGON((-114.95 31.85,-114.4 31.95,-112.8 31.1,-111.0 28.6,-109.3 26.1,"
       "-108.0 24.6,-105.6 23.2,-109.4 22.85,-109.9 23.6,-110.4 24.3,-111.2 25.5,"
       "-112.0 26.8,-112.9 28.0,-113.7 29.3,-114.6 30.6,-114.95 31.85))")
CRYPTO = ["Apogonidae", "Blenniidae", "Bythitidae", "Chaenopsidae", "Dactyloscopidae",
          "Gobiesocidae", "Gobiidae", "Labrisomidae", "Tripterygiidae"]
CONSP = ["Pomacentridae", "Labridae", "Scaridae", "Serranidae", "Epinephelidae",
         "Lutjanidae", "Haemulidae", "Acanthuridae", "Chaetodontidae", "Pomacanthidae",
         "Balistidae", "Tetraodontidae", "Mullidae", "Kyphosidae", "Cirrhitidae"]
# especies de la tesis de El Khattabi Salazar (2014) en familias cripticas sin registros OBIS
EXTRA = {"Lythrypnus rhizophora": "Gobiidae", "Stathmonotus culebrai": "Chaenopsidae",
         "Stathmonotus lugubris": "Chaenopsidae", "Paraclinus tanygnathus": "Labrisomidae",
         "Axoclinus multicinctus": "Tripterygiidae"}
MIN_RECORDS = 3

PRIMERS = {  # 5'->3'
    "MiFishU_F": "GTCGGTAAAACTCGTGCCAGC",          # Miya et al. 2015
    "MiFishU_R": "CATAGTGGGGTATCTAATCCCAGTTTG",
    "teleo_F":   "ACACCGCCCGTCACTCT",              # Valentini et al. 2016
    "teleo_R":   "CTTCCGGTACACTTACCATG",
    "Riaz_F":    "ACTGGGATTAGATACCCC",             # Riaz et al. 2011 (12S-V5)
    "Riaz_R":    "TAGAACAGGCTCCTCTAG",
}
PAIRS = {"MiFishU": ("MiFishU_F", "MiFishU_R"), "teleo": ("teleo_F", "teleo_R"),
         "Riaz": ("Riaz_F", "Riaz_R")}
REF_ACC = "NC_027594.1"   # Holacanthus passer, mitogenoma de referencia para coordenadas
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# ---------------------------------------------------------------- utilidades
COMP = str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")
def rc(s): return s.translate(COMP)[::-1]

def best_site(seq, q, maxmm=9):
    best = (99, -1, [])
    L = len(q)
    for i in range(len(seq) - L + 1):
        mm, pos = 0, []
        for j in range(L):
            if seq[i + j] != q[j]:
                mm += 1; pos.append(j)
                if mm > maxmm: break
        if mm < best[0]:
            best = (mm, i, pos)
    return best

def primer_mm(seq, name):
    """Desajustes del cebador contra la secuencia (ambas hebras).
    Devuelve (desajustes, desajustes en las 5 bases del extremo 3')."""
    p = PRIMERS[name]; rev = name.endswith("_R")
    q = rc(p) if rev else p
    out = (99, 99)
    for s in (seq, rc(seq)):
        mm, _, pos = best_site(s, q)
        d3 = [j if rev else len(q) - 1 - j for j in pos]
        cand = (mm, sum(d < 5 for d in d3))
        if cand < out: out = cand
    return out

def amplicon(seq, pair, maxmm=6):
    f, r = (PRIMERS[x] for x in PAIRS[pair])
    for s in (seq, rc(seq)):
        mf, i, _ = best_site(s, f); mr, k, _ = best_site(s, rc(r))
        if mf <= maxmm and mr <= maxmm and k > i:
            return s[i + len(f):k]
    return None

def levenshtein(a, b):
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]

class NCBI:
    def __init__(self, email, api_key=None):
        self.p = {"email": email, "tool": "cripticos_insilico"}
        if api_key: self.p["api_key"] = api_key
        self.delay = 0.11 if api_key else 0.35
    def _get(self, ep, **kw):
        for k in range(4):
            try:
                r = requests.get(EUTILS + ep, params={**self.p, **kw}, timeout=60)
                r.raise_for_status(); time.sleep(self.delay); return r
            except Exception:
                time.sleep(2 + 3 * k)
        raise RuntimeError(f"NCBI fallo: {ep} {kw}")
    def search(self, term, retmax=0):
        j = self._get("esearch.fcgi", db="nuccore", retmode="json", retmax=retmax, term=term).json()
        return int(j["esearchresult"]["count"]), j["esearchresult"]["idlist"]
    def fasta(self, ids):
        out = {}
        for i in range(0, len(ids), 50):
            t = self._get("efetch.fcgi", db="nuccore", rettype="fasta", retmode="text",
                          id=",".join(ids[i:i + 50])).text
            for rec in t.split(">")[1:]:
                h, *s = rec.split("\n"); out[h] = "".join(s).upper()
        return out

# ---------------------------------------------------------------- mapeo a referencia
class RefMap:
    """Ubica cada secuencia en coordenadas del 12S de referencia por k-meros unicos."""
    def __init__(self, ref, K=8, span=1300):
        self.K = K; self.idx = defaultdict(list)
        for i in range(span - K):
            self.idx[ref[i:i + K]].append(i)
        self.coords = {}
        for n, p in PRIMERS.items():
            q = rc(p) if n.endswith("_R") else p
            mm, i, _ = best_site(ref[:span], q); self.coords[n] = (i, i + len(q))
    def map(self, s):
        best = (0, None, None)
        for seq in (s, rc(s)):
            hits = [(i, self.idx[seq[i:i + self.K]][0]) for i in range(len(seq) - self.K + 1)
                    if len(self.idx.get(seq[i:i + self.K], [])) == 1]
            if not hits: continue
            votes = defaultdict(int)
            for i, p in hits: votes[round((p - i) / 15)] += 1
            off = max(votes, key=votes.get) * 15
            hs = [p for i, p in hits if abs((p - i) - off) <= 30]
            if len(hs) > best[0]: best = (len(hs), min(hs), max(hs) + self.K)
        return best
    def covers(self, m, pair):
        n, a, b = m
        if not n: return 0
        f, r = PAIRS[pair]
        s0, e0 = self.coords[f][1], self.coords[r][0]   # inserto entre cebadores
        if a <= s0 + 20 and b >= e0 - 20: return 2       # completo
        if a <= e0 and b >= s0 and (min(b, e0) - max(a, s0)) >= 0.5 * (e0 - s0): return 1  # parcial
        return 0

# ---------------------------------------------------------------- pasos
def checklist_from_csv(path):
    """Lista regional propia (columna Species2). Familia y nombre aceptado desde GBIF."""
    import csv as _csv
    sp = {}
    with open(path, encoding="utf-8-sig") as fh:
        names = [r["Species2"].strip() for r in _csv.DictReader(fh) if r.get("Species2", "").strip()]
    names = [n for n in names if len(n.split()) == 2 and "species" not in n]
    for n in names:
        j = requests.get("https://api.gbif.org/v1/species/match", params={"name": n}, timeout=60).json()
        fam = j.get("family", "")
        sp[n] = {"familia": fam, "nombre_gbif": j.get("species", ""),
                 "grupo": "cripto" if fam in CRYPTO else ("conspicuo_arrecifal" if fam in REEF_CONSP else "otros"),
                 "registros_obis": ""}
    return sp

REEF_CONSP = ["Pomacentridae", "Labridae", "Scaridae", "Serranidae", "Lutjanidae", "Haemulidae",
              "Acanthuridae", "Chaetodontidae", "Pomacanthidae", "Balistidae", "Tetraodontidae",
              "Mullidae", "Kyphosidae", "Cirrhitidae", "Muraenidae", "Scorpaenidae", "Holocentridae",
              "Diodontidae", "Ostraciidae", "Monacanthidae", "Zanclidae", "Aulostomidae",
              "Fistulariidae", "Antennariidae"]

def checklist():
    sp = {}
    for fam in CRYPTO + CONSP:
        r = requests.get("https://api.obis.org/v3/checklist",
                         params={"scientificname": fam, "geometry": WKT, "size": 500}, timeout=60).json()
        for x in r["results"]:
            if x.get("taxonRank") == "Species" and x.get("records", 0) >= MIN_RECORDS:
                sp.setdefault(x["scientificName"], {"familia": fam,
                    "grupo": "cripto" if fam in CRYPTO else "conspicuo", "registros_obis": x["records"]})
    for s, f in EXTRA.items():
        sp.setdefault(s, {"familia": f, "grupo": "cripto", "registros_obis": 0})
    return sp

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", required=True); ap.add_argument("--api-key")
    ap.add_argument("--outdir", default="resultados_insilico")
    ap.add_argument("--species-csv", help="Lista propia de especies (columna Species2); si se omite se usa OBIS")
    a = ap.parse_args(); os.makedirs(a.outdir, exist_ok=True)
    ncbi = NCBI(a.email, a.api_key)

    if a.species_csv:
        print("1. Lista regional (archivo propio + familias GBIF)"); sp = checklist_from_csv(a.species_csv)
    else:
        print("1. Lista regional (OBIS)"); sp = checklist()
    print(f"   {len(sp)} especies")
    ref = ncbi.fasta([REF_ACC]); ref = next(iter(ref.values()))
    rm = RefMap(ref)

    print("2. Referencias en GenBank")
    seqs = {}
    for s, d in sp.items():
        alt = d.get("nombre_gbif")
        o = f'("{s}"[Organism] OR "{alt}"[Organism])' if alt and alt != s else f'"{s}"[Organism]'
        n12, ids = ncbi.search(f'{o} AND (12S OR "small subunit ribosomal RNA" OR "complete genome" '
                               f'OR mitogenome OR "partial genome") NOT PREDICTED[Title]', 30)
        d["coi"], _ = ncbi.search(f'{o} AND (COI OR COX1 OR "cytochrome oxidase subunit I" OR '
                                  f'"cytochrome c oxidase subunit 1" OR "cytochrome oxidase subunit 1")')
        d["mitogenomas"], _ = ncbi.search(f'{o} AND mitochondrion[Title] AND "complete genome"[Title]')
        seqs[s] = {h: q for h, q in ncbi.fasta(ids).items() if len(q) < 100000} if ids else {}
        for pair in PAIRS:
            d[f"ref_{pair}"] = max([rm.covers(rm.map(q), pair) for q in seqs[s].values()] or [0])
        d["ref_12S_cualquiera"] = int(any(rm.map(q)[0] for q in seqs[s].values()))
        for name in PRIMERS:   # desajustes en secuencias que contienen el sitio
            vals = [primer_mm(q, name) for q in seqs[s].values()]
            vals = [v for v in vals if v[0] <= 6]
            d[f"mm_{name}"] = "" if not vals else "%d/%d" % min(vals)
        for pair in PAIRS:
            amps = [amplicon(q, pair) for q in seqs[s].values()]
            d[f"amp_{pair}"] = sorted({x for x in amps if x}, key=len)
        print(f"   {s}: 12S={d['ref_12S_cualquiera']} MiFish={d['ref_MiFishU']} teleo={d['ref_teleo']}")

    print("3. Resolucion")
    for pair in PAIRS:
        have = {s: d[f"amp_{pair}"] for s, d in sp.items() if d[f"amp_{pair}"]}
        names = list(have)
        for s in names:
            dmin, vec = 10**9, ""
            for t in names:
                if t == s: continue
                dd = min(levenshtein(x, y) for x in have[s] for y in have[t])
                if dd < dmin: dmin, vec = dd, t
            sp[s][f"dmin_{pair}"] = dmin; sp[s][f"vecino_{pair}"] = vec

    cols = ["especie", "familia", "grupo", "registros_obis", "ref_12S_cualquiera",
            "ref_MiFishU", "ref_teleo", "ref_Riaz", "coi", "mitogenomas"] + \
           [f"mm_{n}" for n in PRIMERS] + \
           [f"{k}_{p}" for p in PAIRS for k in ("dmin", "vecino")]
    with open(os.path.join(a.outdir, "especies_golfo_insilico.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(cols)
        for s, d in sorted(sp.items(), key=lambda x: (x[1]["grupo"], x[1]["familia"], x[0])):
            w.writerow([s] + [d.get(c, "") for c in cols[1:]])
    print("Listo:", a.outdir)

if __name__ == "__main__":
    main()
