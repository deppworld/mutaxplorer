#!/usr/bin/env python3
"""MutaXplorer - gene BED (chrom, start0, end, gene_name, strand, biotype) from a GENCODE/Ensembl GTF.
Contig names are matched to the reference .fai (adds or strips 'chr'; MT <-> chrM)."""
import argparse, gzip, re, sys

ap = argparse.ArgumentParser()
ap.add_argument("--gtf", required=True); ap.add_argument("--fai", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
ref = [l.split("\t")[0] for l in open(a.fai)]
refset = set(ref)
order = {c: i for i, c in enumerate(ref)}


def fix(c):
    if c in refset:
        return c
    for alt in (("chr" + c) if not c.startswith("chr") else c[3:], "chrM" if c == "MT" else None, "MT" if c == "chrM" else None):
        if alt and alt in refset:
            return alt
    return None


attr = re.compile(r'(\S+) "([^"]*)"')
rows, skipped = [], 0
op = gzip.open if a.gtf.endswith(".gz") else open
with op(a.gtf, "rt") as fh:
    for l in fh:
        if l.startswith("#"):
            continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        c = fix(f[0])
        if c is None:
            skipped += 1; continue
        d = dict(attr.findall(f[8]))
        name = d.get("gene_name") or d.get("gene_id", "")
        bt = d.get("gene_type") or d.get("gene_biotype", "")
        rows.append((c, int(f[3]) - 1, int(f[4]), name, f[6], bt))
rows.sort(key=lambda r: (order.get(r[0], 10**9), r[1]))
with open(a.out, "w") as o:
    for r in rows:
        o.write("\t".join(map(str, r)) + "\n")
print(f"{len(rows)} genes written, {skipped} on contigs absent from the reference", file=sys.stderr)
