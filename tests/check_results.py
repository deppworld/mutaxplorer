#!/usr/bin/env python3
"""Compare MutaXplorer outputs with the simulated truth (tests/simulate_test_data.py).
  --somatic S1,S2      every truth somatic variant PASS; no germline/common SNP PASS; no other PASS calls
  --tumor-only S2      tumour-only samples: the RARE germline variant (absent from gnomAD) is expected to be
                       called somatic (documented limitation); common germline variants must be removed
  --germline N         both germline truth variants PASS; no somatic variant
  --cnv SAMPLE:GENE:CALL,...      gene-level CNV calls expected in <outdir>/somatic/cohort/cnv_gene_calls.tsv
  --fusion PATIENT:FUSION:EVIDENCE,...   rows expected in fusion_summary.tsv (gene order ignored)
  --germline-sv SAMPLE:SVTYPE:CHROM:POS  a germline SV within 200 bp in germline/cohort/sv_calls.tsv
  --identity A:B:STATUS,...       status prefix expected in qc sample_identity.tsv"""
import argparse, csv, glob, gzip, os, sys

ap = argparse.ArgumentParser()
ap.add_argument("truth"); ap.add_argument("outdir")
for o in ("somatic", "tumor-only", "germline", "cnv", "fusion", "germline-sv", "identity"):
    ap.add_argument(f"--{o}", default="")
a = ap.parse_args()
truth = {}
for r in csv.DictReader(open(a.truth), delimiter="\t"):
    truth.setdefault(r["class"], set()).add((r["chrom"], int(r["pos"]), r["ref"], r["alt"]))
som, germ, chet = truth.get("somatic", set()), truth.get("germline", set()), truth.get("common_het", set())
ok = True


def report(name, cond, msg):
    global ok
    ok &= bool(cond)
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}: {msg}")


def pass_calls(v):
    out = set()
    with gzip.open(v, "rt") as fh:
        for l in fh:
            if l[0] != "#":
                f = l.split("\t")
                if f[6] == "PASS":
                    out.add((f[0], int(f[1]), f[3], f[4]))
    return out


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t")) if os.path.exists(path) else None


tonly = set(filter(None, a.tumor_only.split(",")))
for sid in filter(None, a.somatic.split(",")):
    v = f"{a.outdir}/somatic/vcf/{sid}.norm.vcf.gz"
    if not os.path.exists(v):
        report(sid, False, "somatic VCF missing"); continue
    calls = pass_calls(v)
    report(sid, som <= calls, f"somatic sensitivity {len(som & calls)}/{len(som)}")
    rare_leak = calls & germ
    if sid in tonly:
        # Mutect2 may pass a common germline SNP in a tumour LOH region (VAF far from 50%); the cohort
        # population-AF filter must keep it out of the REPORTABLE calls
        rep = {(r["chrom"], int(r["pos"]), r["ref"], r["alt"]) for r in csv.DictReader(
            open(f"{a.outdir}/somatic/cohort/variants_all_pass.tsv"), delimiter="\t")
            if r["sample_id"] == sid and r["filter_reason"] == "REPORTABLE"}
        print(f"  [INFO] {sid}: common germline SNPs passing Mutect2 (LOH regions): {len(calls & chet)}")
        report(sid, not (rep & (chet | {g for g in germ if g in chet})), f"no common germline SNP reportable ({len(rep & chet)})")
        print(f"  [INFO] {sid}: tumour-only - germline variants called somatic (expected for the rare one): {len(rare_leak)}")
        fp = calls - som - germ - chet
    else:
        report(sid, not (calls & (germ | chet)), f"no germline variant called somatic ({len(calls & (germ | chet))})")
        fp = calls - som
    report(sid, not fp, f"other false-positive PASS calls: {len(fp)} {sorted(fp)[:3]}")
for sid in filter(None, a.germline.split(",")):
    v = f"{a.outdir}/germline/vcf/{sid}.norm.vcf.gz"
    if not os.path.exists(v):
        report(sid, False, "germline VCF missing"); continue
    calls = pass_calls(v)
    report(sid, germ <= calls, f"germline sensitivity {len(germ & calls)}/{len(germ)}")
    report(sid, len(chet & calls) >= 0.9 * len(chet), f"common heterozygous SNPs {len(chet & calls)}/{len(chet)}")
    report(sid, not (som & calls), f"no somatic variants in normal ({len(som & calls)})")
if a.cnv:
    rows = tsv(f"{a.outdir}/somatic/cohort/cnv_gene_calls.tsv") or []
    got = {(r["sample_id"], r["gene"]): r["call"] for r in rows}
    for item in a.cnv.split(","):
        s, g, c = item.split(":")
        report(f"{s} CNV", got.get((s, g)) == c, f"{g} expected {c}, got {got.get((s, g), 'none')}")
if a.fusion:
    rows = tsv(f"{a.outdir}/somatic/cohort/fusion_summary.tsv") or []
    got = {(r["patient"], frozenset(r["fusion"].split("--"))): r["evidence"] for r in rows}
    for item in a.fusion.split(","):
        p, fu, ev = item.split(":")
        g = got.get((p, frozenset(fu.split("--"))))
        report(f"{p} fusion", g == ev, f"{fu} expected '{ev}', got '{g}'")
if a.germline_sv:
    rows = tsv(f"{a.outdir}/germline/cohort/sv_calls.tsv") or []
    s, t, c, p = a.germline_sv.split(":")
    hit = [r for r in rows if r["sample_id"] == s and r["svtype"] == t and r["chrom1"] == c and abs(int(r["pos1"]) - int(p)) <= 200]
    report(f"{s} SV", hit, f"{t} at {c}:{p} {'found' if hit else 'missing'} ({len(rows)} PASS SVs)")
if a.identity:
    rows = tsv(f"{a.outdir}/somatic/cohort/sample_identity.tsv") or []
    for item in a.identity.split(","):
        x, y, st = item.split(":")
        r = [r for r in rows if {r["sample_a"], r["sample_b"]} == {x, y}]
        report("identity", r and r[0]["status"].startswith(st), f"{x} vs {y}: expected {st}, got {r[0]['status'] if r else 'none'}")
print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
sys.exit(0 if ok else 1)
