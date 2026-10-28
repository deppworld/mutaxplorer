#!/usr/bin/env python3
"""
Build a SYNTHETIC 120-patient targeted-panel cohort in the exact per-sample formats that MutaXplorer's
callers write (VEP-annotated Mutect2 VCF, mosdepth, MSIsensor2, CNVkit, Manta, Arriba), on the small
simulated genome from tests/simulate_test_data.py. Running the pipeline's own report scripts on it
produces the example figures in docs/example_report/. NOT real patient data.

usage: python3 tests/make_example_cohort.py --data tests/data --out example_calls [--n 120 --seed 7]
       then see docs/example_report/README.md for the three report commands.
"""
import argparse, gzip, math, os, random, re
import numpy as np

CTX_BASES = "ACGT"
COMP = str.maketrans("ACGT", "TGCA")
AA = ["Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met",
      "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"]
PROT_LEN = dict(KRAS=189, TP53=393, PIK3CA=1068, EML4=981, PTEN=403, CDKN2A=156, APC=2843, ERBB2=1255,
                MYC=439, BRCA2=3418, NF1=2839, CDH1=882, ALK=1620, EGFR=1210, BRAF=766, SMAD4=552,
                ATM=3056, RB1=928)
# hotspots: gene, exon index, offset in exon, protein change, consequence, cohort frequency
HOT = [("KRAS", 1, 40, "p.Gly12Asp", "missense_variant", .14), ("KRAS", 1, 41, "p.Gly12Val", "missense_variant", .08),
       ("KRAS", 1, 43, "p.Gly13Asp", "missense_variant", .04), ("KRAS", 2, 60, "p.Gln61His", "missense_variant", .03),
       ("TP53", 4, 50, "p.Arg175His", "missense_variant", .09), ("TP53", 6, 80, "p.Arg248Gln", "missense_variant", .08),
       ("TP53", 6, 120, "p.Arg273His", "missense_variant", .07), ("TP53", 7, 30, "p.Arg306Ter", "stop_gained", .04),
       ("PIK3CA", 3, 70, "p.Glu545Lys", "missense_variant", .07), ("PIK3CA", 7, 90, "p.His1047Arg", "missense_variant", .08),
       ("BRAF", 5, 100, "p.Val600Glu", "missense_variant", .09), ("EGFR", 4, 60, "p.Leu858Arg", "missense_variant", .05),
       ("APC", 6, 40, "p.Arg1450Ter", "stop_gained", .08), ("APC", 6, 140, "p.Glu1309fs", "frameshift_variant", .06),
       ("PTEN", 3, 50, "p.Arg130Ter", "stop_gained", .03), ("SMAD4", 5, 70, "p.Arg361Cys", "missense_variant", .03),
       ("ERBB2", 4, 110, "p.Val777Leu", "missense_variant", .02)]
BG_GENES = {"TP53": 6, "APC": 5, "ATM": 4, "NF1": 4, "BRCA2": 4, "PIK3CA": 2, "RB1": 3, "CDH1": 3, "SMAD4": 2,
            "PTEN": 3, "EGFR": 2, "ALK": 2, "MYC": 1, "ERBB2": 2, "BRAF": 1, "EML4": 1, "CDKN2A": 1, "KRAS": 0.5}
CNV_EVENTS = [("ERBB2", 10, .10), ("MYC", 7, .08), ("EGFR", 9, .06), ("CDKN2A", 0, .12), ("PTEN", 1, .10),
              ("RB1", 1, .07), ("SMAD4", 0, .05), ("KRAS", 6, .04), ("CDH1", 1, .05), ("ATM", 1, .04)]
FUSIONS = [("EML4", "ALK", .10, True), ("NF1", "CDH1", .02, False), ("BRAF", "SMAD4", .02, False)]

CSQ_FMT = "Allele|Consequence|IMPACT|SYMBOL|Feature|HGVSc|HGVSp|Existing_variation|CANONICAL|gnomADe_AF|CLIN_SIG"
IMPACT = dict(missense_variant="MODERATE", stop_gained="HIGH", frameshift_variant="HIGH", inframe_deletion="MODERATE",
              splice_donor_variant="HIGH", synonymous_variant="LOW")


def read_fasta(p):
    seq, name = {}, None
    for l in open(p):
        if l.startswith(">"):
            name = l[1:].split()[0]; seq[name] = []
        else:
            seq[name].append(l.strip().upper())
    return {k: "".join(v) for k, v in seq.items()}


def signatures():
    """Four illustrative SBS96 profiles (synthetic, shaped like clock-like / tobacco / APOBEC / MMR)."""
    ctx = [f"{a}[{m}]{b}" for m in ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G") for a in "ACGT" for b in "ACGT"]
    S = {}
    def mk(f):
        v = np.array([f(c) for c in ctx], float) + 1e-3; return v / v.sum()
    S["Clock-like"] = mk(lambda c: 12 if c[2:5] == "C>T" and c[-1] == "G" else 1.5 if c[2:5] in ("C>T", "T>C") else .3)
    S["Tobacco-like"] = mk(lambda c: 5 if c[2:5] == "C>A" else 1 if c[2:5] in ("C>T", "T>A") else .4)
    S["APOBEC-like"] = mk(lambda c: 9 if c[0] == "T" and c[2:5] in ("C>T", "C>G") and c[-1] in "AT" else .25)
    S["MMR-like"] = mk(lambda c: 3 if c[2:5] in ("C>T", "T>C") and c[-1] in "AG" else .6)
    return ctx, S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="tests/data from simulate_test_data.py")
    ap.add_argument("--out", required=True); ap.add_argument("--n", type=int, default=120)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    rnd = random.Random(a.seed); nrnd = np.random.default_rng(a.seed)
    os.makedirs(a.out, exist_ok=True)
    fa = read_fasta(os.path.join(a.data, "genome.fa"))
    genes = [l.split("\t")[:4] for l in open(os.path.join(a.data, "genes.bed")) if not l.startswith("#")]
    genes = {g: (c, int(s), int(e)) for c, s, e, g in genes}
    targets = [tuple(l.split("\t")[:3]) for l in open(os.path.join(a.data, "targets.bed"))]
    targets = [(c, int(s), int(e)) for c, s, e in targets]
    exons = {g: sorted((s, e) for c, s, e in targets if c == genes[g][0] and genes[g][1] <= s < genes[g][2]) for g in genes}

    # trinucleotide index of coding positions -> sample mutations by signature
    ctx, S = signatures()
    by_ctx = {c: [] for c in ctx}
    for g, ex in exons.items():
        chrom = genes[g][0]
        for s, e in ex:
            for p in range(s + 2, e - 1):           # 0-based
                tri = fa[chrom][p - 1:p + 2]
                if len(tri) < 3 or "N" in tri:
                    continue
                ref = tri[1]
                key_tri = tri if ref in "CT" else tri.translate(COMP)[::-1]
                for m in (("C>A", "C>G", "C>T") if key_tri[1] == "C" else ("T>A", "T>C", "T>G")):
                    by_ctx[f"{key_tri[0]}[{m}]{key_tri[2]}"].append((g, chrom, p + 1, ref, m))
    with open(os.path.join(a.out, "example_signatures.tsv"), "w") as f:
        f.write("Type\t" + "\t".join(S) + "\n")
        for i, c in enumerate(ctx):
            f.write(c + "\t" + "\t".join(f"{S[k][i]:.6f}" for k in S) + "\n")

    def alt_of(ref, m):
        alt = m[2]
        return alt if ref in "CT" else alt.translate(COMP)

    meta, chrom_len = [], {c: len(s) for c, s in fa.items()}
    for i in range(1, a.n + 1):
        pid = f"PT{i:03d}"; sid = f"{pid}_T"
        meta.append((sid, pid, "tumor", "dna"))
        msi_high = rnd.random() < .06
        purity = rnd.uniform(.3, .9)
        smoker = rnd.random() < .2; apobec = rnd.random() < .15
        recs = []
        # hotspots
        for g, ei, off, hp, csq, freq in HOT:
            if rnd.random() < freq:
                s, _ = exons[g][ei]; pos = s + off + 1; chrom = genes[g][0]
                ref = fa[chrom][pos - 1]
                if csq == "frameshift_variant":
                    recs.append((chrom, pos, ref + fa[chrom][pos:pos + 1], ref, g, csq, hp, True))
                else:
                    alt = {"A": "G", "C": "T", "G": "A", "T": "C"}[ref]
                    recs.append((chrom, pos, ref, alt, g, csq, hp, True))
        # background SNVs from signatures
        n_bg = int(nrnd.poisson(110 if msi_high else 16 if smoker else 1.5))
        w = (S["MMR-like"] * .8 + S["Clock-like"] * .2 if msi_high else
             S["Tobacco-like"] * .6 + S["Clock-like"] * .4 if smoker else
             S["APOBEC-like"] * .6 + S["Clock-like"] * .4 if apobec else S["Clock-like"])
        gw = np.array([BG_GENES[g] for g in genes], float); gl = list(genes)
        for _ in range(n_bg):
            for _try in range(20):
                c = ctx[nrnd.choice(96, p=w / w.sum())]
                g = gl[nrnd.choice(len(gl), p=gw / gw.sum())]
                cand = [x for x in by_ctx[c] if x[0] == g]
                if cand:
                    break
            if not cand:
                continue
            g, chrom, pos, ref, m = rnd.choice(cand)
            r = rnd.random()
            csq = "stop_gained" if r < .1 else "synonymous_variant" if r < .3 else "splice_donor_variant" if r < .33 else "missense_variant"
            pp = rnd.randint(1, PROT_LEN[g]); a1 = rnd.choice(AA)
            hp = f"p.{a1}{pp}{'Ter' if csq == 'stop_gained' else a1 if csq == 'synonymous_variant' else rnd.choice(AA)}"
            if csq == "splice_donor_variant":
                hp = ""
            recs.append((chrom, pos, ref, alt_of(ref, m), g, csq, hp, False))
        # indels (more in MSI-high)
        for _ in range(int(nrnd.poisson(12 if msi_high else .6))):
            g = rnd.choice(gl); s, e = rnd.choice(exons[g]); pos = rnd.randint(s + 5, e - 10); chrom = genes[g][0]
            L = rnd.choice([1, 1, 1, 2, 3]); ref = fa[chrom][pos - 1:pos + L]
            csq = "inframe_deletion" if L == 3 else "frameshift_variant"
            recs.append((chrom, pos, ref, ref[0], g, csq, f"p.{rnd.choice(AA)}{rnd.randint(1, PROT_LEN[g])}"
                         + ("del" if L == 3 else "fs"), False))
        # filtered-out examples: common SNPs, low VAF, low depth
        noise = []
        for _ in range(rnd.randint(2, 5)):
            g = rnd.choice(gl); s, e = rnd.choice(exons[g]); pos = rnd.randint(s + 5, e - 5); chrom = genes[g][0]
            ref = fa[chrom][pos - 1]; alt = {"A": "G", "C": "T", "G": "A", "T": "C"}[ref]
            noise.append((chrom, pos, ref, alt, g, rnd.choice(["missense_variant", "synonymous_variant"]),
                          f"p.{rnd.choice(AA)}{rnd.randint(1, PROT_LEN[g])}{rnd.choice(AA)}", rnd.choice(["snp", "lowvaf", "lowdp"])))
        write_vcf(a.out, sid, recs, noise, purity, rnd, chrom_len)
        write_cov(a.out, sid, targets, rnd, nrnd)
        with open(os.path.join(a.out, f"{sid}.msi"), "w") as f:
            tot = rnd.randint(80, 140); pct = rnd.uniform(28, 55) if msi_high else max(0, nrnd.normal(3, 2.5))
            f.write("Total_Number_of_Sites\tNumber_of_Somatic_Sites\t%\n"
                    f"{tot}\t{int(round(tot * pct / 100))}\t{pct:.2f}\n")
        cn_events = {g: cn for g, cn, fr in CNV_EVENTS if rnd.random() < fr}
        write_cnv(a.out, sid, genes, targets, chrom_len, cn_events, purity, rnd, nrnd)
        fus = [(g1, g2, k) for g1, g2, fr, k in FUSIONS if rnd.random() < fr]
        write_manta(a.out, sid, genes, fus, chrom_len, rnd)
        if rnd.random() < .45:                       # RNA available for ~45 % of patients
            rid = f"{pid}_RNA"; meta.append((rid, pid, "tumor", "rna"))
            write_arriba(a.out, rid, genes, [(g1, g2) for g1, g2, k in fus if rnd.random() < .9], rnd)
    with open(os.path.join(a.out, "samples.tsv"), "w") as f:
        for r in meta:
            f.write("\t".join(r) + "\n")
    with open(os.path.join(a.out, "known_fusions.tsv"), "w") as f:
        f.write("EML4\tALK\nBCR\tABL1\n")
    print(f"wrote {a.n} synthetic tumours to {a.out}")


def write_vcf(out, sid, recs, noise, purity, rnd, chrom_len):
    hdr = ["##fileformat=VCFv4.2", f"##tumor_sample={sid}", "##source=MutaXplorer synthetic example (not real data)"]
    hdr += [f"##contig=<ID={c},length={n}>" for c, n in chrom_len.items()]
    hdr += ['##FILTER=<ID=PASS,Description="All filters passed">',
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
            f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: {CSQ_FMT}">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
            '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sid]
    rows, seen = [], set()
    for chrom, pos, ref, alt, g, csq, hp, hot in recs:
        if (chrom, pos) in seen:
            continue
        seen.add((chrom, pos))
        dp = rnd.randint(250, 1100); vaf = min(.95, max(.04, rnd.gauss(purity / 2 * (1.0 if hot else .8), .06)))
        rows.append((chrom, pos, ref, alt, g, csq, hp, dp, vaf, "PASS", "", "COSV" + str(abs(hash((g, hp))) % 10**8) if hot else "",
                     "pathogenic" if hot else ""))
    for chrom, pos, ref, alt, g, csq, hp, kind in noise:
        if (chrom, pos) in seen:
            continue
        seen.add((chrom, pos))
        if kind == "snp":
            rows.append((chrom, pos, ref, alt, g, csq, hp, rnd.randint(300, 900), rnd.uniform(.42, .55), "PASS",
                         f"{rnd.uniform(.05, .45):.3f}", "rs" + str(rnd.randint(10**5, 10**8)), "benign"))
        elif kind == "lowvaf":
            rows.append((chrom, pos, ref, alt, g, csq, hp, rnd.randint(300, 900), rnd.uniform(.01, .035), "weak_evidence", "", "", ""))
        else:
            rows.append((chrom, pos, ref, alt, g, csq, hp, rnd.randint(30, 90), rnd.uniform(.08, .3), "PASS", "", "", ""))
    rows.sort(key=lambda r: (r[0], r[1]))
    lines = []
    for chrom, pos, ref, alt, g, csq, hp, dp, vaf, filt, gaf, ev, clin in rows:
        alt_n = max(1, int(round(dp * vaf))); hgvsc = f"c.{pos % 10000}{ref[0]}>{alt[0]}" if len(ref) == len(alt) else f"c.{pos % 10000}del{ref[1:]}"
        csqs = "|".join([alt, csq, IMPACT[csq], g, f"ENST_{g}", f"ENST_{g}:{hgvsc}", f"ENSP_{g}:{hp}" if hp else "",
                         ev, "YES", gaf, clin])
        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{filt}\tDP={dp};CSQ={csqs}\tGT:AD:AF:DP\t"
                     f"0/1:{dp - alt_n},{alt_n}:{vaf:.3f}:{dp}")
    body = "\n".join(hdr + lines) + "\n"
    for suf in (".annotated.vcf.gz", ".norm.vcf.gz"):
        with gzip.open(os.path.join(out, sid + suf), "wt") as f:
            f.write(body)


def write_cov(out, sid, targets, rnd, nrnd):
    mean = max(150, nrnd.normal(780, 220))
    size = sum(e - s for c, s, e in targets)
    with open(os.path.join(out, f"{sid}.mosdepth.summary.txt"), "w") as f:
        f.write("chrom\tlength\tbases\tmean\tmin\tmax\n"
                f"total_region\t{size}\t{int(size * mean)}\t{mean:.2f}\t{int(mean * .2)}\t{int(mean * 1.8)}\n")
    with gzip.open(os.path.join(out, f"{sid}.thresholds.bed.gz"), "wt") as f:
        f.write("#chrom\tstart\tend\tregion\t100X\t250X\t500X\n")
        for c, s, e in targets:
            L = e - s; d = mean * rnd.uniform(.5, 1.4)
            f.write(f"{c}\t{s}\t{e}\tunknown\t{L if d > 100 else int(L * .6)}\t{L if d > 250 else int(L * .4)}"
                    f"\t{L if d > 500 else int(L * .2)}\n")


def write_cnv(out, sid, genes, targets, chrom_len, events, purity, rnd, nrnd):
    def l2(cn):
        return math.log2(max(.02, (purity * cn + 2 * (1 - purity)) / 2))
    arm = {}                                          # occasional broad gain / loss of a chromosome half
    for c in chrom_len:
        if rnd.random() < .15:
            arm[c] = (rnd.choice([0, 1]), rnd.choice([1, 3]))
    segs = []
    for c, n in chrom_len.items():
        tg = [(s, e) for cc, s, e in targets if cc == c]
        lo, hi = tg[0][0], tg[-1][1]
        bps = {lo, hi, (lo + hi) // 2}
        ev = [(genes[g][1] - 25, genes[g][2] + 25, cn) for g, cn in events.items() if genes[g][0] == c]
        for s, e, _ in ev:
            bps |= {s, e}
        bps = sorted(bps)
        for s, e in zip(bps[:-1], bps[1:]):
            mid = (s + e) / 2; cn = 2
            if c in arm:
                half, acn = arm[c]
                if (mid < (lo + hi) / 2) == (half == 0):
                    cn = acn
            for es, ee, ecn in ev:
                if es <= mid < ee:
                    cn = ecn
            probes = [(ts, te) for ts, te in tg if s <= ts < e]
            if not probes:
                continue
            segs.append((c, s, e, cn, l2(cn) + nrnd.normal(0, .03), probes))
    with open(os.path.join(out, f"{sid}.call.cns"), "w") as f, open(os.path.join(out, f"{sid}.cnr"), "w") as g:
        f.write("chromosome\tstart\tend\tgene\tlog2\tci_hi\tci_lo\tcn\tdepth\tp_ttest\tprobes\tweight\n")
        g.write("chromosome\tstart\tend\tgene\tdepth\tlog2\tweight\n")
        for c, s, e, cn, lr, probes in segs:
            f.write(f"{c}\t{s}\t{e}\t-\t{lr:.4f}\t{lr + .02:.4f}\t{lr - .02:.4f}\t{cn}\t{600 * 2 ** lr:.1f}\t0.01\t{len(probes)}\t{len(probes)}\n")
            for ts, te in probes:
                r = lr + nrnd.normal(0, .12)
                g.write(f"{c}\t{ts}\t{te}\t-\t{600 * 2 ** r:.1f}\t{r:.4f}\t0.99\n")


def write_manta(out, sid, genes, fus, chrom_len, rnd):
    hdr = ["##fileformat=VCFv4.1", "##source=MutaXplorer synthetic example (not real data)"]
    hdr += [f"##contig=<ID={c},length={n}>" for c, n in chrom_len.items()]
    hdr += ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\t" + sid]
    L, k = [], 0
    def sup():
        return f"PR:SR\t{rnd.randint(300, 600)},0:{rnd.randint(500, 900)},0\t{rnd.randint(300, 600)},{rnd.randint(8, 70)}:{rnd.randint(500, 900)},{rnd.randint(10, 110)}"
    for g1, g2, _ in fus:
        c1, s1, e1 = genes[g1]; c2, s2, e2 = genes[g2]
        p1, p2 = rnd.randint(s1 + 500, e1 - 500), rnd.randint(s2 + 500, e2 - 500); k += 1
        if c1 != c2:
            L.append(f"{c1}\t{p1}\tMantaBND:{k}:0\tN\tN[{c2}:{p2}[\t.\tPASS\tSVTYPE=BND;MATEID=MantaBND:{k}:1;SOMATIC;SOMATICSCORE={rnd.randint(60, 600)}\t{sup()}")
        else:
            L.append(f"{c1}\t{min(p1, p2)}\tMantaINV:{k}\tN\t<INV>\t.\tPASS\tSVTYPE=INV;END={max(p1, p2)};SOMATIC;SOMATICSCORE={rnd.randint(60, 600)}\t{sup()}")
    for _ in range(rnd.choice([0, 0, 1, 1, 2, 3, 5])):
        c = rnd.choice(list(chrom_len)); p = rnd.randint(5000, chrom_len[c] - 20000); t = rnd.choice(["DEL", "DEL", "DUP", "INV", "INS"])
        k += 1; end = p + rnd.randint(300, 15000)
        L.append(f"{c}\t{p}\tManta{t}:{k}\tN\t<{t}>\t.\tPASS\tSVTYPE={t};END={end};SOMATIC;SOMATICSCORE={rnd.randint(30, 300)}\t{sup()}")
    with gzip.open(os.path.join(out, f"{sid}.manta.vcf.gz"), "wt") as f:
        f.write("\n".join(hdr + sorted(L)) + "\n")


def write_arriba(out, rid, genes, fus, rnd):
    with open(os.path.join(out, f"{rid}.arriba.tsv"), "w") as f:
        f.write("#gene1\tgene2\tstrand1(gene/fusion)\tstrand2(gene/fusion)\tbreakpoint1\tbreakpoint2\tsite1\tsite2\ttype\t"
                "split_reads1\tsplit_reads2\tdiscordant_mates\tcoverage1\tcoverage2\tconfidence\treading_frame\n")
        for g1, g2 in fus:
            c1, s1, _ = genes[g1]; c2, s2, _ = genes[g2]
            f.write(f"{g1}\t{g2}\t+/+\t+/+\t{c1}:{s1 + 1150}\t{c2}:{s2 + 1001}\tsplice-site\tsplice-site\t"
                    f"{'translocation' if c1 != c2 else 'inversion'}\t{rnd.randint(10, 300)}\t{rnd.randint(10, 300)}\t"
                    f"{rnd.randint(2, 60)}\t{rnd.randint(200, 1200)}\t{rnd.randint(200, 1200)}\thigh\tin-frame\n")


if __name__ == "__main__":
    main()
