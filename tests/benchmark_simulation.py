#!/usr/bin/env python3
"""
MutaXplorer simulation benchmark (panel-like design).

Genome: 3 x 100 kb chromosomes; 30 targets of 2 kb.
Planted per tumor:  60 somatic SNVs + 20 somatic indels (1-12 bp) at ONE VAF level
Planted germline :  20 rare het (absent from the population VCF) + 20 common het (in it, AF 0.05-0.5)
Samples: one normal (P), tumors at VAF 2/5/10/20 %, two replicates each (different read seeds).

usage:
  benchmark_simulation.py make <dir> [--depth 500]
  benchmark_simulation.py eval <dir> <results_paired> <results_tumoronly> <out_prefix>
"""
import gzip, os, random, subprocess, sys, csv

VAFS = [0.02, 0.05, 0.10, 0.20]
REPS = 2
comp = lambda b: {"A": "C", "C": "G", "G": "T", "T": "A"}[b]


def sh(c):
    subprocess.run(c, shell=True, check=True, executable="/bin/bash")


def make(d, depth=500):
    random.seed(2026)
    os.makedirs(f"{d}/fastq", exist_ok=True); os.makedirs(f"{d}/tmp", exist_ok=True)
    chroms = {f"chr{i}": 100000 for i in (1, 2, 3)}
    seq = {c: "".join(random.choice("ACGT") for _ in range(n)) for c, n in chroms.items()}
    regs = [(c, s, s + 2000) for c in chroms for s in range(5000, 95000, 9000)]       # 30 targets
    used = set()

    def pick(n, margin=60):
        out = []
        while len(out) < n:
            c, s, e = random.choice(regs); p = random.randint(s + margin, e - margin)
            if all(abs(p - q) > 40 for cc, q in used if cc == c):
                used.add((c, p)); out.append((c, p))
        return out
    snv = [("snv", c, p, seq[c][p - 1], comp(seq[c][p - 1])) for c, p in pick(60)]
    indel = []
    for c, p in pick(20):
        L = random.choice([1, 2, 3, 4, 6, 9, 12])
        if random.random() < 0.5:   # deletion: REF = anchor + L bases
            indel.append(("del", c, p, seq[c][p - 1:p + L], seq[c][p - 1]))
        else:                       # insertion
            ins = "".join(random.choice("ACGT") for _ in range(L))
            indel.append(("ins", c, p, seq[c][p - 1], seq[c][p - 1] + ins))
    rare = [("rare_germ", c, p, seq[c][p - 1], comp(seq[c][p - 1])) for c, p in pick(20)]
    common = [("common_germ", c, p, seq[c][p - 1], comp(seq[c][p - 1])) for c, p in pick(20)]

    def apply(vs):
        s = {c: list(v) for c, v in seq.items()}
        for t, c, p, r, a in sorted(vs, key=lambda x: -x[2]):          # right-to-left keeps coords valid
            s[c][p - 1:p - 1 + len(r)] = list(a)
        return {c: "".join(v) for c, v in s.items()}
    haps = {"A": apply([]), "B": apply(rare + common), "T": apply(rare + common + snv + indel)}
    with open(f"{d}/genome.fa", "w") as f:
        for c, v in seq.items():
            f.write(f">{c}\n{v}\n")
    for k, h in haps.items():
        with open(f"{d}/tmp/{k}.fa", "w") as f:
            for c, s, e in regs:
                f.write(f">{c}_{s}\n{h[c][s - 400:e + 400]}\n")
    with open(f"{d}/targets.bed", "w") as f:
        f.writelines(f"{c}\t{s}\t{e}\n" for c, s, e in regs)
    # population resource: the common germline variants + background sites
    hdr = ("##fileformat=VCFv4.2\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n"
           + "".join(f"##contig=<ID={c},length={n}>\n" for c, n in chroms.items())
           + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    rows = [(c, p, r, a, round(random.uniform(0.05, 0.5), 3)) for _, c, p, r, a in common]
    bg = [(c, p, seq[c][p - 1], comp(seq[c][p - 1]), round(random.uniform(1e-4, 0.3), 4)) for c, p in pick(150)]
    with open(f"{d}/gnomad.vcf", "w") as f:
        f.write(hdr)
        for c, p, r, a, af in sorted(rows + bg, key=lambda x: (x[0], x[1])):
            f.write(f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\tAF={af}\n")
    sh(f"bgzip -f {d}/gnomad.vcf && tabix -f -p vcf {d}/gnomad.vcf.gz")
    # truth, left-aligned exactly like the pipeline (bcftools norm) so indels compare correctly
    with open(f"{d}/truth.raw.vcf", "w") as f:
        f.write(hdr.replace("#CHROM", '##INFO=<ID=CLS,Number=1,Type=String,Description="class">\n#CHROM'))
        for t, c, p, r, a in sorted(snv + indel + rare + common, key=lambda x: (x[1], x[2])):
            f.write(f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\tCLS={t}\n")
    sh(f"bcftools norm -f {d}/genome.fa -o {d}/truth.vcf {d}/truth.raw.vcf 2>/dev/null")
    with open(f"{d}/truth.tsv", "w") as f, open(f"{d}/truth.vcf") as v:
        f.write("class\tchrom\tpos\tref\talt\n")
        for l in v:
            if l[0] != "#":
                x = l.split("\t"); f.write(f"{x[7].split('CLS=')[1].split(';')[0].strip()}\t{x[0]}\t{x[1]}\t{x[3]}\t{x[4]}\n")

    span = sum(e - s + 800 for _, s, e in regs)
    pairs_total = int(depth * span / 300)

    def sample(name, fr, seed):
        parts = []
        for i, (k, frac) in enumerate(fr.items()):
            n = int(pairs_total * frac)
            if n == 0:
                continue
            sh(f"wgsim -S {seed * 10 + i} -N {n} -1 150 -2 150 -d 300 -s 30 -e 0.002 -r 0 -R 0 -h "
               f"{d}/tmp/{k}.fa {d}/tmp/{name}{k}_1.fq {d}/tmp/{name}{k}_2.fq > /dev/null 2>&1")
            parts.append(k)
        for r in (1, 2):
            sh(f"cat " + " ".join(f"{d}/tmp/{name}{k}_{r}.fq" for k in parts) + f" | gzip -1 > {d}/fastq/{name}_R{r}.fastq.gz")
        sh(f"rm -f {d}/tmp/{name}*_?.fq")
    sample("N", {"A": 0.5, "B": 0.5}, 1)
    sheet_p = ["patient,sample,status,fastq_1,fastq_2,bam", f"P,N,normal,{d}/fastq/N_R1.fastq.gz,{d}/fastq/N_R2.fastq.gz,"]
    sheet_t = [sheet_p[0]]
    for v in VAFS:
        for rep in range(1, REPS + 1):
            nm = f"T{int(v * 100):02d}_r{rep}"
            sample(nm, {"A": 0.5, "B": 0.5 - v, "T": v}, int(v * 1000) + rep)
            row = f"{d}/fastq/{nm}_R1.fastq.gz,{d}/fastq/{nm}_R2.fastq.gz,"
            sheet_p.append(f"P,{nm},tumor,{row}")
            sheet_t.append(f"{nm},{nm},tumor,{row}")
    open(f"{d}/samplesheet_paired.csv", "w").write("\n".join(sheet_p) + "\n")
    open(f"{d}/samplesheet_tumoronly.csv", "w").write("\n".join(sheet_t) + "\n")
    sh(f"rm -rf {d}/tmp")
    print(f"benchmark data in {d}: {len(regs)} targets, ~{depth}x, {pairs_total} read pairs per sample")


def pass_calls(vcf):
    out = set()
    with gzip.open(vcf, "rt") as fh:
        for l in fh:
            if l[0] != "#":
                f = l.split("\t")
                if f[6] == "PASS":
                    out.add((f[0], int(f[1]), f[3], f[4]))
    return out


def evaluate(d, res_p, res_t, prefix):
    truth = {}
    for r in csv.DictReader(open(f"{d}/truth.tsv"), delimiter="\t"):
        truth[(r["chrom"], int(r["pos"]), r["ref"], r["alt"])] = r["class"]
    som = {k for k, c in truth.items() if c in ("snv", "del", "ins")}
    rare = {k for k, c in truth.items() if c == "rare_germ"}
    common = {k for k, c in truth.items() if c == "common_germ"}
    rows = []
    for mode, res in (("tumor_normal", res_p), ("tumor_only", res_t)):
        for v in VAFS:
            for rep in range(1, REPS + 1):
                nm = f"T{int(v * 100):02d}_r{rep}"
                calls = pass_calls(f"{res}/somatic/vcf/{nm}.norm.vcf.gz")
                tp_snv = sum(1 for k in calls & som if truth[k] == "snv")
                tp_ind = sum(1 for k in calls & som if truth[k] != "snv")
                fp_other = len(calls - som - rare - common)
                rows.append(dict(mode=mode, vaf=v, replicate=rep,
                                 snv_recall=round(tp_snv / 60, 4), indel_recall=round(tp_ind / 20, 4),
                                 somatic_recall=round((tp_snv + tp_ind) / 80, 4),
                                 rare_germline_called=len(calls & rare), common_germline_called=len(calls & common),
                                 other_false_positives=fp_other,
                                 precision=round((tp_snv + tp_ind) / max(1, len(calls)), 4)))
    with open(f"{prefix}.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    # summary figure
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    import pandas as pd
    df = pd.DataFrame(rows)
    g = df.groupby(["mode", "vaf"]).mean(numeric_only=True).reset_index()
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6), gridspec_kw=dict(wspace=0.3))
    for mode, col, mk in (("tumor_normal", "#2a78d6", "o"), ("tumor_only", "#eb6834", "s")):
        s = g[g["mode"] == mode]
        ax[0].plot(s["vaf"] * 100, s["snv_recall"] * 100, marker=mk, color=col, lw=2, label=f"{mode.replace('_', '-')} SNV")
        ax[0].plot(s["vaf"] * 100, s["indel_recall"] * 100, marker=mk, color=col, lw=2, ls="--", label=f"{mode.replace('_', '-')} indel")
        ax[1].bar([x + (-0.18 if mode == "tumor_normal" else 0.18) for x in range(len(s))],
                  s["rare_germline_called"], width=0.35, color=col, label=mode.replace("_", "-"))
    ax[0].set_xlabel("Variant allele frequency (%)"); ax[0].set_ylabel("Recall (%)"); ax[0].set_ylim(0, 105)
    ax[0].set_xticks([v * 100 for v in VAFS]); ax[0].set_title("Somatic recall vs VAF", loc="left", fontweight="bold")
    ax[0].legend(frameon=False, fontsize=7.5); ax[0].grid(color="#e1e0d9", lw=0.6)
    ax[1].set_xticks(range(len(VAFS))); ax[1].set_xticklabels([f"{int(v*100)}%" for v in VAFS])
    ax[1].set_ylabel("Rare germline variants called somatic (of 20)"); ax[1].set_xlabel("Tumor VAF level")
    ax[1].set_title("Germline leakage", loc="left", fontweight="bold"); ax[1].legend(frameon=False, fontsize=8)
    ax[1].grid(axis="y", color="#e1e0d9", lw=0.6)
    for a_ in ax:
        a_.spines[["top", "right"]].set_visible(False)
    fig.savefig(f"{prefix}.png", dpi=200, bbox_inches="tight")
    print(g.to_string())


if __name__ == "__main__":
    if sys.argv[1] == "make":
        make(os.path.abspath(sys.argv[2]), int(sys.argv[4]) if len(sys.argv) > 4 else 500)
    else:
        evaluate(*sys.argv[2:6])
