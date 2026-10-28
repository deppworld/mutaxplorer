#!/usr/bin/env python3
"""
MutaXplorer test data: a 600-kb, 3-chromosome genome with 18 genes and KNOWN events.

  SNV/indel : 6 somatic (5 SNVs + 9-bp deletion, VAF ~30%); germline G1 (common, in gnomAD),
              G2 (rare, absent from gnomAD) and 15 common heterozygous SNPs (identity markers)
  CNV       : ERBB2 amplification (tumour CN 8), CDKN2A homozygous deletion (CN 0), PTEN loss (CN 1);
              tumour purity 60%
  SV/fusion : t(1;3) joining EML4 intron 2 to ALK intron 1 -> EML4--ALK in tumour DNA and RNA
  germline SV: 2-kb heterozygous deletion chr3:115001-117000 in the WGS normal
Samples: P1_T (FASTQ tumour) + P1_N (FASTQ normal) + P1_RNA (RNA FASTQ); P2_T (imported BAM, no read
group, tumour-only); P3_T (UMI FASTQ, every molecule read 3x); WGS_N (whole-genome normal).
Needs wgsim, bwa-mem2, samtools, bgzip, tabix, bcftools.
"""
import gzip, os, random, subprocess, sys

out = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "tests/data")
for d in ("fastq", "bam", "vep_cache", "tmp"):
    os.makedirs(f"{out}/{d}", exist_ok=True)
random.seed(7)
CL = 200000
seq = {c: "".join(random.choice("ACGT") for _ in range(CL)) for c in ("chr1", "chr2", "chr3")}
GENES = {"chr1": [("KRAS", 10000), ("TP53", 40000), ("PIK3CA", 70000), ("EML4", 100000), ("PTEN", 140000), ("CDKN2A", 170000)],
         "chr2": [("APC", 10000), ("ERBB2", 40000), ("MYC", 70000), ("BRCA2", 100000), ("NF1", 130000), ("CDH1", 160000)],
         "chr3": [("ALK", 10000), ("EGFR", 40000), ("BRAF", 70000), ("SMAD4", 100000), ("ATM", 130000), ("RB1", 160000)]}
EXON = [(i * 1000, i * 1000 + 150) for i in range(8)]   # 8 exons of 150 bp per gene (0-based, half-open)
gpos = {g: (c, s) for c, gl in GENES.items() for g, s in gl}
exon = lambda g, i: (gpos[g][0], gpos[g][1] + EXON[i][0], gpos[g][1] + EXON[i][1])
comp = lambda b: {"A": "C", "C": "G", "G": "T", "T": "A"}[b]

# ---- variants ------------------------------------------------------------------
def snv(g, i, off):
    c, s, e = exon(g, i); p = s + off + 1
    return (c, p, seq[c][p - 1], comp(seq[c][p - 1]))
SOM = [snv("KRAS", 1, 75), snv("TP53", 0, 60), snv("PIK3CA", 5, 100), snv("APC", 3, 40), snv("BRAF", 6, 90)]
c, s, e = exon("EGFR", 0); DEL = (c, s + 50, seq[c][s + 49:s + 59], seq[c][s + 49])   # 9-bp deletion, VCF style
G1, G2 = snv("TP53", 4, 70), snv("BRCA2", 2, 50)
COMMON = []
for g in gpos:                       # 2 common heterozygous SNPs per gene (identity markers)
    COMMON += [snv(g, 7, 30), snv(g, 4, 110)]

def apply(base, vs):
    s = {k: list(v) for k, v in base.items()}
    for c, p, r, a in sorted(vs, key=lambda x: -x[1]):
        s[c][p - 1:p - 1 + len(r)] = list(a)
    return {k: "".join(v) for k, v in s.items()}
hapA = dict(seq)
hapB = apply(seq, [G1, G2] + COMMON)
hapT = apply(seq, [G1, G2] + COMMON + SOM + [DEL])        # tumour allele (somatic on the B haplotype)
hapB_wgs = {k: v for k, v in hapB.items()}
hapB_wgs["chr3"] = hapB["chr3"][:115000] + hapB["chr3"][117000:]   # germline 2-kb deletion

# ---- targets: exons + fusion baits ------------------------------------------------
BP1 = gpos["EML4"][1] + 1600          # in EML4 intron 2 (chr1)
BP2 = gpos["ALK"][1] + 600            # in ALK intron 1 (chr3)
targets = [exon(g, i) for c, gl in GENES.items() for g, _ in gl for i in range(len(EXON))]
targets += [("chr1", BP1 - 350, BP1 + 350), ("chr3", BP2 - 350, BP2 + 350)]
targets.sort()
CN = {}   # tumour copy number per target (default 2)
for t in targets:
    c, s, e = t
    if c == "chr2" and 35000 <= s < 55000: CN[t] = 8
    if c == "chr1" and 165000 <= s < 180000: CN[t] = 0
    if c == "chr1" and 135000 <= s < 150000: CN[t] = 1

def region(h, c, s, e, pad=400):
    return h[c][max(0, s - pad):e + pad]

def write_mix(path, entries):
    with open(path, "w") as f:
        for i, (name, sq) in enumerate(entries):
            f.write(f">{name}_{i}\n{sq}\n")

def dna_mix(kind):
    """units: normal allele copy = 2 units (40% normal), tumour copy = 3 units (60% purity)"""
    E = []
    for t in targets:
        c, s, e = t
        if kind == "normal":
            E += [("A", region(hapA, c, s, e))] * 5 + [("B", region(hapB, c, s, e))] * 5
            continue
        cn = CN.get(t, 2)
        E += [("A", region(hapA, c, s, e))] * 2 + [("B", region(hapB, c, s, e))] * 2
        tum = {2: [hapA, hapT], 8: [hapA] + [hapT] * 7, 1: [hapT], 0: []}[cn]
        for h in tum:
            E += [("T", region(h, c, s, e))] * 3
    if kind == "tumor":   # translocation junctions (both derivatives), one tumour copy each
        der1 = seq["chr1"][BP1 - 1400:BP1] + seq["chr3"][BP2:BP2 + 1400]
        der2 = seq["chr3"][BP2 - 1400:BP2] + seq["chr1"][BP1:BP1 + 1400]
        E += [("D1", der1)] * 3 + [("D2", der2)] * 3
    return E

def sh(cmd):
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")

def wgsim(fa, n, seed, tag, rl=150):
    sh(f"wgsim -S {seed} -N {n} -1 {rl} -2 {rl} -d 300 -s 30 -e 0.002 -r 0 -R 0 -h {fa} "
       f"{out}/tmp/{tag}_1.fq {out}/tmp/{tag}_2.fq > /dev/null 2>&1")
    return f"{out}/tmp/{tag}_1.fq", f"{out}/tmp/{tag}_2.fq"

def to_fastq(parts, name, umi_copies=0):
    o1, o2 = gzip.open(f"{out}/fastq/{name}_R1.fastq.gz", "wt"), gzip.open(f"{out}/fastq/{name}_R2.fastq.gz", "wt")
    k = 0
    for f1, f2 in parts:
        with open(f1) as a, open(f2) as b:
            while True:
                r1 = [a.readline() for _ in range(4)]; r2 = [b.readline() for _ in range(4)]
                if not r1[0]:
                    break
                k += 1; umi = "".join(random.choice("ACGT") for _ in range(8))
                for cp in range(umi_copies or 1):
                    nm = f"@R{k}x{cp}" + (f":{umi}" if umi_copies else "")
                    for rec, o in ((r1, o1), (r2, o2)):
                        sq = list(rec[1].strip())
                        if umi_copies:
                            i = random.randrange(len(sq)); sq[i] = random.choice("ACGT")
                        o.write(f"{nm}\n{''.join(sq)}\n+\n{rec[3]}")
    o1.close(); o2.close()

def simulate_dna(kind, name, depth, seed, umi_copies=0):
    E = dna_mix(kind); fa = f"{out}/tmp/{name}.fa"; write_mix(fa, E)
    total = sum(len(s) for _, s in E)
    n = int(depth * total / (10 * 300) / max(1, umi_copies or 1))
    to_fastq([wgsim(fa, n, seed, name)], name, umi_copies)

# ---- reference, annotation, resources -----------------------------------------------
with open(f"{out}/genome.fa", "w") as f:
    for c, v in seq.items():
        f.write(f">{c}\n" + "\n".join(v[i:i + 60] for i in range(0, len(v), 60)) + "\n")
with open(f"{out}/targets.bed", "w") as f:
    f.writelines(f"{c}\t{s}\t{e}\n" for c, s, e in targets)
with open(f"{out}/genes.gtf", "w") as f:
    for c, gl in GENES.items():
        for g, s in gl:
            a = f'gene_id "{g}_ID"; gene_name "{g}"; gene_type "protein_coding";'
            t = a + f' transcript_id "{g}_T1"; transcript_name "{g}-201"; transcript_type "protein_coding";'
            f.write(f"{c}\tsim\tgene\t{s+1}\t{s+7150}\t.\t+\t.\t{a}\n{c}\tsim\ttranscript\t{s+1}\t{s+7150}\t.\t+\t.\t{t}\n")
            for i, (x, y) in enumerate(EXON):
                f.write(f'{c}\tsim\texon\t{s+x+1}\t{s+y}\t.\t+\t.\t{t} exon_number "{i+1}";\n')
                f.write(f'{c}\tsim\tCDS\t{s+x+1}\t{s+y}\t.\t+\t0\t{t} exon_number "{i+1}";\n')
sh(f"gzip -c {out}/genes.gtf > {out}/genes.gtf.gz")
hdr = ("##fileformat=VCFv4.2\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n"
       + "".join(f"##contig=<ID={c},length={CL}>\n" for c in seq) + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
def vcf(path, rows):
    with open(path[:-3], "w") as f:
        f.write(hdr)
        for c, p, r, a, af in sorted(set(rows)):
            f.write(f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\tAF={af}\n")
    sh(f"bgzip -f {path[:-3]} && tabix -f -p vcf {path}")
bg = []
for c, s, e in targets:
    for p in range(s + 25, e, 70):
        bg.append((c, p + 1, seq[c][p], comp(seq[c][p]), round(random.uniform(1e-4, 0.3), 4)))
common_rows = [(c, p, r, a, round(random.uniform(0.2, 0.5), 3)) for c, p, r, a in [G1] + COMMON]
vcf(f"{out}/gnomad.vcf.gz", common_rows + bg)
vcf(f"{out}/common.vcf.gz", common_rows + [(c, p, r, a, 0.3) for c, p, r, a, _ in bg[::3]])
vcf(f"{out}/known_sites.vcf.gz", common_rows + bg[::2])
with open(f"{out}/known_fusions.tsv", "w") as f:
    f.write("EML4\tALK\nBCR\tABL1\n")
ctx = [f"{x}[{m}]{y}" for m in ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"] for x in "ACGT" for y in "ACGT"]
with open(f"{out}/toy_signatures.tsv", "w") as f:
    f.write("Type\tSBS_toyA\tSBS_toyB\n")
    wa = [random.random() for _ in ctx]; wb = [random.random() ** 3 for _ in ctx]
    for k, x, y in zip(ctx, wa, wb):
        f.write(f"{k}\t{x/sum(wa):.6f}\t{y/sum(wb):.6f}\n")

# ---- reads -------------------------------------------------------------------------------------
simulate_dna("tumor", "P1_T", 400, 11)
simulate_dna("normal", "P1_N", 400, 12)
simulate_dna("tumor", "P3_T_umi", 400, 13, umi_copies=3)
simulate_dna("tumor", "P2_T_src", 400, 14)
# whole-genome normal (germline SV) at ~30x
with open(f"{out}/tmp/wgs.fa", "w") as f:
    for c in seq:
        f.write(f">A_{c}\n{hapA[c]}\n>B_{c}\n{hapB_wgs[c]}\n")
to_fastq([wgsim(f"{out}/tmp/wgs.fa", int(30 * 2 * 3 * CL / 300 / 2), 15, "wgs")], "WGS_N")
# tumour RNA: one transcript per gene + the EML4--ALK fusion transcript
mrna = lambda h, g: "".join(h[gpos[g][0]][gpos[g][1] + x:gpos[g][1] + y] for x, y in EXON)
with open(f"{out}/tmp/rna.fa", "w") as f:
    for g in gpos:
        for k in range(random.randint(6, 12)):
            f.write(f">{g}_{k}\n{mrna(hapT, g)}\n")
    fusion = "".join(seq["chr1"][gpos["EML4"][1] + x:gpos["EML4"][1] + y] for x, y in EXON[:2]) + \
             "".join(seq["chr3"][gpos["ALK"][1] + x:gpos["ALK"][1] + y] for x, y in EXON[1:])   # EML4 ex1-2 + ALK ex2-8
    for k in range(10):
        f.write(f">EML4_ALK_{k}\n{fusion}\n")
to_fastq([wgsim(f"{out}/tmp/rna.fa", 40000, 16, "rna")], "P1_RNA")

# imported tumour BAM: aligned here WITHOUT a read group
sh(f"samtools faidx {out}/genome.fa && mkdir -p {out}/bwamem2 && cp {out}/genome.fa {out}/bwamem2/ && "
   f"bwa-mem2 index {out}/bwamem2/genome.fa > /dev/null 2>&1")
sh(f"bwa-mem2 mem -t 4 {out}/bwamem2/genome.fa {out}/fastq/P2_T_src_R1.fastq.gz {out}/fastq/P2_T_src_R2.fastq.gz 2>/dev/null "
   f"| samtools sort -o {out}/bam/P2_T.bam - 2>/dev/null && samtools index {out}/bam/P2_T.bam")
sh(f"rm -rf {out}/tmp {out}/fastq/P2_T_src_* {out}/bwamem2/genome.fa")

# ---- sample sheets + truth ---------------------------------------------------------------------------
fq = lambda n: f"{out}/fastq/{n}_R1.fastq.gz,{out}/fastq/{n}_R2.fastq.gz"
H = "patient,sample,status,type,fastq_1,fastq_2,bam\n"
open(f"{out}/samplesheet.csv", "w").write(H + f"P1,P1_T,tumor,dna,{fq('P1_T')},\nP1,P1_N,normal,dna,{fq('P1_N')},\n"
                                          f"P1,P1_RNA,tumor,rna,{fq('P1_RNA')},\nP2,P2_T,tumor,dna,,,{out}/bam/P2_T.bam\n")
open(f"{out}/samplesheet_umi.csv", "w").write(H + f"P3,P3_T,tumor,dna,{fq('P3_T_umi')},\n")
open(f"{out}/samplesheet_wgs.csv", "w").write(H + f"W1,WGS_N,normal,dna,{fq('WGS_N')},\n")
with open(f"{out}/truth.tsv", "w") as f:
    f.write("class\tchrom\tpos\tref\talt\n")
    for v in [G1, G2]:
        f.write("germline\t" + "\t".join(map(str, v)) + "\n")
    for v in COMMON:
        f.write("common_het\t" + "\t".join(map(str, v)) + "\n")
    for v in SOM:
        f.write("somatic\t" + "\t".join(map(str, v)) + "\n")
# the deletion is left-aligned with bcftools norm, as the pipeline does
with open(f"{out}/del.vcf", "w") as f:
    f.write(hdr.replace("##INFO=<ID=AF", "##INFO=<ID=AF") + "\t".join(map(str, [DEL[0], DEL[1], ".", DEL[2], DEL[3], ".", "PASS", "AF=0.3"])) + "\n")
sh(f"bcftools norm -f {out}/genome.fa {out}/del.vcf 2>/dev/null | grep -v '^#' | cut -f1,2,4,5 | sed 's/^/somatic\t/' >> {out}/truth.tsv && rm {out}/del.vcf")
print(f"test data written to {out}")
