#!/usr/bin/env python3
"""
MutaXplorer - copy-number, structural-variant, fusion, signature and sample-identity summaries.

Reads (flat --calls-dir, files named <sample>.<suffix>):
  <s>.call.cns / <s>.cnr       CNVkit segments (integer CN) / bins
  <s>.manta.vcf.gz             Manta SVs (somaticSV, tumorSV or diploidSV)
  <s>.arriba.tsv               Arriba RNA fusions
  <s>.norm.vcf.gz              SNVs (for SBS96 mutational signatures)
  identity_genotypes.vcf.gz    genotypes of all DNA samples at common SNPs
Writes tables to --tables and figures + sections.html to --report:
  cnv_segments.tsv, cnv_gene_calls.tsv, sv_calls.tsv, dna_fusions.tsv, rna_fusions.tsv,
  fusion_summary.tsv, oncoplot_events.tsv (gene-level CNV/fusion events for the oncoplot),
  sbs96_counts.tsv, signature_exposures.tsv, sample_identity.tsv
"""
import argparse, base64, glob, gzip, html, io, math, os, re, subprocess, sys
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.ticker
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
BLUE, RED, GREY = "#2a78d6", "#e34948", "#f0efec"
SBS_COL = dict(zip(["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"],
                   ["#2a78d6", "#0b0b0b", "#e34948", "#898781", "#1baf7a", "#e87ba4"]))
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": AXIS,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF,
                     "axes.facecolor": SURF, "savefig.facecolor": SURF, "axes.titlelocation": "left",
                     "axes.titleweight": "bold", "axes.titlesize": 11, "xtick.color": MUTED, "ytick.color": INK2})


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def samples_with(d, suffix):
    return sorted((f[: -len(suffix)], os.path.join(d, f)) for f in os.listdir(d) if f.endswith(suffix))


def save(fig, outdir, name, figs, caption):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"), dpi=200 if ext == "png" else None, bbox_inches="tight")
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=110, bbox_inches="tight"); plt.close(fig)
    figs.append((name, caption, base64.b64encode(buf.getvalue()).decode()))


# ------------------------------------------------------------------ genes
def read_genes(path):
    g = pd.read_csv(path, sep="\t", header=None, comment="#",
                    names=["chrom", "start", "end", "gene", "strand", "biotype"][: None])
    g["biotype"] = g["biotype"].fillna("") if "biotype" in g else ""
    return g


class GeneIndex:
    def __init__(self, genes):
        self.by_chr = {c: d.sort_values("start").reset_index(drop=True) for c, d in genes.groupby("chrom")}

    def at(self, chrom, pos):
        d = self.by_chr.get(chrom)
        if d is None:
            return ""
        hit = d[(d.start <= pos) & (d.end >= pos)]
        if hit.empty:
            return ""
        pc = hit[hit.biotype == "protein_coding"]
        return (pc if len(pc) else hit).iloc[0].gene

    def overlap(self, chrom, s, e):
        d = self.by_chr.get(chrom)
        if d is None:
            return d
        return d[(d.start < e) & (d.end > s)]


# ------------------------------------------------------------------ CNV
def cnv_calls(d, genes, targets, amp_cn, figs, outdir, tables):
    segs = []
    for sid, f in samples_with(d, ".call.cns"):
        s = pd.read_csv(f, sep="\t"); s.insert(0, "sample_id", sid); segs.append(s)
    if not segs:
        return pd.DataFrame(), ""
    seg = pd.concat(segs, ignore_index=True)
    seg = seg[[c for c in ["sample_id", "chromosome", "start", "end", "log2", "cn", "depth", "probes"] if c in seg]]
    seg.to_csv(os.path.join(tables, "cnv_segments.tsv"), sep="\t", index=False)

    # genes informative for this assay = genes overlapping the calling intervals
    gi = GeneIndex(genes)
    if targets is not None:
        keep = set()
        for c, s, e in targets.itertuples(index=False):
            ov = gi.overlap(c, s, e)
            if ov is not None:
                keep.update(ov.gene)
        genes = genes[genes.gene.isin(keep)]
    rows = []
    for sid, sd in seg.groupby("sample_id"):
        for c, gd in genes.groupby("chrom"):
            cs = sd[sd.chromosome == c]
            if cs.empty:
                continue
            for g in gd.itertuples(index=False):
                ov = cs[(cs.start < g.end) & (cs.end > g.start)]
                if ov.empty:
                    continue
                ovl = np.minimum(ov.end, g.end) - np.maximum(ov.start, g.start)
                best = ov.loc[ovl.idxmax()]
                cn = int(best.cn) if "cn" in best and pd.notna(best.cn) else int(round(2 * 2 ** best.log2))
                call = ("Amplification" if cn >= amp_cn else "Gain" if cn > 2 else
                        "Deep deletion" if cn == 0 else "Loss" if cn < 2 else "")
                if call:
                    rows.append(dict(sample_id=sid, gene=g.gene, chrom=c, start=g.start, end=g.end,
                                     cn=cn, log2=round(float(best.log2), 3), call=call))
    gc = pd.DataFrame(rows, columns=["sample_id", "gene", "chrom", "start", "end", "cn", "log2", "call"])
    gc.to_csv(os.path.join(tables, "cnv_gene_calls.tsv"), sep="\t", index=False)

    # figure: genome-wide segment heatmap
    chroms = list(dict.fromkeys(seg.chromosome))
    clen = seg.groupby("chromosome").end.max().reindex(chroms)
    off = dict(zip(chroms, np.concatenate([[0], np.cumsum(clen.values)[:-1]])))
    sids = sorted(seg.sample_id.unique())
    fig, ax = plt.subplots(figsize=(11, 0.28 * len(sids) + 1.6))
    cmap = LinearSegmentedColormap.from_list("div", [BLUE, GREY, RED])
    norm = plt.Normalize(-1.5, 1.5)
    for i, s in enumerate(sids):
        for r in seg[seg.sample_id == s].itertuples():
            ax.add_patch(plt.Rectangle((off[r.chromosome] + r.start, i), r.end - r.start, 0.9,
                                       color=cmap(norm(max(-1.5, min(1.5, r.log2)))), lw=0))
    for c in chroms:
        ax.axvline(off[c], color=AXIS, lw=0.6)
    ax.set_xlim(0, float(clen.sum())); ax.set_ylim(len(sids), 0)
    ax.set_yticks(np.arange(len(sids)) + 0.45); ax.set_yticklabels(sids, fontsize=7 if len(sids) > 40 else 8)
    ax.set_xticks([off[c] + clen[c] / 2 for c in chroms]); ax.set_xticklabels([re.sub("^chr", "", c) for c in chroms], fontsize=7)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); cb = fig.colorbar(sm, ax=ax, shrink=0.8, pad=0.01)
    cb.set_label("segment log2 copy ratio", color=INK2)
    ax.set_title("Genome-wide copy number (red gain, blue loss)")
    save(fig, outdir, "cnv_genome_heatmap", figs, "Copy-number segments per sample (CNVkit)")

    if not gc.empty:
        top = (gc[gc.call.isin(["Amplification", "Deep deletion", "Loss"])]
               .groupby(["gene", "call"]).sample_id.nunique().unstack(fill_value=0))
        if not top.empty:
            top = top.loc[top.sum(axis=1).sort_values(ascending=False).index[:25]].iloc[::-1]
            fig, ax = plt.subplots(figsize=(7, 0.3 * len(top) + 1.4)); left = np.zeros(len(top))
            for call, col in (("Amplification", RED), ("Loss", "#86b6ef"), ("Deep deletion", BLUE)):
                if call in top:
                    ax.barh(range(len(top)), top[call], left=left, color=col, height=0.7, label=call); left += top[call].values
            ax.set_yticks(range(len(top))); ax.set_yticklabels(top.index, fontstyle="italic", color=INK)
            ax.set_xlabel("Samples"); ax.set_title("Most frequent gene-level copy-number events")
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
            ax.legend(frameon=False, fontsize=8, loc="lower right"); ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
            save(fig, outdir, "cnv_recurrent_genes", figs, "Genes with amplification / loss across samples")
    return gc, f"{len(sids)} samples, {len(gc)} gene-level events"


# ------------------------------------------------------------------ SV + DNA fusions
BND_ALT = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def read_manta(path, sid):
    rows, seen = [], set()
    with opener(path) as fh:
        names = []
        for l in fh:
            if l.startswith("#CHROM"):
                names = l.rstrip("\n").split("\t")[9:]
            if l.startswith("#"):
                continue
            f = l.rstrip("\n").split("\t")
            if f[6] != "PASS":
                continue
            info = dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in f[7].split(";"))
            svt = info.get("SVTYPE", "")
            c1, p1 = f[0], int(f[1])
            if svt == "BND":
                m = BND_ALT.search(f[4])
                if not m:
                    continue
                c2, p2 = m.group(1), int(m.group(2))
            else:
                c2, p2 = c1, int(info.get("END", p1))
            key = tuple(sorted([(c1, p1), (c2, p2)]))
            if key in seen:
                continue
            seen.add(key)
            fmt = f[8].split(":") if len(f) > 8 else []
            col = 9 + (names.index(sid) if sid in names else len(names) - 1)
            smp = dict(zip(fmt, f[col].split(":"))) if len(f) > col else {}
            alt = lambda k: int(smp.get(k, "0,0").split(",")[-1]) if smp.get(k) else 0
            rows.append(dict(sample_id=sid, svtype=svt, chrom1=c1, pos1=p1, chrom2=c2, pos2=p2,
                             size=abs(p2 - p1) if c1 == c2 else None, PR_alt=alt("PR"), SR_alt=alt("SR"),
                             somatic_score=info.get("SOMATICSCORE", "")))
    return rows


def sv_and_fusions(d, gi, known, figs, outdir, tables, min_support):
    rows = []
    for sid, f in samples_with(d, ".manta.vcf.gz"):
        rows += read_manta(f, sid)
    sv = pd.DataFrame(rows, columns=["sample_id", "svtype", "chrom1", "pos1", "chrom2", "pos2", "size",
                                     "PR_alt", "SR_alt", "somatic_score"])
    if sv.empty:
        sv.to_csv(os.path.join(tables, "sv_calls.tsv"), sep="\t", index=False)
        return sv, pd.DataFrame(columns=["sample_id", "gene1", "gene2"])
    sv["support"] = sv.PR_alt + sv.SR_alt
    sv["gene1"] = [gi.at(c, p) for c, p in zip(sv.chrom1, sv.pos1)]
    sv["gene2"] = [gi.at(c, p) for c, p in zip(sv.chrom2, sv.pos2)]
    sv.to_csv(os.path.join(tables, "sv_calls.tsv"), sep="\t", index=False)
    fu = sv[(sv.gene1 != "") & (sv.gene2 != "") & (sv.gene1 != sv.gene2) & (sv.support >= min_support)].copy()
    fu["fusion"] = [f"{a}--{b}" for a, b in zip(fu.gene1, fu.gene2)]
    fu["known_fusion"] = [(a, b) in known or (b, a) in known for a, b in zip(fu.gene1, fu.gene2)]
    fu.to_csv(os.path.join(tables, "dna_fusions.tsv"), sep="\t", index=False)
    cnt = sv.groupby(["sample_id", "svtype"]).size().unstack(fill_value=0)
    fig, ax = plt.subplots(figsize=(min(14, 3 + 0.18 * len(cnt)), 3.4)); bottom = np.zeros(len(cnt))
    cols = dict(zip(["DEL", "DUP", "INV", "INS", "BND"], [BLUE, RED, "#eda100", "#1baf7a", "#6250d6"]))
    for t in cnt.columns:
        ax.bar(range(len(cnt)), cnt[t], bottom=bottom, color=cols.get(t, MUTED), label=t, width=0.8); bottom += cnt[t].values
    ax.set_xticks(range(len(cnt))); ax.set_xticklabels(cnt.index, rotation=90, fontsize=7 if len(cnt) > 30 else 8)
    ax.set_ylabel("PASS SVs"); ax.set_title("Structural variants per sample (Manta)")
    ax.legend(frameon=False, fontsize=8, ncol=5); ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
    save(fig, outdir, "sv_per_sample", figs, "PASS structural variants by type")
    return sv, fu


def rna_fusions(d, conf, tables):
    rows = []
    for sid, f in samples_with(d, ".arriba.tsv"):
        t = pd.read_csv(f, sep="\t")
        t.columns = [c.lstrip("#") for c in t.columns]
        if t.empty:
            continue
        t = t[t.confidence.isin(conf)]
        for r in t.itertuples(index=False):
            g1, g2 = str(r.gene1).split(",")[0].split("(")[0], str(r.gene2).split(",")[0].split("(")[0]
            rows.append(dict(sample_id=sid, gene1=g1, gene2=g2, fusion=f"{g1}--{g2}", breakpoint1=r.breakpoint1,
                             breakpoint2=r.breakpoint2, type=r.type, split_reads=int(r.split_reads1) + int(r.split_reads2),
                             discordant_mates=int(r.discordant_mates), confidence=r.confidence,
                             reading_frame=r.reading_frame))
    rf = pd.DataFrame(rows, columns=["sample_id", "gene1", "gene2", "fusion", "breakpoint1", "breakpoint2", "type",
                                     "split_reads", "discordant_mates", "confidence", "reading_frame"])
    rf.to_csv(os.path.join(tables, "rna_fusions.tsv"), sep="\t", index=False)
    return rf


def fusion_summary(fu, rf, meta, known, figs, outdir, tables):
    pat = dict(zip(meta.id, meta.patient))
    recs = defaultdict(lambda: dict(dna=set(), rna=set(), name=None))
    for r in rf.itertuples():      # RNA first: Arriba reports the 5' gene as gene1
        k = (pat.get(r.sample_id, r.sample_id), frozenset([r.gene1, r.gene2]))
        recs[k]["rna"].add(r.sample_id); recs[k]["name"] = recs[k]["name"] or (r.gene1, r.gene2)
    for r in fu.itertuples():
        k = (pat.get(r.sample_id, r.sample_id), frozenset([r.gene1, r.gene2]))
        recs[k]["dna"].add(r.sample_id); recs[k]["name"] = recs[k]["name"] or (r.gene1, r.gene2)
    rows = [dict(patient=p, fusion="--".join(v["name"]), evidence="DNA + RNA" if v["dna"] and v["rna"] else "DNA" if v["dna"] else "RNA",
                 dna_samples=",".join(sorted(v["dna"])), rna_samples=",".join(sorted(v["rna"])),
                 known_fusion=(v["name"] in known or v["name"][::-1] in known)) for (p, g), v in recs.items()]
    fs = pd.DataFrame(rows, columns=["patient", "fusion", "evidence", "dna_samples", "rna_samples", "known_fusion"])
    fs.to_csv(os.path.join(tables, "fusion_summary.tsv"), sep="\t", index=False)
    if not fs.empty:
        c = fs.groupby(["fusion", "evidence"]).patient.nunique().unstack(fill_value=0)
        c = c.loc[c.sum(axis=1).sort_values(ascending=False).index[:20]].iloc[::-1]
        fig, ax = plt.subplots(figsize=(7, 0.32 * len(c) + 1.4)); left = np.zeros(len(c))
        for ev, col in (("DNA + RNA", BLUE), ("DNA", "#86b6ef"), ("RNA", MUTED)):
            if ev in c:
                ax.barh(range(len(c)), c[ev], left=left, color=col, height=0.7, label=ev); left += c[ev].values
        ax.set_yticks(range(len(c))); ax.set_yticklabels(c.index, color=INK); ax.set_xlabel("Patients")
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        ax.set_title("Recurrent fusions and their evidence"); ax.legend(frameon=False, fontsize=8, loc="lower right")
        ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
        save(fig, outdir, "fusions_recurrent", figs, "Gene fusions by patient (Manta DNA, Arriba RNA)")
    return fs


# ------------------------------------------------------------------ signatures
CTX = [f"{a}[{m}]{b}" for m in SBS_COL for a in "ACGT" for b in "ACGT"]
COMP = str.maketrans("ACGT", "TGCA")


def sbs96(d, fasta, min_dp, figs, outdir, tables):
    sites = []
    for sid, f in samples_with(d, ".norm.vcf.gz"):
        with gzip.open(f, "rt") as fh:
            for l in fh:
                if l[0] == "#":
                    continue
                x = l.split("\t")
                if x[6] == "PASS" and len(x[3]) == 1 and len(x[4]) == 1 and x[3] in "ACGT" and x[4] in "ACGT":
                    sites.append((sid, x[0], int(x[1]), x[3], x[4]))
    if not sites or not fasta:
        return pd.DataFrame()
    reg = os.path.join(tables, ".ctx_regions.txt")
    with open(reg, "w") as f:
        for s in sites:
            f.write(f"{s[1]}:{s[2]-1}-{s[2]+1}\n")
    out = subprocess.run(["samtools", "faidx", "-r", reg, fasta], capture_output=True, text=True, check=True).stdout
    os.remove(reg)
    seqs = [s for s in re.split(r">[^\n]*\n", out) if s.strip()]
    counts = defaultdict(lambda: dict.fromkeys(CTX, 0))
    for (sid, c, p, r, a), sq in zip(sites, seqs):
        tri = sq.replace("\n", "").upper()
        if len(tri) != 3 or "N" in tri:
            continue
        if r in "GA":
            tri, r, a = tri.translate(COMP)[::-1], r.translate(COMP), a.translate(COMP)
        counts[sid][f"{tri[0]}[{r}>{a}]{tri[2]}"] += 1
    m = pd.DataFrame(counts).reindex(CTX).fillna(0).astype(int)
    m.index.name = "Type"; m.to_csv(os.path.join(tables, "sbs96_counts.tsv"), sep="\t")
    tot = m.sum(axis=1); frac = 100 * tot / max(tot.sum(), 1)
    fig, ax = plt.subplots(figsize=(12, 3))
    ax.bar(range(96), frac, color=[SBS_COL[k[2:5]] for k in CTX], width=0.8)
    for i, mt in enumerate(SBS_COL):
        ax.text(16 * i + 7.5, frac.max() * 1.08 + 0.2, mt, ha="center", fontsize=9, color=INK, fontweight="bold")
    ax.set_xticks(range(96)); ax.set_xticklabels([k[0] + k[2] + k[-1] for k in CTX], rotation=90, fontsize=5, family="monospace")
    ax.set_ylabel("% of SNVs"); ax.set_ylim(0, frac.max() * 1.2 + 0.5)
    ax.set_title(f"Cohort SBS96 mutational profile ({int(tot.sum())} SNVs, {m.shape[1]} samples)", pad=18)
    ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
    save(fig, outdir, "signatures_sbs96", figs, "Trinucleotide-context profile of somatic SNVs")
    return m


def fit_signatures(m, sigfile, min_snvs, figs, outdir, tables):
    from scipy.optimize import nnls
    S = pd.read_csv(sigfile, sep="\t", index_col=0)
    S = S.reindex(CTX)
    if S.isna().any().any():
        print("[warn] signature matrix rows do not match SBS96 'A[C>A]A' labels - fit skipped", file=sys.stderr)
        return
    rows = []
    for sid in m.columns:
        v = m[sid].values.astype(float)
        if v.sum() < min_snvs:
            rows.append(dict(sample_id=sid, n_snvs=int(v.sum()), cosine=np.nan)); continue
        w, _ = nnls(S.values, v / v.sum())
        rec = S.values @ w
        cos = float(rec @ (v / v.sum()) / (np.linalg.norm(rec) * np.linalg.norm(v / v.sum()) + 1e-12))
        r = dict(sample_id=sid, n_snvs=int(v.sum()), cosine=round(cos, 3))
        r.update({k: round(x / max(w.sum(), 1e-12), 3) for k, x in zip(S.columns, w)})
        rows.append(r)
    ex = pd.DataFrame(rows); ex.to_csv(os.path.join(tables, "signature_exposures.tsv"), sep="\t", index=False)
    sig_cols = [c for c in S.columns if c in ex and ex[c].fillna(0).max() >= 0.05]
    ok = ex.dropna(subset=["cosine"])
    if ok.empty or not sig_cols:
        return
    top = ok[sig_cols].mean().sort_values(ascending=False).index[:8]
    pal = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
    fig, ax = plt.subplots(figsize=(min(14, 3 + 0.2 * len(ok)), 3.6)); bottom = np.zeros(len(ok))
    for c, col in zip(top, pal):
        ax.bar(range(len(ok)), ok[c].fillna(0), bottom=bottom, color=col, label=c, width=0.85); bottom += ok[c].fillna(0).values
    ax.bar(range(len(ok)), 1 - bottom, bottom=bottom, color=GREY, label="other", width=0.85)
    ax.set_xticks(range(len(ok))); ax.set_xticklabels(ok.sample_id, rotation=90, fontsize=7)
    ax.set_ylabel("Relative contribution"); ax.set_title("Signature contributions (NNLS fit)")
    ax.legend(frameon=False, fontsize=7.5, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.35))
    save(fig, outdir, "signatures_exposures", figs, "Fitted mutational-signature contributions")


# ------------------------------------------------------------------ identity
def identity(d, meta, tables, min_sites=20, min_dp=8, mismatch=0.03, match=0.01):
    """Pairwise genotype comparison at common SNPs. LOH in a tumour turns het into hom, so the
    decision uses the opposite-homozygote rate (0/0 vs 1/1) over informative sites, which stays ~0
    for two samples of one person and is ~0.1 for unrelated people at common SNPs."""
    f = os.path.join(d, "identity_genotypes.vcf.gz")
    if not os.path.exists(f):
        return pd.DataFrame()
    names, G = [], []
    with gzip.open(f, "rt") as fh:
        for l in fh:
            if l.startswith("#CHROM"):
                names = l.rstrip("\n").split("\t")[9:]
            if l[0] == "#":
                continue
            x = l.rstrip("\n").split("\t")
            if "," in x[4]:
                continue
            fmt = x[8].split(":"); gi = fmt.index("GT"); di = fmt.index("DP") if "DP" in fmt else None
            row = []
            for v in x[9:]:
                p = v.split(":")
                dp = int(p[di]) if di is not None and len(p) > di and p[di].isdigit() else 0
                gt = p[gi].replace("|", "/")
                row.append({"0/0": 0, "0/1": 1, "1/1": 2}.get(gt, -1) if dp >= min_dp else -1)
            G.append(row)
    if not G:
        return pd.DataFrame()
    G = np.array(G); pat = dict(zip(meta.id, meta.patient)); rows = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a_, b_ = G[:, i], G[:, j]
            ok = (a_ >= 0) & (b_ >= 0) & ((a_ > 0) | (b_ > 0))
            n = int(ok.sum())
            opp = int((ok & (np.abs(a_ - b_) == 2)).sum())
            disc = int((ok & (a_ != b_)).sum())
            rate = opp / n if n else np.nan
            x, y = names[i], names[j]; same = pat.get(x, x) == pat.get(y, y)
            if n < min_sites:
                st = "too few sites"
            elif same and rate > mismatch:
                st = "MISMATCH: same patient but different individuals"
            elif not same and rate < match:
                st = "UNEXPECTED MATCH: different patients look like one individual"
            else:
                st = "OK"
            rows.append(dict(sample_a=x, sample_b=y, patient_a=pat.get(x, ""), patient_b=pat.get(y, ""),
                             informative_sites=n, opposite_hom_rate=round(rate, 4),
                             genotype_discordance=round(disc / n, 4) if n else np.nan, same_patient=same, status=st))
    t = pd.DataFrame(rows).sort_values(["status", "opposite_hom_rate"])
    t.to_csv(os.path.join(tables, "sample_identity.tsv"), sep="\t", index=False)
    return t


# ------------------------------------------------------------------ html
def table_html(df, tid, n=300):
    if df is None or df.empty:
        return "<p class='sub'>No records.</p>"
    df = df.head(n)
    head = "".join(f"<th onclick=\"sortT('{tid}',{i})\">{html.escape(str(c))}</th>" for i, c in enumerate(df.columns))
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(v))}</td>" for v in r) + "</tr>" for r in df.astype(str).values)
    return (f"<input placeholder='Filter' oninput=\"filt('{tid}',this.value)\">"
            f"<div class='tw'><table id='{tid}'><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls-dir", required=True); ap.add_argument("--tables", required=True)
    ap.add_argument("--report", required=True); ap.add_argument("--mode", default="somatic")
    ap.add_argument("--samples", help="TSV id,patient,status,type"); ap.add_argument("--gene-bed")
    ap.add_argument("--targets"); ap.add_argument("--fasta"); ap.add_argument("--known-fusions")
    ap.add_argument("--signatures"); ap.add_argument("--signature-min-snvs", type=int, default=50)
    ap.add_argument("--amp-cn", type=int, default=4); ap.add_argument("--sv-min-support", type=int, default=4)
    ap.add_argument("--rna-confidence", default="high,medium"); ap.add_argument("--min-dp", type=int, default=20)
    a = ap.parse_args()
    os.makedirs(a.tables, exist_ok=True); os.makedirs(a.report, exist_ok=True)
    meta = (pd.read_csv(a.samples, sep="\t", names=["id", "patient", "status", "type"]) if a.samples and os.path.exists(a.samples)
            else pd.DataFrame(columns=["id", "patient", "status", "type"]))
    genes = read_genes(a.gene_bed) if a.gene_bed and os.path.getsize(a.gene_bed) else pd.DataFrame(columns=["chrom", "start", "end", "gene", "strand", "biotype"])
    gi = GeneIndex(genes)
    targets = (pd.read_csv(a.targets, sep="\t", header=None, usecols=[0, 1, 2], names=["c", "s", "e"])
               if a.targets and os.path.exists(a.targets) else None)
    known = set()
    if a.known_fusions and os.path.exists(a.known_fusions):
        for l in opener(a.known_fusions):
            x = l.rstrip("\n").split("\t")
            if len(x) >= 2 and not x[0].startswith("#"):
                known.add((x[0], x[1]))
    figs, sections, events = [], [], []

    gc, _ = cnv_calls(a.calls_dir, genes, targets, a.amp_cn, figs, a.report, a.tables)
    if not gc.empty:
        events += [dict(sample_id=r.sample_id, gene=r.gene, event=r.call) for r in gc.itertuples()
                   if r.call in ("Amplification", "Loss", "Deep deletion")]
        sections.append(("Copy-number alterations", ["cnv_genome_heatmap", "cnv_recurrent_genes"],
                         table_html(gc.sort_values(["call", "sample_id"]), "cnvt")))
    sv, fu = sv_and_fusions(a.calls_dir, gi, known, figs, a.report, a.tables, a.sv_min_support)
    rf = rna_fusions(a.calls_dir, a.rna_confidence.split(","), a.tables)
    if not sv.empty or not rf.empty:
        fs = fusion_summary(fu, rf, meta, known, figs, a.report, a.tables)
        dna_tumours = defaultdict(list)
        for r in meta.itertuples():
            if r.type == "dna" and r.status == "tumor":
                dna_tumours[r.patient].append(r.id)
        pat = dict(zip(meta.id, meta.patient))
        for r in list(fu.itertuples()) + list(rf.itertuples()):
            targets_ = dna_tumours.get(pat.get(r.sample_id), []) if r.sample_id not in set(fu.sample_id) else [r.sample_id]
            for sid in targets_ or [r.sample_id]:
                events += [dict(sample_id=sid, gene=r.gene1, event="Fusion"), dict(sample_id=sid, gene=r.gene2, event="Fusion")]
        body = ("<h3>Fusions (DNA + RNA, per patient)</h3>" + table_html(fs, "fust") +
                "<h3>Structural variants (Manta, PASS)</h3>" + table_html(sv, "svt") +
                ("<h3>RNA fusions (Arriba)</h3>" + table_html(rf, "rnat") if not rf.empty else ""))
        sections.append(("Structural variants and fusions", ["sv_per_sample", "fusions_recurrent"], body))
    if a.mode == "somatic":
        m = sbs96(a.calls_dir, a.fasta, a.min_dp, figs, a.report, a.tables)
        if not m.empty:
            if a.signatures and os.path.exists(a.signatures):
                fit_signatures(m, a.signatures, a.signature_min_snvs, figs, a.report, a.tables)
            ex = os.path.join(a.tables, "signature_exposures.tsv")
            sections.append(("Mutational signatures", ["signatures_sbs96", "signatures_exposures"],
                             table_html(pd.read_csv(ex, sep="\t"), "sigt") if os.path.exists(ex) else
                             "<p class='sub'>Supply <code>--signatures</code> (COSMIC SBS matrix) to fit signature contributions.</p>"))
    idt = identity(a.calls_dir, meta, a.tables)
    if not idt.empty:
        n_bad = int((~idt.status.isin(["OK", "too few sites"])).sum())
        sections.append(("Sample identity", [], f"<p class='sub'>{n_bad} flagged pair(s). Same-patient pairs should match "
                                                 f"(discordance &le; 0.10); different patients should not (&ge; 0.05).</p>" + table_html(idt, "idt")))
    pd.DataFrame(events, columns=["sample_id", "gene", "event"]).drop_duplicates().to_csv(
        os.path.join(a.tables, "oncoplot_events.tsv"), sep="\t", index=False)

    figmap = {n: (c, b) for n, c, b in figs}
    out = []
    for title, fnames, body in sections:
        imgs = "".join(f"<figure><img src='data:image/png;base64,{figmap[n][1]}' alt='{n}'><figcaption>{html.escape(figmap[n][0])}</figcaption></figure>"
                       for n in fnames if n in figmap)
        out.append(f"<h2>{html.escape(title)}</h2>{imgs}{body}")
    with open(os.path.join(a.report, "sections.html"), "w") as fh:
        fh.write("\n".join(out))
    print(f"structural summary: {len(sections)} sections, {len(figs)} figures", file=sys.stderr)


if __name__ == "__main__":
    main()
