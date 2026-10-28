#!/usr/bin/env python3
"""
MutaXplorer - aggregate per-sample variant calls into cohort tables.

Input (--calls-dir): either one folder per sample, or a flat folder of files named
  <sid>.annotated.vcf.gz   VEP-annotated PASS calls (fallback <sid>.norm.vcf.gz)
  <sid>.mosdepth.summary.txt / <sid>.thresholds.bed.gz      coverage (mosdepth)
  <sid>.msi                                                  MSIsensor2 (somatic only)
DRAGEN outputs (<sid>.hard-filtered.vcf.gz + Nirvana JSON, TMB/MSI files) are also understood.

--mode somatic : Mutect2 calls; computes TMB and MSI status.
--mode germline: HaplotypeCaller calls; reports rare protein-altering variants + genotype.

Writes to --out: variants_all_pass.tsv, variants_filtered.tsv, top_variants.tsv, top_genes.tsv,
gene_sample_matrix.tsv, sample_metrics.tsv, filter_summary.tsv, cohort_results.xlsx
"""
import argparse, csv, glob, gzip, json, os, re, sys
from collections import Counter, defaultdict
import numpy as np
import pandas as pd

# ---------------------------------------------------------------- consequences
SEVERITY = [  # most -> least severe (Sequence Ontology terms)
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "transcript_amplification", "inframe_insertion", "inframe_deletion",
    "missense_variant", "protein_altering_variant", "splice_region_variant",
    "incomplete_terminal_codon_variant", "start_retained_variant",
    "stop_retained_variant", "synonymous_variant", "coding_sequence_variant",
    "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant", "intron_variant",
    "NMD_transcript_variant", "non_coding_transcript_variant",
    "upstream_gene_variant", "downstream_gene_variant", "intergenic_variant",
]
RANK = {c: i for i, c in enumerate(SEVERITY)}
CLASS_MAP = {
    "stop_gained": "Nonsense", "frameshift_variant": "Frameshift",
    "splice_acceptor_variant": "Splice site", "splice_donor_variant": "Splice site",
    "missense_variant": "Missense", "protein_altering_variant": "Missense",
    "inframe_insertion": "In-frame indel", "inframe_deletion": "In-frame indel",
    "start_lost": "Start/Stop lost", "stop_lost": "Start/Stop lost",
    "transcript_ablation": "Nonsense",
}
NONSILENT = set(CLASS_MAP)


def worst(csqs):
    csqs = [c for c in csqs if c]
    return min(csqs, key=lambda c: RANK.get(c, 999)) if csqs else ""


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


# ---------------------------------------------------------------- VCF reader
def read_vcf(path):
    """Yield one dict per ALT allele with sample-level metrics + any VEP/snpEff annotation."""
    csq_fields = ann_fields = None
    tumor_name, col = None, 9
    with opener(path) as fh:
        for line in fh:
            if line.startswith("##tumor_sample="):
                tumor_name = line.strip().split("=", 1)[1]
            if line.startswith("#CHROM"):
                names = line.rstrip("\n").split("\t")[9:]
                if tumor_name in names:
                    col = 9 + names.index(tumor_name)
            if line.startswith("##INFO=<ID=CSQ"):
                csq_fields = re.search(r"Format: ([^\"]+)", line).group(1).split("|")
            elif line.startswith("##INFO=<ID=ANN"):
                ann_fields = ["Allele", "Annotation", "Impact", "Gene_Name", "Gene_ID",
                              "Feature_Type", "Feature_ID", "Biotype", "Rank",
                              "HGVS.c", "HGVS.p"]
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alts, _q, filt, info = f[:8]
            info_d = {}
            for kv in info.split(";"):
                k, _, v = kv.partition("=")
                info_d[k] = v if v else True
            fmt = f[8].split(":") if len(f) > 8 else []
            smp = dict(zip(fmt, f[col].split(":"))) if len(f) > col else {}
            alts = alts.split(",")
            ad = [int(x) if x.isdigit() else 0 for x in smp.get("AD", "").split(",")] if "AD" in smp else []
            dp = smp.get("DP") or info_d.get("DP")
            dp = int(dp) if dp and str(dp).isdigit() else (sum(ad) if ad else 0)
            afs = smp.get("AF", "").split(",") if "AF" in smp else []
            germ = any("germline" in k.lower() or "germline" in str(v).lower()
                       for k, v in info_d.items() if k not in ("CSQ", "ANN"))
            for i, alt in enumerate(alts):
                alt_reads = ad[i + 1] if len(ad) > i + 1 else 0
                try:
                    vaf = float(afs[i])
                except (IndexError, ValueError):
                    vaf = alt_reads / dp if dp else 0.0
                rec = dict(chrom=chrom, pos=int(pos), ref=ref, alt=alt, filter=filt,
                           depth=dp, alt_reads=alt_reads, vaf=round(vaf, 4),
                           germline_tag=germ, genotype=smp.get("GT", ""))
                # optional VEP / snpEff annotations present in the VCF
                if csq_fields and "CSQ" in info_d:
                    rec.update(pick_vep(info_d["CSQ"], csq_fields, alt))
                elif ann_fields and "ANN" in info_d:
                    rec.update(pick_snpeff(info_d["ANN"], alt))
                yield rec


def pick_vep(csq, fields, alt):
    best = None
    for entry in csq.split(","):
        d = dict(zip(fields, entry.split("|")))
        cons = d.get("Consequence", "").split("&")
        score = (0 if d.get("CANONICAL") == "YES" else 1, RANK.get(worst(cons), 999))
        if best is None or score < best[0]:
            best = (score, d, cons)
    if not best:
        return {}
    d, cons = best[1], best[2]
    af = []
    for k in ("gnomADe_AF", "gnomAD_AF", "gnomADg_AF", "AF", "MAX_AF"):
        for x in (d.get(k) or "").split("&"):
            try:
                af.append(float(x))
            except ValueError:
                pass
    af = max(af or [0.0])
    return dict(gene=d.get("SYMBOL", ""), consequence=worst(cons),
                transcript=d.get("Feature", ""),
                hgvsc=d.get("HGVSc", "").split(":")[-1], hgvsp=d.get("HGVSp", "").split(":")[-1].replace("%3D", "="),
                pop_af=af, clinvar=d.get("CLIN_SIG", "").replace("&", ";"),
                cosmic=";".join(x for x in d.get("Existing_variation", "").split("&") if x.startswith(("COSV", "COSM"))))


def pick_snpeff(ann, alt):
    best = None
    for entry in ann.split(","):
        p = entry.split("|")
        if len(p) < 11:
            continue
        cons = p[1].split("&")
        r = RANK.get(worst(cons), 999)
        if best is None or r < best[0]:
            best = (r, p, cons)
    if not best:
        return {}
    p = best[1]
    return dict(gene=p[3], consequence=worst(best[2]), transcript=p[6],
                hgvsc=p[9], hgvsp=p[10], pop_af=0.0, clinvar="", cosmic="")


# ---------------------------------------------------------------- Nirvana JSON
def read_nirvana(path):
    """Return {(chrom,pos,ref,alt): annotation dict}. Nirvana writes one position per line."""
    out = {}
    with opener(path) as fh:
        for line in fh:
            line = line.strip().rstrip(",")
            if not line.startswith('{"chromosome"'):
                continue
            try:
                p = json.loads(line)
            except json.JSONDecodeError:
                continue
            for alt, v in zip(p.get("altAlleles", []), p.get("variants", [])):
                out[(p["chromosome"], int(p["position"]), p["refAllele"], alt)] = nirvana_variant(v)
    return out


def nirvana_variant(v):
    tx = v.get("transcripts", []) or []

    def key(t):
        return (0 if t.get("isManeSelect") else 1,
                0 if t.get("isCanonical") else 1,
                0 if t.get("source") == "RefSeq" else 1,
                RANK.get(worst(t.get("consequence", [])), 999))
    t = min(tx, key=key) if tx else {}
    pops = []
    for k in ("gnomad", "gnomad-exome", "oneKg"):
        d = v.get(k)
        if isinstance(d, dict) and d.get("allAf") is not None:
            pops.append(float(d["allAf"]))
    clin = sorted({s for c in v.get("clinvar", []) or [] for s in c.get("significance", [])})
    cosmic_n = sum(int(c.get("numSamples", 0) or 0) for c in v.get("cosmic", []) or [])
    hgvsp = t.get("hgvsp", "") or ""
    return dict(gene=t.get("hgnc", ""), consequence=worst(t.get("consequence", [])),
                transcript=t.get("transcript", ""),
                hgvsc=(t.get("hgvsc", "") or "").split(":")[-1],
                hgvsp=hgvsp.split(":")[-1].replace("%3D", "="),
                pop_af=max(pops) if pops else 0.0,
                clinvar=";".join(clin), cosmic=str(cosmic_n) if cosmic_n else "")


# ---------------------------------------------------------------- metrics
def read_kv_csv(path):
    d = {}
    if path and os.path.exists(path):
        with open(path) as fh:
            for row in csv.reader(fh):
                if len(row) >= 4:
                    d[row[2].strip()] = row[3].strip()
    return d


def sample_metrics(sdir, sid):
    m = {"sample_id": sid}
    # --- DRAGEN outputs (if present) ---
    tmb = read_kv_csv(os.path.join(sdir, f"{sid}.tmb.metrics.csv"))
    for k in ("Nonsyn TMB", "TMB"):
        if k in tmb:
            try:
                m["TMB"] = float(tmb[k]); break
            except ValueError:
                pass
    msi = os.path.join(sdir, f"{sid}.microsat_output.json")
    if os.path.exists(msi):
        with open(msi) as fh:
            m["MSI_pct_unstable"] = json.load(fh).get("PercentageUnstableSites")
    # --- MSIsensor2 ---
    msi2 = os.path.join(sdir, f"{sid}.msi")
    if os.path.exists(msi2):
        with open(msi2) as fh:
            lines = [l.split() for l in fh if l.strip()]
        if len(lines) >= 2 and len(lines[1]) >= 3:
            m["MSI_sites"] = int(lines[1][0]); m["MSI_unstable_sites"] = int(lines[1][1])
            m["MSI_pct_unstable"] = float(lines[1][2])
    # --- mosdepth ---
    summ = os.path.join(sdir, f"{sid}.mosdepth.summary.txt")
    if os.path.exists(summ):
        sm = pd.read_csv(summ, sep="\t")
        r = sm[sm["chrom"] == "total_region"]
        if len(r):
            m["mean_coverage"] = round(float(r["mean"].iloc[0]), 1)
    thr = os.path.join(sdir, f"{sid}.thresholds.bed.gz")
    if os.path.exists(thr):
        t = pd.read_csv(thr, sep="\t")
        size = (t.iloc[:, 2] - t.iloc[:, 1]).sum()
        for col in t.columns[4:]:
            x = re.sub(r"\D", "", col)
            m[f"pct_target_{x}x"] = round(100 * t[col].sum() / size, 2)
        if len(t.columns) > 4:
            m["bases_min_cov_Mb"] = round(t[t.columns[4]].sum() / 1e6, 4)  # lowest threshold
    # --- DRAGEN coverage csv ---
    covf = os.path.join(sdir, f"{sid}.qc-coverage-region-1_coverage_metrics.csv")
    for k, v in read_kv_csv(covf).items():
        if k.startswith("Average alignment coverage"):
            m["mean_coverage"] = float(v)
        if re.search(r"\[\s*100x", k):
            m["pct_target_100x"] = float(v)
    return m


def find_vcf(sdir, sid):
    for name in (f"{sid}.annotated.vcf.gz", f"{sid}.norm.vcf.gz", f"{sid}.mutect2.filtered.vcf.gz",
                 f"{sid}.hard-filtered.vcf.gz"):
        p = os.path.join(sdir, name)
        if os.path.exists(p):
            return p
    return None


VCF_SUFFIXES = (".annotated.vcf.gz", ".norm.vcf.gz", ".mutect2.filtered.vcf.gz", ".hard-filtered.vcf.gz")


def discover_samples(calls_dir):
    """[(folder, sample_id)] from per-sample sub-folders or a flat folder."""
    subs = sorted(d for d in glob.glob(os.path.join(calls_dir, "*")) if os.path.isdir(d))
    found = [(d, os.path.basename(d)) for d in subs if find_vcf(d, os.path.basename(d))]
    if found:
        return found
    sids = set()
    for f in os.listdir(calls_dir):
        for suf in VCF_SUFFIXES:
            if f.endswith(suf):
                sids.add(f[: -len(suf)])
    return [(calls_dir, sid) for sid in sorted(sids)]


TMB_CODING = NONSILENT | {"synonymous_variant", "stop_retained_variant", "start_retained_variant",
                          "coding_sequence_variant", "incomplete_terminal_codon_variant"}


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls-dir", "--dragen-dir", dest="calls_dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", choices=["somatic", "germline"], default="somatic")
    ap.add_argument("--min-dp", type=int, default=100)
    ap.add_argument("--min-alt", type=int, default=5)
    ap.add_argument("--min-vaf", type=float, default=0.05)
    ap.add_argument("--max-pop-af", type=float, default=0.01)
    ap.add_argument("--msi-high", type=float, default=20.0, help="%% unstable sites to call MSI-High")
    ap.add_argument("--msi-min-sites", type=int, default=20,
                    help="fewer evaluable MSIsensor2 sites than this -> MSI 'Insufficient sites'")
    ap.add_argument("--keep-germline-tagged", action="store_true")
    ap.add_argument("--tmb-region-mb", type=float, default=None,
                    help="fixed TMB denominator (Mb); default = per-sample target bases at the lowest coverage threshold")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    rows, metrics, missing_ann = [], [], []
    for sdir, sid in discover_samples(a.calls_dir):
        vcf = find_vcf(sdir, sid)
        if not vcf:
            print(f"[warn] {sid}: no VCF, skipped", file=sys.stderr); continue
        nj = os.path.join(sdir, f"{sid}.hard-filtered.annotated.json.gz")
        ann = read_nirvana(nj) if os.path.exists(nj) else {}
        has_ann = False
        for r in read_vcf(vcf):
            if r["filter"] not in ("PASS", "."):
                continue
            r.update(ann.get((r["chrom"], r["pos"], r["ref"], r["alt"]), {}))
            has_ann = has_ann or bool(r.get("gene"))
            r["sample_id"] = sid
            rows.append(r)
        if not has_ann:
            missing_ann.append(sid)
        metrics.append(sample_metrics(sdir, sid))
    if not rows:
        print("[warn] no PASS SNV/indel calls found - SNV tables left empty", file=sys.stderr)
        pd.DataFrame(metrics or [{"sample_id": sid} for _, sid in discover_samples(a.calls_dir)]).to_csv(
            os.path.join(a.out, "sample_metrics.tsv"), sep="\t", index=False)
        return
    if missing_ann:
        print(f"[warn] {len(missing_ann)} samples have no annotation (Nirvana JSON / CSQ / ANN). "
              "Check the VEP step (VEP_CACHE in config.sh).", file=sys.stderr)

    df = pd.DataFrame(rows)
    for c, dflt in [("genotype", ""), ("gene", ""), ("consequence", ""), ("hgvsc", ""), ("hgvsp", ""),
                    ("transcript", ""), ("pop_af", 0.0), ("clinvar", ""), ("cosmic", "")]:
        if c not in df:
            df[c] = dflt
        df[c] = df[c].fillna(dflt)
    df["variant_class"] = df["consequence"].map(CLASS_MAP).fillna("Other")
    df["variant_id"] = df["chrom"] + ":" + df["pos"].astype(str) + ":" + df["ref"] + ">" + df["alt"]
    df["protein_change"] = df["hgvsp"].str.replace(r"^p\.\(?|\)$", "", regex=True)
    df["label"] = df["gene"] + " " + df["protein_change"].where(df["protein_change"] != "", df["hgvsc"])
    df["label"] = df["label"].where(df["gene"] != "", df["variant_id"])
    df["var_type"] = ["SNV" if len(r) == 1 and len(x) == 1 else ("MNV" if len(r) == len(x) else "Indel")
                      for r, x in zip(df["ref"], df["alt"])]

    # ---- filtering cascade (first failing reason is recorded) ----------------
    reason = pd.Series("", index=df.index)
    checks = [
        ("low_depth", df["depth"] < a.min_dp),
        ("low_alt_reads", df["alt_reads"] < a.min_alt),
        ("low_vaf", df["vaf"] < a.min_vaf),
        ("population_af", df["pop_af"] > a.max_pop_af),
        ("germline_tagged", df["germline_tag"] & (a.mode == "somatic") & (not a.keep_germline_tagged)),
        ("non_coding_or_silent", ~df["consequence"].isin(NONSILENT)),
    ]
    for name, mask in checks:
        reason[(reason == "") & mask] = name
    df["filter_reason"] = reason.replace("", "REPORTABLE")
    fs = df["filter_reason"].value_counts().rename_axis("filter").reset_index(name="n_calls")
    fs.to_csv(os.path.join(a.out, "filter_summary.tsv"), sep="\t", index=False)

    cols = ["sample_id", "gene", "label", "variant_id", "chrom", "pos", "ref", "alt", "var_type",
            "consequence", "variant_class", "transcript", "hgvsc", "hgvsp", "protein_change",
            "vaf", "depth", "alt_reads", "genotype", "pop_af", "clinvar", "cosmic", "germline_tag", "filter_reason"]
    df = df[cols].sort_values(["sample_id", "chrom", "pos"])
    df.to_csv(os.path.join(a.out, "variants_all_pass.tsv"), sep="\t", index=False)
    flt = df[df["filter_reason"] == "REPORTABLE"].copy()
    flt.drop(columns="filter_reason").to_csv(os.path.join(a.out, "variants_filtered.tsv"), sep="\t", index=False)

    n_samples = len(metrics)
    # ---- top variants -------------------------------------------------------
    tv = (flt.groupby(["variant_id", "label", "gene", "variant_class", "consequence", "hgvsc", "hgvsp"])
             .agg(n_samples=("sample_id", "nunique"), median_vaf=("vaf", "median"),
                  mean_depth=("depth", lambda s: round(s.mean())), clinvar=("clinvar", "first"), cosmic=("cosmic", "first"),
                  samples=("sample_id", lambda s: ",".join(sorted(set(s)))))
             .reset_index())
    tv["pct_samples"] = (100 * tv["n_samples"] / n_samples).round(1)
    tv = tv.sort_values(["n_samples", "median_vaf"], ascending=False)
    tv.to_csv(os.path.join(a.out, "top_variants.tsv"), sep="\t", index=False)

    # ---- top genes ------------------------------------------------------------
    tg = (flt[flt["gene"] != ""].groupby("gene")
             .agg(n_samples=("sample_id", "nunique"), n_variants=("variant_id", "size"),
                  n_distinct_variants=("variant_id", "nunique"), median_vaf=("vaf", "median"))
             .reset_index())
    tg["pct_samples"] = (100 * tg["n_samples"] / n_samples).round(1)
    tg = tg.sort_values(["n_samples", "n_variants"], ascending=False)
    tg.to_csv(os.path.join(a.out, "top_genes.tsv"), sep="\t", index=False)

    # ---- oncoplot matrix ------------------------------------------------------
    g = flt[flt["gene"] != ""].groupby(["gene", "sample_id"])["variant_class"].agg(
        lambda s: "Multi hit" if s.nunique() > 1 or len(s) > 1 else s.iloc[0])
    mat = g.unstack(fill_value="")
    for m in metrics:
        if m["sample_id"] not in mat.columns:
            mat[m["sample_id"]] = ""
    mat = mat.loc[tg["gene"]] if len(tg) else mat
    mat.to_csv(os.path.join(a.out, "gene_sample_matrix.tsv"), sep="\t")

    # ---- per-sample metrics ---------------------------------------------------
    sm = pd.DataFrame(metrics)
    # TMB (computed unless DRAGEN already provided it)
    tm = df[(df["depth"] >= a.min_dp) & (df["vaf"] >= a.min_vaf) & (df["pop_af"] <= a.max_pop_af)
            & df["consequence"].isin(TMB_CODING) & ~df["germline_tag"]]
    tmb_n = tm.groupby("sample_id").size()
    if a.mode == "somatic" and ("TMB" not in sm or sm["TMB"].isna().all()):
        denom = (pd.Series(a.tmb_region_mb, index=sm.index) if a.tmb_region_mb
                 else sm.get("bases_min_cov_Mb", pd.Series(np.nan, index=sm.index)))
        sm["TMB_variants"] = sm["sample_id"].map(tmb_n).fillna(0).astype(int)
        sm["TMB"] = (sm["TMB_variants"] / denom.replace(0, np.nan)).round(2)
    cnt = flt.groupby("sample_id").size().rename("n_reportable")
    cnt_all = df.groupby("sample_id").size().rename("n_pass_calls")
    sm = sm.merge(cnt_all, on="sample_id", how="left").merge(cnt, on="sample_id", how="left")
    sm[["n_pass_calls", "n_reportable"]] = sm[["n_pass_calls", "n_reportable"]].fillna(0).astype(int)
    if a.mode == "somatic" and "MSI_pct_unstable" in sm:
        def msi_status(r):
            if pd.isna(r["MSI_pct_unstable"]):
                return ""
            if "MSI_sites" in r and pd.notna(r["MSI_sites"]) and r["MSI_sites"] < a.msi_min_sites:
                return "Insufficient sites"
            return "MSI-High" if r["MSI_pct_unstable"] >= a.msi_high else "MSS"
        sm["MSI_status"] = sm.apply(msi_status, axis=1)
    if "TMB" in sm and sm["TMB"].notna().any():
        sm["TMB_class"] = pd.cut(sm["TMB"], [-1, 10, 20, 1e9], labels=["Low (<10)", "Intermediate (10-20)", "High (>=20)"]).astype(str)
    sm.to_csv(os.path.join(a.out, "sample_metrics.tsv"), sep="\t", index=False)

    try:
        with pd.ExcelWriter(os.path.join(a.out, "cohort_results.xlsx")) as xw:
            for name, d in [("top_variants", tv), ("top_genes", tg), ("variants_filtered", flt),
                            ("sample_metrics", sm), ("filter_summary", fs)]:
                d.to_excel(xw, sheet_name=name, index=False)
    except ImportError:
        print("[info] openpyxl missing - skipping xlsx", file=sys.stderr)

    print(f"samples={n_samples}  PASS calls={len(df)}  reportable={len(flt)}  "
          f"genes={tg.shape[0]}  -> {a.out}")


if __name__ == "__main__":
    main()
