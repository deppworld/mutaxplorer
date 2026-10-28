#!/usr/bin/env python3
"""
STEP 04 - Visualise cohort results (top variants / genes) + HTML report.

Input : results/cohort/*.tsv from 03_aggregate_cohort.py
Output: results/report/
    fig1_oncoplot.png/.pdf          top genes x samples, TMB on top, % mutated on right
    fig2_top_variants.png/.pdf      most recurrent variants (#samples)
    fig3_vaf_top_variants.png/.pdf  VAF spread of recurrent variants
    fig4_variant_summary.png/.pdf   variant classes + SNV substitution spectrum
    fig5_tmb_msi.png/.pdf           TMB distribution and MSI vs TMB
    fig6_lollipop_<GENE>.png/.pdf   protein-position hotspots for top genes
    fig7_cooccurrence.png/.pdf      co-occurrence / mutual exclusivity of top genes
    TSO500_cohort_report.html       self-contained report (figures + searchable tables)
"""
import argparse, base64, html, io, os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Patch
from scipy.stats import fisher_exact

# ---- palette (validated categorical order, fixed; never cycled) --------------
CLASS_ORDER = ["Missense", "Nonsense", "Frameshift", "Splice site",
               "In-frame indel", "Start/Stop lost", "Multi hit"]
CLASS_COL = dict(zip(CLASS_ORDER, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                                   "#e87ba4", "#008300", "#6250d6"]))
EMPTY = "#f0efec"; INK = "#0b0b0b"; INK2 = "#52514e"; MUTED = "#898781"
GRID = "#e1e0d9"; AXIS = "#c3c2b7"; SURF = "#fcfcfb"; BLUE = "#2a78d6"; RED = "#e34948"
SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": AXIS,
    "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF,
    "axes.facecolor": SURF, "savefig.facecolor": SURF, "axes.titlecolor": INK,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
})


def save(fig, out, name, figs):
    for ext in ("png", "pdf"):  # noqa
        fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=200 if ext == "png" else None,
                    bbox_inches="tight")
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    figs.append((name, base64.b64encode(buf.getvalue()).decode()))
    plt.close(fig)


def grid_x(ax):
    ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)


# ------------------------------------------------------------------ fig 1
EVENT_MARK = {"Amplification": ("^", "Amplification"), "Loss": ("v", "Loss"),
              "Deep deletion": ("v", "Deep deletion"), "Fusion": ("D", "Fusion")}


def oncoplot(mat, sm, top_n, out, figs, events=None):
    """SNV/indel classes as colour; CNV and fusion events as markers on the same cell."""
    ev = events if events is not None and len(events) else pd.DataFrame(columns=["sample_id", "gene", "event"])
    samples = list(dict.fromkeys(list(mat.columns) + list(sm.get("sample_id", [])) + list(ev.sample_id)))
    altered = {}
    for g in mat.index:
        altered[g] = set(mat.columns[(mat.loc[g] != "").values])
    for r in ev.itertuples():
        altered.setdefault(r.gene, set()).add(r.sample_id)
    genes = sorted(altered, key=lambda g: (-len(altered[g]), g))[:top_n]
    genes = [g for g in genes if altered[g]]
    if not genes:
        return
    m = mat.reindex(index=genes, columns=samples).replace("", np.nan)
    binm = pd.DataFrame({s: [int(s in altered[g]) for g in genes] for s in samples}, index=genes)
    # memo-sort samples: lexicographic on mutation pattern of ranked genes
    order = sorted(binm.columns, key=lambda s: tuple(-binm[s].values), reverse=False)
    m, binm = m[order], binm[order]
    ns = len(order)
    codes = {c: i + 1 for i, c in enumerate(CLASS_ORDER)}
    arr = m.apply(lambda col: col.map(lambda v: codes.get(v, 0) if isinstance(v, str) else 0)).values
    cmap = ListedColormap([EMPTY] + [CLASS_COL[c] for c in CLASS_ORDER])

    h = max(4, 0.22 * len(genes) + 2)
    fig = plt.figure(figsize=(min(22, 4 + ns * 0.085), h))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 6], width_ratios=[8, 1.2], hspace=0.04, wspace=0.03)
    ax = fig.add_subplot(gs[1, 0])
    ax.pcolormesh(arr, cmap=cmap, vmin=0, vmax=len(CLASS_ORDER), edgecolors=SURF, linewidth=0.4)
    ax.set_ylim(len(genes), 0)
    ax.set_yticks(np.arange(len(genes)) + 0.5)
    ax.set_yticklabels(genes, fontstyle="italic", color=INK)
    ax.set_xticks([]); ax.set_xlabel(f"{ns} samples (sorted by alteration pattern)")
    gidx, sidx = {g: i for i, g in enumerate(genes)}, {s: j for j, s in enumerate(order)}
    used = set()
    for r in ev.itertuples():
        if r.gene in gidx and r.sample_id in sidx and r.event in EVENT_MARK:
            mk, lab = EVENT_MARK[r.event]
            off = {"Fusion": 0.22}.get(r.event, -0.22)
            ax.scatter(sidx[r.sample_id] + 0.5, gidx[r.gene] + 0.5 + off * 0, marker=mk, s=18 if ns < 80 else 9,
                       facecolor=INK if r.event != "Loss" else "none", edgecolor=INK, linewidth=0.8, zorder=3)
            used.add(r.event)
    for s in ax.spines.values():
        s.set_visible(False)

    axt = fig.add_subplot(gs[0, 0], sharex=ax)
    if "TMB" in sm:
        tmb = sm.set_index("sample_id").reindex(order)["TMB"].fillna(0).values
        axt.bar(np.arange(ns) + 0.5, tmb, width=0.8, color=MUTED)
        axt.set_ylabel("TMB\n(mut/Mb)", fontsize=8)
    axt.tick_params(labelbottom=False, bottom=False)
    axt.set_title(f"Top {len(genes)} altered genes across the cohort", pad=8)
    axt.spines["bottom"].set_visible(False)

    axr = fig.add_subplot(gs[1, 1], sharey=ax)
    pct = 100 * binm.sum(axis=1).values / max(ns, 1)
    left = np.zeros(len(genes))
    for c in CLASS_ORDER:  # stacked by class
        v = 100 * (m == c).sum(axis=1).values / ns
        axr.barh(np.arange(len(genes)) + 0.5, v, left=left, height=0.7, color=CLASS_COL[c],
                 edgecolor=SURF, linewidth=0.5)
        left += v
    axr.tick_params(labelleft=False, left=False)
    axr.set_xlabel("% samples"); axr.set_xlim(0, max(left.max() * 1.35, 5))
    grid_x(axr)

    for i, g in enumerate(genes):   # bar = any alteration; the coloured part = SNV/indel classes
        rest = 100 * len(altered[g]) / ns - left[i]
        if rest > 0:
            axr.barh(i + 0.5, rest, left=left[i], height=0.7, color=MUTED, edgecolor=SURF, linewidth=0.5)
            left[i] += rest
    for i, p in enumerate(pct):
        axr.text(left[i] + 0.8, i + 0.5, f"{p:.0f}%", va="center", fontsize=7, color=INK2)
    axr.set_xlim(0, max(left.max() * 1.35, 5))
    present = [c for c in CLASS_ORDER if (m == c).values.any()]
    from matplotlib.lines import Line2D
    handles = [Patch(color=CLASS_COL[c], label=c) for c in present]
    handles += [Line2D([], [], marker=EVENT_MARK[e][0], ls="", markerfacecolor=INK if e != "Loss" else "none",
                       markeredgecolor=INK, label=e) for e in ["Amplification", "Loss", "Deep deletion", "Fusion"] if e in used]
    fig.legend(handles=handles,
               loc="lower center", ncol=min(len(handles), 6), frameon=False,
               bbox_to_anchor=(0.45, 0.0 if len(genes) > 15 else -0.08))
    save(fig, out, "fig1_oncoplot", figs)


# ------------------------------------------------------------------ fig 2
def top_variants_bar(tv, n, ns, out, figs):
    d = tv.head(n).iloc[::-1]
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 0.3 * len(d) + 1.4))
    y = np.arange(len(d))
    ax.barh(y, d["n_samples"], height=0.7, color=[CLASS_COL.get(c, MUTED) for c in d["variant_class"]])
    ax.set_yticks(y); ax.set_yticklabels(d["label"], color=INK)
    for yi, (k, p) in enumerate(zip(d["n_samples"], d["pct_samples"])):
        ax.text(k + d["n_samples"].max() * 0.01, yi, f"{k}  ({p}%)", va="center", fontsize=7.5, color=INK2)
    ax.set_xlabel(f"Samples carrying the variant (of {ns})")
    ax.set_xlim(0, d["n_samples"].max() * 1.18)
    ax.set_title(f"Top {len(d)} recurrent variants")
    present = [c for c in CLASS_ORDER if c in set(d["variant_class"])]
    ax.legend(handles=[Patch(color=CLASS_COL[c], label=c) for c in present], frameon=False,
              loc="lower right", fontsize=8)
    grid_x(ax)
    save(fig, out, "fig2_top_variants", figs)


# ------------------------------------------------------------------ fig 3
def vaf_plot(flt, tv, n, out, figs, lod=5.0):
    labels = list(tv.head(n)["variant_id"])
    d = flt[flt["variant_id"].isin(labels)]
    if d.empty:
        return
    lab = dict(zip(tv["variant_id"], tv["label"]))
    cls = dict(zip(tv["variant_id"], tv["variant_class"]))
    labels = labels[::-1]
    fig, ax = plt.subplots(figsize=(8, 0.3 * len(labels) + 1.4))
    rng = np.random.default_rng(1)
    for i, v in enumerate(labels):
        x = d.loc[d["variant_id"] == v, "vaf"].values * 100
        ax.boxplot(x, positions=[i], orientation="horizontal", widths=0.55, showfliers=False,
                   medianprops=dict(color=INK, lw=1.4), boxprops=dict(color=AXIS),
                   whiskerprops=dict(color=AXIS), capprops=dict(color=AXIS))
        ax.scatter(x, i + rng.uniform(-0.18, 0.18, len(x)), s=12, color=CLASS_COL.get(cls[v], MUTED),
                   edgecolor=SURF, linewidth=0.5, zorder=3)
    ax.set_yticks(range(len(labels))); ax.set_yticklabels([lab[v] for v in labels], color=INK)
    ax.axvline(lod, color=MUTED, ls="--", lw=0.8)
    ax.text(lod + 0.5, len(labels) - 0.4, f"{lod:g}% min VAF", fontsize=7, color=MUTED)
    ax.set_xlabel("Variant allele frequency (%)"); ax.set_xlim(0, 100)
    ax.set_title("VAF of the most recurrent variants (each dot = one sample)")
    grid_x(ax)
    save(fig, out, "fig3_vaf_top_variants", figs)


# ------------------------------------------------------------------ fig 4
def summary_plot(flt, out, figs):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 3.6), gridspec_kw=dict(wspace=0.45))
    vc = flt["variant_class"].value_counts().reindex([c for c in CLASS_ORDER if c != "Multi hit"]).dropna()
    vc = vc.iloc[::-1]
    if vc.empty:
        plt.close(fig); return
    a1.barh(range(len(vc)), vc.values, color=[CLASS_COL[c] for c in vc.index], height=0.7)
    a1.set_yticks(range(len(vc))); a1.set_yticklabels(vc.index, color=INK)
    for i, v in enumerate(vc.values):
        a1.text(v + vc.max() * 0.01, i, f"{int(v)}", va="center", fontsize=8, color=INK2)
    a1.set_xlim(0, vc.max() * 1.15)
    a1.set_title("Reportable variants by class"); a1.set_xlabel("Variant calls"); grid_x(a1)

    snv = flt[flt["var_type"] == "SNV"]
    comp = str.maketrans("ACGT", "TGCA")

    def norm(r, a):
        return f"{r}>{a}" if r in "CT" else f"{r.translate(comp)}>{a.translate(comp)}"
    sub = snv.apply(lambda r: norm(r["ref"], r["alt"]), axis=1).value_counts().reindex(SUBS, fill_value=0)
    frac = 100 * sub / max(sub.sum(), 1)
    a2.bar(SUBS, frac.values, color=BLUE, width=0.65)
    for i, v in enumerate(frac.values):
        a2.text(i, v + 0.8, f"{v:.0f}%", ha="center", fontsize=8, color=INK2)
    a2.set_ylabel("% of SNVs"); a2.set_title("SNV substitution spectrum")
    a2.grid(axis="y", color=GRID, lw=0.6); a2.set_axisbelow(True)
    a2.text(0, -0.2, "High C>T can indicate FFPE deamination; high C>A, oxidation / tobacco signature",
            transform=a2.transAxes, fontsize=7, color=MUTED)
    save(fig, out, "fig4_variant_summary", figs)


# ------------------------------------------------------------------ fig 5
def tmb_msi_plot(sm, msi_high, out, figs):
    has_tmb = "TMB" in sm and sm["TMB"].notna().any()
    has_msi = "MSI_pct_unstable" in sm and sm["MSI_pct_unstable"].notna().any()
    if not has_tmb:
        return
    fig, axes = plt.subplots(1, 2 if has_msi else 1, figsize=(11 if has_msi else 6, 3.6),
                             gridspec_kw=dict(wspace=0.3))
    axes = np.atleast_1d(axes)
    t = sm["TMB"].dropna()
    bins = np.arange(0, max(t.max(), 25) + 2.5, 2.5)
    axes[0].hist(t, bins=bins, color=BLUE, edgecolor=SURF, linewidth=1.5)
    for x, lbl in [(10, "10"), (20, "20")]:
        axes[0].axvline(x, color=MUTED, ls="--", lw=0.8)
    hi = (t >= 10).sum()
    axes[0].set_title("TMB distribution")
    axes[0].text(0.98, 0.95, f"median {t.median():.1f} mut/Mb\n{hi} samples >= 10", transform=axes[0].transAxes, ha="right", va="top", fontsize=8, color=INK2)
    axes[0].set_xlabel("TMB (mutations / Mb)"); axes[0].set_ylabel("Samples")
    axes[0].grid(axis="y", color=GRID, lw=0.6); axes[0].set_axisbelow(True)
    if has_msi:
        d = sm.dropna(subset=["TMB", "MSI_pct_unstable"])
        if "MSI_status" in d:
            d = d[d["MSI_status"] != "Insufficient sites"]
        high = d["MSI_pct_unstable"] >= msi_high
        axes[1].scatter(d.loc[~high, "TMB"], d.loc[~high, "MSI_pct_unstable"], s=22, color=BLUE,
                        edgecolor=SURF, lw=0.6, label=f"MSS ({(~high).sum()})")
        axes[1].scatter(d.loc[high, "TMB"], d.loc[high, "MSI_pct_unstable"], s=26, color=RED,
                        edgecolor=SURF, lw=0.6, marker="D", label=f"MSI-High ({high.sum()})")
        axes[1].axhline(msi_high, color=MUTED, ls="--", lw=0.8)
        axes[1].set_xlabel("TMB (mutations / Mb)"); axes[1].set_ylabel("MSI: % unstable sites")
        axes[1].set_title("MSI vs TMB"); axes[1].legend(frameon=False, fontsize=8)
        axes[1].grid(color=GRID, lw=0.6); axes[1].set_axisbelow(True)
    save(fig, out, "fig5_tmb_msi", figs)


# ------------------------------------------------------------------ fig 6
def lollipops(flt, genes, out, figs):
    for g in genes:
        d = flt[(flt["gene"] == g) & flt["protein_change"].ne("")].copy()
        d["aa_pos"] = d["protein_change"].str.extract(r"[A-Za-z*]{1,3}(\d+)")[0].astype(float)
        d = d.dropna(subset=["aa_pos"])
        if d.empty:
            continue
        agg = (d.groupby(["aa_pos", "protein_change", "variant_class"])["sample_id"].nunique()
                 .reset_index(name="n").sort_values("n", ascending=False))
        fig, ax = plt.subplots(figsize=(10, 3))
        L = agg["aa_pos"].max() * 1.08 + 10
        ax.add_patch(plt.Rectangle((0, -0.35), L, 0.7, color=EMPTY, zorder=1))
        for _, r in agg.iterrows():
            ax.plot([r.aa_pos, r.aa_pos], [0.35, r.n], color=AXIS, lw=1, zorder=2)
            ax.scatter(r.aa_pos, r.n, s=40 + 12 * r.n, color=CLASS_COL.get(r.variant_class, MUTED),
                       edgecolor=SURF, lw=0.8, zorder=3)
        for _, r in agg.head(5).iterrows():
            ax.annotate(r.protein_change, (r.aa_pos, r.n), xytext=(0, 9), textcoords="offset points",
                        ha="center", fontsize=7.5, color=INK)
        ax.set_xlim(0, L); ax.set_ylim(-0.5, agg["n"].max() * 1.3 + 0.5)
        ax.set_xlabel(f"{g} protein position (aa)"); ax.set_ylabel("Samples")
        ax.set_title(f"{g}: mutation positions ({d['sample_id'].nunique()} mutated samples)")
        present = [c for c in CLASS_ORDER if c in set(agg["variant_class"])]
        ax.legend(handles=[Patch(color=CLASS_COL[c], label=c) for c in present],
                  frameon=False, fontsize=7.5, loc="upper right")
        ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
        save(fig, out, f"fig6_lollipop_{re.sub(r'[^A-Za-z0-9]', '_', g)}", figs)


# ------------------------------------------------------------------ fig 7
def cooccurrence(mat, n, out, figs):
    genes = list(mat.index[:n])
    if len(genes) < 3:
        return
    b = (mat.loc[genes] != "").astype(int)
    k = len(genes); S = np.zeros((k, k)); P = np.ones((k, k))
    for i in range(k):
        for j in range(i + 1, k):
            a, c = b.iloc[i].values, b.iloc[j].values
            t = [[np.sum(a & c), np.sum(a & (1 - c))], [np.sum((1 - a) & c), np.sum((1 - a) & (1 - c))]]
            orr, p = fisher_exact(t)
            sign = 1 if (t[0][0] * t[1][1]) >= (t[0][1] * t[1][0]) else -1
            S[i, j] = S[j, i] = sign * min(-np.log10(max(p, 1e-10)), 6); P[i, j] = P[j, i] = p
    np.fill_diagonal(S, np.nan)
    cmap = LinearSegmentedColormap.from_list("div", [RED, "#f0efec", BLUE])
    fig, ax = plt.subplots(figsize=(0.42 * k + 2.5, 0.42 * k + 1.5))
    im = ax.imshow(np.ma.masked_invalid(np.tril(S, -1) + np.triu(np.full_like(S, np.nan))),
                   cmap=cmap, vmin=-3, vmax=3)
    for i in range(k):
        for j in range(i):
            if P[i, j] < 0.05:
                ax.text(j, i, "*" if P[i, j] >= 0.01 else "**", ha="center", va="center", fontsize=8, color=INK)
    ax.set_xticks(range(k)); ax.set_xticklabels(genes, rotation=90, fontstyle="italic")
    ax.set_yticks(range(k)); ax.set_yticklabels(genes, fontstyle="italic")
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, shrink=0.6)
    cb.set_label("signed -log10(p)  (blue co-occur / red exclusive)", color=INK2)
    ax.set_title("Pairwise co-occurrence (Fisher; * p<0.05, ** p<0.01)\nhypermutated samples inflate co-occurrence - consider excluding TMB-high", fontsize=10)
    save(fig, out, "fig7_cooccurrence", figs)


# ------------------------------------------------------------------ HTML
def table_html(df, tid):
    head = "".join(f"<th onclick=\"sortT('{tid}',{i})\">{html.escape(str(c))}</th>" for i, c in enumerate(df.columns))
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(v))}</td>" for v in r) + "</tr>"
                   for r in df.astype(str).values)
    return f"<table id='{tid}'><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def write_html(out, figs, tv, tg, sm, fs, n_samples, params, extra=""):
    tiles = [("Samples", n_samples), ("Reportable variant calls", int(tv["n_samples"].sum()) if len(tv) else 0),
             ("Distinct variants", len(tv)), ("Mutated genes", len(tg))]
    if "TMB" in sm and sm["TMB"].notna().any():
        tiles.append(("Median TMB", f"{sm['TMB'].median():.1f}"))
    if "MSI_status" in sm:
        tiles.append(("MSI-High", int((sm["MSI_status"] == "MSI-High").sum())))
    tiles_h = "".join(f"<div class='tile'><div class='v'>{v}</div><div class='k'>{k}</div></div>" for k, v in tiles)
    figs_h = "".join(f"<figure><img src='data:image/png;base64,{b}' alt='{n}'><figcaption>{n}</figcaption></figure>"
                     for n, b in figs)
    tvc = tv.drop(columns=["samples"], errors="ignore").head(200)
    css = """
:root{--bg:#f9f9f7;--surf:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--mut:#898781;--line:#e1e0d9}
@media (prefers-color-scheme:dark){:root{--bg:#0d0d0d;--surf:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--line:#2c2c2a}
 figure img{background:#fcfcfb}}
body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.5 system-ui,-apple-system,'Segoe UI',sans-serif}
main{max-width:1200px;margin:auto;padding:24px 16px}
h1{font-size:24px;margin:0 0 4px}h2{font-size:17px;margin:32px 0 10px}.sub{color:var(--ink2)}
.tiles{display:flex;flex-wrap:wrap;gap:10px;margin:18px 0}
.tile{background:var(--surf);border:1px solid var(--line);border-radius:10px;padding:12px 16px;min-width:130px}
.tile .v{font-size:24px;font-weight:600}.tile .k{color:var(--ink2);font-size:12px}
figure{background:var(--surf);border:1px solid var(--line);border-radius:10px;margin:0 0 14px;padding:12px;overflow-x:auto}
figure img{max-width:100%;display:block}figcaption{color:var(--mut);font-size:12px;margin-top:6px}
.tw{overflow:auto;max-height:520px;border:1px solid var(--line);border-radius:10px;background:var(--surf)}
table{border-collapse:collapse;width:100%;font-size:12.5px;font-variant-numeric:tabular-nums}
th{position:sticky;top:0;background:var(--surf);cursor:pointer;text-align:left;border-bottom:1px solid var(--line);padding:6px 8px;white-space:nowrap}
td{padding:5px 8px;border-bottom:1px solid var(--line);white-space:nowrap}
input{padding:6px 10px;border:1px solid var(--line);border-radius:8px;background:var(--surf);color:var(--ink);margin-bottom:8px;width:260px;max-width:100%}
code{font-size:12px}"""
    js = """
function sortT(id,c){const t=document.getElementById(id),b=t.tBodies[0],r=[...b.rows];
const asc=t.dataset.c==c&&t.dataset.d!='asc';t.dataset.c=c;t.dataset.d=asc?'asc':'desc';
r.sort((x,y)=>{let a=x.cells[c].innerText,z=y.cells[c].innerText,na=parseFloat(a),nz=parseFloat(z);
let v=(!isNaN(na)&&!isNaN(nz))?na-nz:a.localeCompare(z);return asc?v:-v});r.forEach(e=>b.appendChild(e));}
function filt(id,q){q=q.toLowerCase();for(const r of document.getElementById(id).tBodies[0].rows)
r.style.display=r.innerText.toLowerCase().includes(q)?'':'none';}"""
    page = f"""<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>
<title>MutaXplorer Report</title><style>{css}</style></head><body><main>
<h1>{html.escape(params.title)}</h1>
<div class='sub'>{html.escape(params.subtitle)}<br>Reportable = PASS, depth &ge; {params.min_dp}, alt reads &ge; {params.min_alt},
VAF &ge; {params.min_vaf:.0%}, population AF &le; {params.max_pop_af}, protein-altering or splice. Research use only.</div>
<div class='tiles'>{tiles_h}</div>
<h2>SNVs and indels</h2>{figs_h}
<h2>Top variants</h2><input placeholder='Filter (gene, change...)' oninput="filt('tv',this.value)">
<div class='tw'>{table_html(tvc, 'tv')}</div>
<h2>Top genes</h2><div class='tw'>{table_html(tg.head(100), 'tg')}</div>
<h2>Sample metrics</h2><input placeholder='Filter samples' oninput="filt('sm',this.value)">
<div class='tw'>{table_html(sm, 'sm')}</div>
<h2>Filter summary</h2><div class='tw'>{table_html(fs, 'fs')}</div>
{extra}</main><script>{js}</script></body></html>"""
    with open(os.path.join(out, "mutaxplorer_report.html"), "w") as fh:
        fh.write(page)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--mode", choices=["somatic", "germline"], default="somatic")
    ap.add_argument("--title", default=None); ap.add_argument("--subtitle", default="")
    ap.add_argument("--top-genes", type=int, default=30); ap.add_argument("--top-variants", type=int, default=25)
    ap.add_argument("--lollipop-genes", type=int, default=3); ap.add_argument("--msi-high", type=float, default=20.0)
    ap.add_argument("--min-dp", type=int, default=100); ap.add_argument("--min-alt", type=int, default=5)
    ap.add_argument("--min-vaf", type=float, default=0.05); ap.add_argument("--max-pop-af", type=float, default=0.01)
    ap.add_argument("--events", help="oncoplot_events.tsv (CNV / fusion events)")
    ap.add_argument("--extra-html", help="HTML sections to append (CNV, SV, fusions, signatures, identity)")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    if a.title is None:
        a.title = "MutaXplorer - " + ("somatic variant report" if a.mode == "somatic" else "germline variant report")
    def rd(f, **k):
        p = os.path.join(a.cohort, f)
        return pd.read_csv(p, sep="\t", **k) if os.path.exists(p) and os.path.getsize(p) else pd.DataFrame()
    flt = rd("variants_filtered.tsv", keep_default_na=False)
    flt = flt if not flt.empty else pd.DataFrame(columns=["sample_id", "gene", "variant_id", "vaf", "variant_class", "var_type", "ref", "alt", "protein_change"])
    flt["vaf"] = flt["vaf"].astype(float)
    tv, tg = rd("top_variants.tsv", keep_default_na=False), rd("top_genes.tsv")
    tv = tv if not tv.empty else pd.DataFrame(columns=["variant_id", "label", "variant_class", "n_samples", "pct_samples"])
    tg = tg if not tg.empty else pd.DataFrame(columns=["gene", "n_samples"])
    mat = rd("gene_sample_matrix.tsv", index_col=0, keep_default_na=False)
    sm, fs = rd("sample_metrics.tsv"), rd("filter_summary.tsv")
    events = pd.read_csv(a.events, sep="\t") if a.events and os.path.exists(a.events) and os.path.getsize(a.events) else None
    ns = len(sm); figs = []
    if flt.empty:
        print("[warn] no reportable variants (is VEP annotation configured? see filter_summary.tsv) "
              "- only QC figures will be drawn")
    plots = [("oncoplot", lambda: oncoplot(mat, sm, a.top_genes, a.out, figs, events)),
             ("top variants", lambda: top_variants_bar(tv, a.top_variants, ns, a.out, figs)),
             ("VAF", lambda: vaf_plot(flt, tv, min(a.top_variants, 20), a.out, figs, 100 * a.min_vaf)),
             ("variant summary", lambda: summary_plot(flt, a.out, figs)),
             ("TMB/MSI", lambda: a.mode == "somatic" and tmb_msi_plot(sm, a.msi_high, a.out, figs)),
             ("lollipops", lambda: lollipops(flt, list(tg["gene"].head(a.lollipop_genes)), a.out, figs)),
             ("co-occurrence", lambda: cooccurrence(mat, min(20, a.top_genes), a.out, figs))]
    for name, fn in plots:
        try:
            fn()
        except Exception as e:  # one bad figure must never block the report
            plt.close("all")
            print(f"[warn] figure '{name}' skipped: {type(e).__name__}: {e}")
    extra = open(a.extra_html).read() if a.extra_html and os.path.exists(a.extra_html) else ""
    write_html(a.out, figs, tv, tg, sm, fs, ns, a, extra)
    print(f"{len(figs)} figures + mutaxplorer_report.html -> {a.out}")


if __name__ == "__main__":
    main()
