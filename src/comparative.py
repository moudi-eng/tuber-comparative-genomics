"""
Comparative genomics from the literature-sourced dataset (report Sections 4-7):
Table 1-2, genome size, size vs gene number, TE content, CAZyme range, published
cladogram, plus an optional cross-check of genome sizes against NCBI Assembly.
"""
import io

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo

from config import short

plt.rcParams.update({"figure.dpi": 150, "font.size": 10})


def species_table(sources, core_only=True):
    df = pd.DataFrame(sources["species"])
    return df[df["core"]].reset_index(drop=True) if core_only else df


def write_tables(sources, out):
    t1 = species_table(sources)[["name", "genome_size_mb", "genes", "te_percent", "source"]]
    t1.columns = ["Species", "Size (Mb)", "Predicted genes", "TE (%)", "Source"]
    t1.astype(object).where(t1.notna(), "n.s.").to_csv(out / "table1_genomic_features.csv", index=False)

    t2 = pd.DataFrame(sources["divergence_ages"])
    t2["Estimated age (Ma)"] = [f"{a} ± {e}" if e else f"~{a}" for a, e in zip(t2["age_ma"], t2["error_ma"])]
    t2[["event", "Estimated age (Ma)", "source"]].rename(
        columns={"event": "Divergence event", "source": "Source"}).to_csv(out / "table2_divergence_ages.csv", index=False)
    return t1, t2


def _bar(ax, df, col, unit, title, ylabel):
    df = df.dropna(subset=[col]).sort_values(col)
    bars = ax.bar([n for n in df["name"]], df[col], color=df["color"], edgecolor="black", linewidth=0.8)
    for b, v in zip(bars, df[col]):
        ax.text(b.get_x() + b.get_width() / 2, v + df[col].max() * 0.01,
                f"{v:g}{unit}", ha="center", va="bottom", fontsize=9)
    ax.set_ylim(0, df[col].max() * 1.15)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["name"], rotation=15, ha="right", style="italic")


def fig_genome_size(sources, path):
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    _bar(ax, species_table(sources), "genome_size_mb", " Mb",
         "Genome size in four Tuber species", "Assembled genome size (Mb)")
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def fig_te_content(sources, path):
    fig, ax = plt.subplots(figsize=(5.8, 4.2))
    _bar(ax, species_table(sources), "te_percent", "%",
         "Transposable element content", "Transposable elements (% of genome)")
    ax.set_ylim(0, 100)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def pearson(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    return float(np.corrcoef(x, y)[0, 1]) if len(x) >= 3 else float("nan")


def fig_size_vs_genes(df, path, title, published=None):
    """Scatter of genome size vs predicted genes for species with both values."""
    d = df.dropna(subset=["genome_size_mb", "genes"])
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    for _, r in d.iterrows():
        ax.scatter(r["genome_size_mb"], r["genes"], s=110, color=r["color"], edgecolor="black", zorder=3)
        label = short(r["name"]) + (" *" if r.get("size_from_ncbi") else "")
        ax.annotate(label, (r["genome_size_mb"], r["genes"]), xytext=(7, 5),
                    textcoords="offset points", style="italic", fontsize=9)
    if len(d) >= 2:
        k, b = np.polyfit(d["genome_size_mb"], d["genes"], 1)
        xs = np.linspace(d["genome_size_mb"].min() - 5, d["genome_size_mb"].max() + 5, 50)
        ax.plot(xs, k * xs + b, "--", color="gray", linewidth=1)
    r = pearson(d["genome_size_mb"], d["genes"])
    txt = f"this dataset: r = {r:.2f} (n = {len(d)})"
    if published:
        txt += f"\nMurat 2018: r = {published['r']} (n = {published['n']} Pezizomycetes)"
    ax.text(0.98, 0.03, txt, transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            bbox=dict(boxstyle="round", fc="white", ec="lightgray"))
    ax.margins(x=0.18, y=0.15)
    if "size_from_ncbi" in d and d["size_from_ncbi"].any():
        fig.text(0.01, 0.01, "* genome size from NCBI Assembly", fontsize=7, color="gray")
    ax.set_xlabel("Genome size (Mb)")
    ax.set_ylabel("Predicted protein-coding genes")
    ax.set_title(title)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)
    return r, len(d)


def fig_cazyme_range(sources, path):
    c = sources["cazyme_pcwde_range"]
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    ax.hlines(0, c["min"]["count"], c["max"]["count"], color="lightgray", linewidth=8, zorder=1)
    for side, color in (("min", "#7FB37F"), ("max", "#B5523B")):
        e = c[side]
        ax.scatter(e["count"], 0, s=160, color=color, edgecolor="black", zorder=3)
        ax.annotate(f"{e['species']}\n{e['lifestyle']}\n{e['count']} PCWDEs", (e["count"], 0),
                    xytext=(0, 18), textcoords="offset points", ha="center", fontsize=8)
    ax.annotate("Tuber spp.: lower end\n(no individual value in source)", (c["min"]["count"] + 15, 0),
                xytext=(0, -42), textcoords="offset points", ha="center", fontsize=8, color="gray")
    ax.set_xlim(c["min"]["count"] - 30, c["max"]["count"] + 30)
    ax.set_ylim(-1, 1)
    ax.set_yticks([])
    ax.set_xlabel("Secreted plant cell wall-degrading CAZymes")
    ax.set_title(f"Range across 8 Pezizomycetes ({c['source']})")
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def fig_published_cladogram(sources, path):
    """Cladogram redrawn from the published topology and ages (not reconstructed from sequences)."""
    pc = sources["published_cladogram"]
    tree = Phylo.read(io.StringIO(pc["newick"]), "newick")
    colors = {s["name"].replace(" ", "_"): s["color"] for s in sources["species"]}
    root_age = pc["root_age_ma"]

    depths = tree.depths()
    tips = tree.get_terminals()
    ypos = {t: i for i, t in enumerate(reversed(tips))}

    def y_of(clade):
        if clade.is_terminal():
            return ypos[clade]
        ys = [y_of(c) for c in clade.clades]
        ypos[clade] = (min(ys) + max(ys)) / 2
        return ypos[clade]

    y_of(tree.root)
    x = {c: root_age - depths[c] for c in depths}  # Ma before present

    fig, ax = plt.subplots(figsize=(8.5, 5))
    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal():
            ax.scatter(x[clade], ypos[clade], s=90, color=colors.get(clade.name, "white"),
                       edgecolor="black", zorder=3)
            ax.text(x[clade] - 3, ypos[clade], short(clade.name), va="center", style="italic", fontsize=11)
            continue
        ys = [ypos[c] for c in clade.clades]
        ax.plot([x[clade], x[clade]], [min(ys), max(ys)], color="#3B2F2F", linewidth=2)
        for c in clade.clades:
            ax.plot([x[clade], x[c]], [ypos[c], ypos[c]], color="#3B2F2F", linewidth=2)
        label = clade.name or ""
        if label in pc.get("schematic_nodes", []):
            ax.text(x[clade] + 1, ypos[clade] - 0.35, "schematic node\n(no precise age found\nin the cited sources)",
                    ha="right", va="top", fontsize=8, style="italic", color="gray")
        elif label:
            ax.text(x[clade] + 1, ypos[clade] + 0.1, label, ha="right", va="bottom", fontsize=9, color="#555")
    ax.set_xlim(root_age + 20, -25)
    ax.set_ylim(-0.8, len(tips) - 0.2)
    ax.set_yticks([])
    ax.set_xlabel("Million years before present")
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.set_title(f"Schematic cladogram of the genus Tuber\n(topology and ages after {pc['source']})")
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def ncbi_crosscheck(sources, assembly_stats, out):
    """Compare literature genome sizes with NCBI Assembly and fill sizes the report left empty."""
    df = species_table(sources, core_only=False).copy()
    df["ncbi_accession"] = [(assembly_stats.get(n) or {}).get("accession") for n in df["name"]]
    df["ncbi_size_mb"] = [(assembly_stats.get(n) or {}).get("size_mb") for n in df["name"]]
    df["ncbi_level"] = [(assembly_stats.get(n) or {}).get("level") for n in df["name"]]
    df["ncbi_contig_n50"] = [(assembly_stats.get(n) or {}).get("contig_n50") for n in df["name"]]
    df["difference_pct"] = ((df["ncbi_size_mb"] - df["genome_size_mb"]) / df["genome_size_mb"] * 100).round(2)
    df[["name", "genome_size_mb", "source", "ncbi_accession", "ncbi_size_mb", "difference_pct",
        "ncbi_level", "ncbi_contig_n50"]].to_csv(out / "ncbi_assembly_crosscheck.csv", index=False)

    df["size_from_ncbi"] = df["genome_size_mb"].isna() & df["ncbi_size_mb"].notna()
    df["genome_size_mb"] = df["genome_size_mb"].fillna(df["ncbi_size_mb"])
    return df
