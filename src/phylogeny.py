"""
De novo phylogeny from marker sequences: MAFFT alignment -> gap-column trimming ->
IQ-TREE 2 (ModelFinder + ultrafast bootstrap). Falls back to FastTree, and if no
external tool is installed, to a neighbour-joining tree from pairwise alignments
(Biopython only) so the pipeline always produces a tree.
"""
import io
import shutil
import subprocess
from itertools import combinations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import AlignIO, Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from config import short


def _which(*names):
    for n in names:
        if shutil.which(n):
            return n
    return None


def _run(cmd, **kw):
    print("  $", " ".join(map(str, cmd)))
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw)


def split_id(seq_id):
    """'Tuber_melanosporum__OP978675.1' -> ('Tuber_melanosporum', 'OP978675.1'); '|' also accepted."""
    for sep in ("__", "|"):
        if sep in seq_id:
            taxon, acc = seq_id.split(sep, 1)
            return taxon, acc
    return seq_id, ""


def taxon_of(seq_id):
    return split_id(seq_id)[0]


# ----------------------------------------------------------------------
# Alignment
# ----------------------------------------------------------------------
def align(fasta, out_dir):
    mafft = _which("mafft", "mafft.bat")
    if not mafft:
        return None
    aln = out_dir / "alignment.fasta"
    res = _run([mafft, "--auto", "--adjustdirection", "--quiet", str(fasta)])
    # --adjustdirection prefixes reverse-complemented IDs with "_R_"
    aln.write_text(res.stdout.replace(">_R_", ">"))
    return aln


def trim(aln_path, out_dir, gap_threshold):
    aln = AlignIO.read(aln_path, "fasta")
    arr = np.array([list(str(r.seq)) for r in aln])
    keep = (arr == "-").mean(axis=0) <= gap_threshold
    trimmed = MultipleSeqAlignment([SeqRecord(Seq("".join(row)), id=rec.id, description="")
                                    for rec, row in zip(aln, arr[:, keep])])
    path = out_dir / "alignment_trimmed.fasta"
    AlignIO.write(trimmed, path, "fasta")
    return path, {"columns_before": int(arr.shape[1]), "columns_after": int(keep.sum())}


def p_distance_matrix(aln_path):
    """Uncorrected p-distance between aligned sequences (gaps ignored pairwise)."""
    aln = AlignIO.read(aln_path, "fasta")
    names = [r.id for r in aln]
    arr = np.array([list(str(r.seq).upper()) for r in aln])
    n = len(names)
    d = np.zeros((n, n))
    for i, j in combinations(range(n), 2):
        ok = (arr[i] != "-") & (arr[j] != "-")
        d[i, j] = d[j, i] = (arr[i][ok] != arr[j][ok]).mean() if ok.any() else np.nan
    return pd.DataFrame(d, index=names, columns=names)


def pairwise_distance_matrix(fasta):
    """Fallback without MAFFT: p-distance from global pairwise alignments (both strands tried)."""
    recs = list(SeqIO.parse(fasta, "fasta"))
    al = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                         open_gap_score=-4, extend_gap_score=-1)
    names = [r.id for r in recs]
    n = len(recs)
    d = np.zeros((n, n))
    for i, j in combinations(range(n), 2):
        a, b = recs[i].seq, recs[j].seq
        best = max((al.align(a, b)[0], al.align(a, b.reverse_complement())[0]), key=lambda x: x.score)
        s1, s2 = best[0], best[1]
        pairs = [(x, y) for x, y in zip(s1, s2) if x != "-" and y != "-"]
        d[i, j] = d[j, i] = sum(x != y for x, y in pairs) / max(len(pairs), 1)
    return pd.DataFrame(d, index=names, columns=names)


# ----------------------------------------------------------------------
# Tree inference
# ----------------------------------------------------------------------
def infer_tree(aln_path, out_dir, bootstrap, seed=42):
    iqtree = _which("iqtree2", "iqtree3", "iqtree")
    if iqtree:
        prefix = out_dir / "iqtree"
        version = _run([iqtree, "--version"]).stdout
        if iqtree != "iqtree" or "version 2" in version or "version 3" in version:
            cmd = [iqtree, "-s", aln_path, "-m", "MFP", "-B", bootstrap, "-T", "AUTO",
                   "--seed", seed, "--prefix", prefix, "-redo", "--quiet"]
        else:  # IQ-TREE 1.x syntax
            cmd = [iqtree, "-s", aln_path, "-m", "MFP", "-bb", bootstrap, "-nt", "AUTO",
                   "-seed", seed, "-pre", prefix, "-redo", "-quiet"]
        _run([str(c) for c in cmd])
        model = next((l.split(":", 1)[1].strip() for l in (prefix.with_suffix(".iqtree")).read_text().splitlines()
                      if l.startswith("Best-fit model")), "see iqtree.iqtree")
        return Phylo.read(f"{prefix}.treefile", "newick"), f"IQ-TREE ML ({model}, {bootstrap} UFBoot)"

    fasttree = _which("FastTree", "fasttree", "FastTreeMP")
    if fasttree:
        res = _run([fasttree, "-nt", "-gtr", "-quiet", str(aln_path)])
        (out_dir / "fasttree.nwk").write_text(res.stdout)
        return Phylo.read(io.StringIO(res.stdout), "newick"), "FastTree ML (GTR, SH-like support)"
    return None, None


def nj_tree(dist):
    names = list(dist.index)
    lower = [[float(dist.iloc[i, j]) for j in range(i + 1)] for i in range(len(names))]
    tree = DistanceTreeConstructor().nj(DistanceMatrix(names, lower))
    for c in tree.get_nonterminals():
        c.name = None
    return tree


def root_tree(tree, outgroup_taxa):
    out_tips = [t for t in tree.get_terminals() if taxon_of(t.name) in outgroup_taxa]
    if out_tips:
        try:
            tree.root_with_outgroup(out_tips if len(out_tips) > 1 else out_tips[0])
        except ValueError:  # outgroup not monophyletic: root on its first sequence
            tree.root_with_outgroup(out_tips[0])
    else:
        tree.root_at_midpoint()
    tree.ladderize()
    return tree


# ----------------------------------------------------------------------
# Interpretation
# ----------------------------------------------------------------------
def species_monophyly(tree):
    rows = []
    for taxon in sorted({taxon_of(t.name) for t in tree.get_terminals()}):
        tips = [t for t in tree.get_terminals() if taxon_of(t.name) == taxon]
        mono = len(tips) == 1 or tree.is_monophyletic(tips) is not False
        support = None
        if len(tips) > 1 and mono:
            support = tree.common_ancestor(tips).confidence
        rows.append({"taxon": taxon, "n_sequences": len(tips), "monophyletic": bool(mono), "support": support})
    return pd.DataFrame(rows)


def clade_support(tree, taxa):
    tips = [t for t in tree.get_terminals() if taxon_of(t.name) in taxa]
    if not tips:
        return None, None
    mono = tree.is_monophyletic(tips) is not False
    return bool(mono), (tree.common_ancestor(tips).confidence if mono else None)


def closest_species(dist, focal):
    """Mean p-distance from the focal species to every other taxon."""
    taxa = dist.index.map(taxon_of)
    mean = {}
    for t in sorted(set(taxa) - {focal}):
        mean[t] = float(dist.loc[taxa == focal, taxa == t].values.mean())
    return pd.Series(mean).sort_values()


def _splits(tree, taxa):
    """Non-trivial bipartitions of a species-level tree, as frozensets not containing the first taxon."""
    taxa = sorted(taxa)
    anchor, full = taxa[0], frozenset(taxa)
    out = set()
    for c in tree.get_nonterminals():
        s = frozenset(taxon_of(t.name) for t in c.get_terminals()) & full
        if anchor in s:
            s = full - s
        if 1 < len(s) < len(full) - 1:
            out.add(s)
    return out


def species_tree(tree):
    """One tip per taxon (the first sequence of each), for topology comparison."""
    t = Phylo.read(io.StringIO(tree.format("newick")), "newick")
    seen = set()
    for tipc in list(t.get_terminals()):
        tax = taxon_of(tipc.name)
        if tax in seen:
            t.prune(tipc)
        else:
            seen.add(tax)
            tipc.name = tax
    return t


def compare_topology(tree, reference_newick):
    """Robinson-Foulds distance between the inferred species tree and the published cladogram."""
    ref = Phylo.read(io.StringIO(reference_newick), "newick")
    sp = species_tree(tree)
    shared = {t.name for t in ref.get_terminals()} & {t.name for t in sp.get_terminals()}
    if len(shared) < 4:
        return None
    for tr in (ref, sp):
        for t in list(tr.get_terminals()):
            if t.name not in shared:
                tr.prune(t)
    a, b = _splits(ref, shared), _splits(sp, shared)
    max_rf = 2 * (len(shared) - 3)
    return {"shared_taxa": sorted(shared), "rf_distance": len(a ^ b), "max_rf": max_rf,
            "reference_splits_recovered": [sorted(s) for s in a & b],
            "reference_splits_missing": [sorted(s) for s in a - b]}


# ----------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------
def _readable(color):
    """Darken very light colours so tip labels stay legible on white."""
    rgb = matplotlib.colors.to_rgb(color)
    lum = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
    return tuple(c * 0.55 for c in rgb) if lum > 0.75 else rgb


def plot_tree(tree, colors, method, marker, outgroup_taxa, path, min_support=70):
    """`colors` maps taxon ('Tuber_melanosporum') to a colour."""
    colors = {short(k): _readable(v) for k, v in colors.items()}
    n = len(tree.get_terminals())
    fig, ax = plt.subplots(figsize=(9, max(4, 0.32 * n + 1.5)))

    def label(c):
        if c.is_terminal():
            taxon, acc = split_id(c.name)
            return f" {short(taxon)}  {acc}"
        return None

    def support(c):
        conf = c.confidence
        if conf is None or c.is_terminal():
            return None
        conf = conf * 100 if conf <= 1 else conf  # FastTree/NJ report 0-1
        return f"{conf:.0f}" if conf >= min_support else None

    Phylo.draw(tree, axes=ax, do_show=False, label_func=label, branch_labels=support,
               label_colors=lambda lbl: colors.get(lbl.strip().split("  ")[0], "black"))
    for txt in ax.texts:
        if txt.get_text().strip() and not txt.get_text().strip().isdigit():
            txt.set_style("italic")
        else:
            txt.set_fontsize(7); txt.set_color("#555")
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.set_ylabel("")
    ax.set_xlabel("Substitutions per site")
    og = ", ".join(short(o) for o in outgroup_taxa) or "midpoint"
    ax.set_title(f"{marker} phylogeny of Tuber — {method}\n(rooted on {og}; support ≥ {min_support} shown)")
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def plot_distance_heatmap(dist, path):
    taxa = dist.index.map(taxon_of)
    order = sorted(set(taxa))
    mean = pd.DataFrame([[np.nanmean(dist.loc[taxa == a, taxa == b].values) for b in order] for a in order],
                        index=order, columns=order)
    labels = [short(t) for t in order]
    fig, ax = plt.subplots(figsize=(1.0 * len(order) + 2.5, 0.8 * len(order) + 2))
    im = ax.imshow(mean.values * 100, cmap="YlOrBr")
    for i in range(len(order)):
        for j in range(len(order)):
            v = mean.values[i, j] * 100
            ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=8,
                    color="white" if v > np.nanmax(mean.values) * 70 else "black")
    ax.set_xticks(range(len(order))); ax.set_xticklabels(labels, rotation=35, ha="right", style="italic")
    ax.set_yticks(range(len(order))); ax.set_yticklabels(labels, style="italic")
    fig.colorbar(im, ax=ax, fraction=0.046, label="Mean p-distance (%)")
    ax.set_title("Mean pairwise sequence divergence between taxa")
    fig.tight_layout(); fig.savefig(path); plt.close(fig)
    return mean
