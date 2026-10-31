"""
Comparative genomics and phylogeny of the genus Tuber (Perigord black truffle).

Three steps, each can be run on its own:
  comparative  literature-sourced genome features, Tables 1-2, Figures 1-5 (report Sections 4-7)
  ncbi         cross-check genome sizes against NCBI Assembly and fill missing sizes
  phylogeny    de novo marker phylogeny: NCBI download -> MAFFT -> trimming -> IQ-TREE
               (the reconstruction the report lists as future work)

Usage
-----
    python src/pipeline.py --config configs/tuber.json --email you@example.com
    python src/pipeline.py --config configs/tuber.json --steps comparative
    python src/pipeline.py --config configs/tuber.json --marker LSU --name tuber_lsu

External tools for the phylogeny step: MAFFT and IQ-TREE 2 (on Colab / Ubuntu:
`apt-get install mafft iqtree`). Without them the step still runs, using a
neighbour-joining tree from pairwise alignments.
"""
import argparse
import json
import sys
import warnings
from pathlib import Path

import pandas as pd
from Bio import Phylo

sys.path.insert(0, str(Path(__file__).resolve().parent))
import comparative  # noqa: E402
import ncbi  # noqa: E402
import phylogeny as phy  # noqa: E402
from config import Config, load_sources, short, tip  # noqa: E402

STEPS = ("comparative", "ncbi", "phylogeny")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default="configs/tuber.json")
    p.add_argument("--steps", default=",".join(STEPS), help=f"comma-separated subset of {STEPS}")
    p.add_argument("--email", help="contact e-mail sent to NCBI with each request")
    p.add_argument("--marker", help="ITS, LSU, RPB2 or EF1a")
    p.add_argument("--fasta", dest="local_fasta", help="use your own FASTA (headers: Genus_species__id)")
    p.add_argument("--name", help="run name (output goes to results/<name>)")
    return p.parse_args()


def run_comparative(cfg, sources, out, fig):
    print("\n=== Comparative genomics (literature-sourced) ===")
    t1, t2 = comparative.write_tables(sources, out)
    print(t1.to_string(index=False))
    comparative.fig_genome_size(sources, fig / "fig1_genome_size.png")
    corr = next(c for c in sources["published_correlations"] if c["y"] == "genes")
    r, n = comparative.fig_size_vs_genes(comparative.species_table(sources), fig / "fig2_size_vs_genes.png",
                                         "Genome size vs number of predicted genes", corr)
    print(f"Size vs genes in this dataset: r = {r:.2f} (n = {n}); published: r = {corr['r']} (n = {corr['n']})")
    comparative.fig_te_content(sources, fig / "fig3_te_content.png")
    comparative.fig_cazyme_range(sources, fig / "fig4_cazyme_range.png")
    comparative.fig_published_cladogram(sources, fig / "fig5_published_cladogram.png")
    print(f"\nDivergence ages (Table 2):\n{t2[['event', 'age_ma', 'source']].to_string(index=False)}")
    return {"size_vs_genes_r": r, "size_vs_genes_n": n}


def run_ncbi(cfg, sources, out, fig):
    print("\n=== NCBI Assembly cross-check ===")
    names = [s["name"] for s in sources["species"]]
    stats = ncbi.fetch_assembly_stats(names, cfg, cfg.path("data/ncbi"))
    df = comparative.ncbi_crosscheck(sources, stats, out)
    print(pd.read_csv(out / "ncbi_assembly_crosscheck.csv").to_string(index=False))
    corr = next(c for c in sources["published_correlations"] if c["y"] == "genes")
    r, n = comparative.fig_size_vs_genes(df, fig / "fig6_size_vs_genes_extended.png",
                                         "Genome size vs genes, incl. T. brumale and T. indicum", corr)
    print(f"Extended size vs genes: r = {r:.2f} (n = {n})")
    return {"extended_size_vs_genes_r": r, "extended_size_vs_genes_n": n}


def run_phylogeny(cfg, sources, out, fig):
    print(f"\n=== De novo {cfg.marker} phylogeny ===")
    work = out / "phylogeny"
    work.mkdir(parents=True, exist_ok=True)

    if cfg.local_fasta:
        fasta = cfg.path(cfg.local_fasta)
    else:
        fasta, meta = ncbi.fetch_marker_sequences(cfg, cfg.path(f"data/ncbi/{cfg.name}"))
        meta.to_csv(out / "sequences_used.csv", index=False)

    aln = phy.align(fasta, work)
    if aln:
        trimmed, trim_stats = phy.trim(aln, work, cfg.gap_threshold)
        print(f"Alignment: {trim_stats['columns_before']} columns -> {trim_stats['columns_after']} after trimming")
        dist = phy.p_distance_matrix(trimmed)
        tree, method = phy.infer_tree(trimmed, work, cfg.bootstrap)
        if tree is None:
            tree, method = phy.nj_tree(dist), "Neighbour-joining (p-distance; IQ-TREE/FastTree not installed)"
    else:
        print("MAFFT not found: using pairwise alignments + neighbour-joining.")
        trim_stats = None
        dist = phy.pairwise_distance_matrix(fasta)
        tree, method = phy.nj_tree(dist), "Neighbour-joining (pairwise p-distance; MAFFT not installed)"

    outgroup = [tip(o) for o in cfg.outgroup]
    tree = phy.root_tree(tree, outgroup)
    Phylo.write(tree, out / f"{cfg.marker}_tree.nwk", "newick")
    dist.to_csv(work / "p_distance_matrix.csv")

    colors = {tip(s["name"]): s["color"] for s in sources["species"]}
    colors.update({tip(o): "#777777" for o in cfg.outgroup})
    phy.plot_tree(tree, colors, method, cfg.marker, cfg.outgroup, fig / f"fig7_{cfg.marker}_tree.png")
    phy.plot_distance_heatmap(dist, fig / "fig8_distance_heatmap.png")

    mono = phy.species_monophyly(tree)
    mono.to_csv(out / "species_monophyly.csv", index=False)
    tuber_mono, tuber_support = phy.clade_support(tree, [tip(s) for s in cfg.species])
    closest = phy.closest_species(dist, tip(cfg.focal_species))
    topo = phy.compare_topology(tree, sources["published_cladogram"]["newick"])

    print(f"Method: {method}")
    print(f"Sequences: {len(tree.get_terminals())}")
    print(f"\nSpecies monophyly:\n{mono.to_string(index=False)}")
    print(f"\nTuber monophyletic vs outgroup: {tuber_mono} (support: {tuber_support})")
    print(f"\nClosest taxa to {short(cfg.focal_species)} (mean p-distance):")
    print((closest * 100).round(2).astype(str).add(" %").to_string())
    if topo:
        print(f"\nTopology vs published cladogram (shared taxa {topo['shared_taxa']}): "
              f"RF = {topo['rf_distance']} / {topo['max_rf']}")
        print(f"  recovered splits: {topo['reference_splits_recovered']}")
        print(f"  missing splits:   {topo['reference_splits_missing']}")

    return {
        "marker": cfg.marker, "method": method, "n_sequences": len(tree.get_terminals()),
        "alignment": trim_stats,
        "tuber_monophyletic": tuber_mono, "tuber_support": tuber_support,
        "species_monophyletic": dict(zip(mono["taxon"], mono["monophyletic"])),
        "closest_to_focal": {k: round(v, 4) for k, v in closest.items()},
        "topology_vs_published": topo,
    }


def main():
    warnings.filterwarnings("ignore")
    args = parse_args()
    overrides = {k: v for k, v in vars(args).items() if k not in ("config", "steps")}
    cfg = Config.load(args.config, **overrides)
    steps = [s.strip() for s in args.steps.split(",") if s.strip()]
    unknown = set(steps) - set(STEPS)
    if unknown:
        raise SystemExit(f"Unknown steps {unknown}; choose from {STEPS}")
    if {"ncbi", "phylogeny"} & set(steps) and not cfg.local_fasta and cfg.email.endswith("example.com"):
        print("Note: set your e-mail with --email (NCBI asks for a contact address).")

    sources = load_sources(cfg)
    out, fig = cfg.out_dir, cfg.out_dir / "figures"
    fig.mkdir(parents=True, exist_ok=True)

    summary = {"name": cfg.name, "steps": steps}
    if "comparative" in steps:
        summary["comparative"] = run_comparative(cfg, sources, out, fig)
    if "ncbi" in steps:
        summary["ncbi"] = run_ncbi(cfg, sources, out, fig)
    if "phylogeny" in steps:
        summary["phylogeny"] = run_phylogeny(cfg, sources, out, fig)

    with open(out / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\nAll tables and figures saved to: {out}")


if __name__ == "__main__":
    main()
