"""
NCBI access with Biopython Entrez: marker sequences (ITS, LSU, RPB2, EF1a) and
genome-assembly statistics. Everything is cached in data/ncbi/ so re-runs are offline.
"""
import json
import re
import time

import pandas as pd
from Bio import Entrez, SeqIO

from config import tip

MARKER_QUERIES = {
    "ITS":  "(internal transcribed spacer[Title] OR ITS1[Title] OR ITS2[Title])",
    "LSU":  "(large subunit ribosomal RNA[Title] OR 28S[Title] OR LSU[Title])",
    "RPB2": "(RPB2[Title] OR RNA polymerase II second largest subunit[Title])",
    "EF1a": "(translation elongation factor 1-alpha[Title] OR tef1[Title] OR EF1a[Title])",
}


def _entrez(cfg):
    Entrez.email = cfg.email
    Entrez.tool = "tuber-comparative-genomics"


def _pause():
    time.sleep(0.4)  # NCBI allows ~3 requests/second without an API key


# ----------------------------------------------------------------------
# Marker sequences
# ----------------------------------------------------------------------
def _species_records(species, cfg, retmax=60):
    query = (f'"{species}"[Organism] AND {cfg.query or MARKER_QUERIES[cfg.marker]} '
             f'AND {cfg.min_length}:{cfg.max_length}[SLEN]')
    ids = Entrez.read(Entrez.esearch(db="nucleotide", term=query, retmax=retmax))["IdList"]
    _pause()
    if not ids:
        return []
    handle = Entrez.efetch(db="nucleotide", id=",".join(ids), rettype="gb", retmode="text")
    records = list(SeqIO.parse(handle, "genbank"))
    _pause()
    # exact species only (no "cf.", "aff.", subspecies or mislabelled organisms)
    return [r for r in records if r.annotations.get("organism", "").strip() == species]


def _pick(records, n, marker):
    """Prefer records spanning the whole marker (e.g. ITS1 + 5.8S + ITS2), then lengths near the median."""
    if not records:
        return []

    def complete(r):
        d = r.description.lower()
        if marker == "ITS":
            return ("internal transcribed spacer 1" in d or "its1" in d) and \
                   ("internal transcribed spacer 2" in d or "its2" in d)
        return True

    pool = [r for r in records if complete(r)] or records
    pool = [r for r in pool if set(str(r.seq).upper()) <= set("ACGTN")] or pool
    median = sorted(len(r) for r in pool)[len(pool) // 2]
    pool.sort(key=lambda r: (abs(len(r) - median), r.id))
    return pool[:n]


def fetch_marker_sequences(cfg, cache_dir):
    """Download `seqs_per_species` sequences per taxon and write one FASTA (header: Genus_species__accession; IQ-TREE would rename "|")."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    fasta = cache_dir / f"{cfg.marker}_sequences.fasta"
    meta_path = cache_dir / f"{cfg.marker}_sequences.csv"
    wanted = set(tip(s) for s in cfg.all_taxa)
    if fasta.exists() and meta_path.exists():
        meta = pd.read_csv(meta_path)
        if wanted <= set(meta["taxon"]):
            print(f"Using cached sequences: {fasta}")
            return fasta, meta

    _entrez(cfg)
    rows, lines = [], []
    for species in cfg.all_taxa:
        recs = _pick(_species_records(species, cfg), cfg.seqs_per_species, cfg.marker)
        print(f"  {species}: {len(recs)} {cfg.marker} sequence(s) {[r.id for r in recs]}")
        for r in recs:
            lines.append(f">{tip(species)}__{r.id}\n{str(r.seq).upper()}\n")
            rows.append({"taxon": tip(species), "species": species, "accession": r.id,
                         "length": len(r), "description": r.description,
                         "outgroup": species in cfg.outgroup})
    if not rows:
        raise RuntimeError("No sequences downloaded; check species names, marker and network access.")
    fasta.write_text("".join(lines))
    meta = pd.DataFrame(rows)
    meta.to_csv(meta_path, index=False)
    return fasta, meta


# ----------------------------------------------------------------------
# Genome assemblies
# ----------------------------------------------------------------------
def fetch_assembly_stats(species_list, cfg, cache_dir):
    """Size, level and contig N50 of each species' NCBI assembly (reference genome preferred)."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache = cache_dir / "assembly_stats.json"
    stats = json.loads(cache.read_text()) if cache.exists() else {}
    missing = [s for s in species_list if s not in stats]
    if missing:
        _entrez(cfg)
        for species in missing:
            ids = Entrez.read(Entrez.esearch(db="assembly", term=f'"{species}"[Organism]', retmax=20))["IdList"]
            _pause()
            if not ids:
                stats[species] = None
                continue
            docs = Entrez.read(Entrez.esummary(db="assembly", id=",".join(ids)), validate=False)
            _pause()
            best = None
            for d in docs["DocumentSummarySet"]["DocumentSummary"]:
                size = re.search(r'category="total_length" sequence_tag="all">(\d+)', d["Meta"])
                entry = {
                    "accession": str(d["AssemblyAccession"]),
                    "assembly_name": str(d["AssemblyName"]),
                    "level": str(d["AssemblyStatus"]),
                    "size_mb": round(int(size.group(1)) / 1e6, 2) if size else None,
                    "contig_n50": int(d.get("ContigN50") or 0) or None,
                    "refseq_category": str(d.get("RefSeq_category", "")),
                    "submitted": str(d.get("SubmissionDate", ""))[:10],
                }
                rank = (entry["refseq_category"] == "reference genome", entry["submitted"])
                if best is None or rank > best[0]:
                    best = (rank, entry)
            stats[species] = best[1] if best else None
            print(f"  {species}: {stats[species]}")
        cache.write_text(json.dumps(stats, indent=2))
    return stats
