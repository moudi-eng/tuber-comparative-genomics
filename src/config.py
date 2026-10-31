"""Run configuration (configs/*.json) and the literature-sourced dataset (data/verified_sources.json)."""
from dataclasses import dataclass, field, fields
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]


@dataclass
class Config:
    name: str = "tuber"
    email: str = "your.name@example.com"   # NCBI asks for a contact e-mail with every Entrez request

    species: list = field(default_factory=list)
    outgroup: list = field(default_factory=list)

    marker: str = "ITS"                    # ITS, LSU, RPB2 or EF1a (see ncbi.MARKER_QUERIES)
    query: str | None = None               # custom Entrez title query, overrides `marker`
    seqs_per_species: int = 3
    min_length: int = 450
    max_length: int = 1200
    local_fasta: str | None = None         # use your own sequences instead of downloading

    gap_threshold: float = 0.5             # drop alignment columns with more gaps than this
    bootstrap: int = 1000                  # IQ-TREE ultrafast bootstrap replicates
    focal_species: str = "Tuber melanosporum"

    sources_file: str = "data/verified_sources.json"
    output_dir: str | None = None          # default: results/<name>

    def path(self, value):
        if value is None:
            return None
        p = Path(value)
        return p if p.is_absolute() else ROOT / p

    @property
    def out_dir(self):
        return self.path(self.output_dir or f"results/{self.name}")

    @property
    def all_taxa(self):
        return list(self.species) + list(self.outgroup)

    @classmethod
    def load(cls, path, **overrides):
        with open(path, encoding="utf-8") as f:
            raw = json.load(f)
        raw.update({k: v for k, v in overrides.items() if v is not None})
        unknown = set(raw) - {f.name for f in fields(cls)}
        if unknown:
            raise ValueError(f"Unknown config keys: {sorted(unknown)}")
        return cls(**raw)


def load_sources(cfg):
    with open(cfg.path(cfg.sources_file), encoding="utf-8") as f:
        return json.load(f)


def tip(name):
    """Species name as used in FASTA headers and tree tips: 'Tuber melanosporum' -> 'Tuber_melanosporum'."""
    return name.strip().replace(" ", "_")


def short(name):
    """'Tuber melanosporum' -> 'T. melanosporum'."""
    parts = name.replace("_", " ").split()
    return f"{parts[0][0]}. {' '.join(parts[1:])}" if len(parts) > 1 else name
