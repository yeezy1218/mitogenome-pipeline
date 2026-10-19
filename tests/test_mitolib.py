import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts"))

from mitolib import (  # noqa: E402
    GENES,
    build_supermatrix,
    clean_protein,
    gene_from_header,
    gene_from_product,
    genes_by_taxon,
    normalize_gene,
    partitions_nexus,
    read_fasta,
    seq_stats,
    write_fasta,
)


@pytest.mark.parametrize("name,expected", [
    ("COX1", "COX1"), ("COI", "COX1"), ("co1", "COX1"), ("MT-CO1", "COX1"),
    ("COIII", "COX3"), ("cob", "CYTB"), ("Cyt-b", "CYTB"), ("nad4L", "ND4L"),
    ("ND4", "ND4"), ("NADH6", "ND6"), ("MT-ATP8", "ATP8"), ("ATPase6", "ATP6"),
    ("rrnL", None), ("trnK", None), ("", None), (None, None),
])
def test_normalize_gene(name, expected):
    assert normalize_gene(name) == expected


def test_every_canonical_name_maps_to_itself():
    assert all(normalize_gene(g) == g for g in GENES)


@pytest.mark.parametrize("product,expected", [
    ("cytochrome c oxidase subunit I", "COX1"),
    ("cytochrome c oxidase subunit III", "COX3"),
    ("cytochrome b", "CYTB"),
    ("NADH dehydrogenase subunit 4L", "ND4L"),
    ("NADH dehydrogenase subunit 4", "ND4"),
    ("ATP synthase F0 subunit 8", "ATP8"),
    ("16S ribosomal RNA", None),
])
def test_gene_from_product(product, expected):
    assert gene_from_product(product) == expected


@pytest.mark.parametrize("header,expected", [
    (">Homo_sapiens@COX1", "COX1"),
    ("Sample1@nad4l", "ND4L"),
    ("contig_12 ND5", "ND5"),
    ("Sample1_cob", "CYTB"),
    ("Sample1@rrnS", None),
])
def test_gene_from_header(header, expected):
    assert gene_from_header(header) == expected


def test_clean_protein():
    assert clean_protein("mfa*") == "MFA"
    assert clean_protein("MF*AB-J") == "MFXAXX"


def test_genes_by_taxon_keeps_longest_copy():
    records = [("s@COX1", "MFA"), ("s@COI", "MFAWL*"), ("s@rrnL", "ACGT")]
    assert genes_by_taxon(records) == {"COX1": "MFAWL"}


def test_build_supermatrix_fills_missing_genes_with_gaps():
    alignments = {
        "COX1": {"A": "MF-A", "B": "MFWA"},
        "CYTB": {"B": "LL"},
        "ND6": {},
    }
    matrix, parts = build_supermatrix(alignments, ["A", "B"])
    assert matrix == {"A": "MF-A--", "B": "MFWALL"}
    assert parts == [("COX1", 1, 4), ("CYTB", 5, 6)]
    assert "charset CYTB = 5-6;" in partitions_nexus(parts)


def test_build_supermatrix_rejects_unaligned_input():
    with pytest.raises(ValueError, match="different lengths"):
        build_supermatrix({"COX1": {"A": "MFA", "B": "MF"}}, ["A", "B"])


def test_fasta_round_trip(tmp_path):
    path = tmp_path / "x.fasta"
    records = [("one", "A" * 130), ("two", "MF")]
    write_fasta(path, records)
    assert read_fasta(path) == records


def test_seq_stats():
    stats = seq_stats([("c1", "GGCC"), ("c2", "AATT")])
    assert stats == {"contigs": 2, "length_bp": 8, "gc_pct": 50.0}
    assert seq_stats([])["gc_pct"] is None
