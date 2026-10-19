"""Pure-Python helpers shared by the pipeline scripts (and unit-tested in tests/)."""

from __future__ import annotations

import re
from pathlib import Path

GENES = ["ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB",
         "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"]

# Names seen across MitoFinder, NCBI and older literature, after cleaning
# (upper case, "MT" prefix and punctuation removed).
_SYNONYMS = {
    "ATP6": "ATP6", "ATPASE6": "ATP6", "ATPA6": "ATP6",
    "ATP8": "ATP8", "ATPASE8": "ATP8", "ATPA8": "ATP8",
    "COX1": "COX1", "CO1": "COX1", "COI": "COX1", "COXI": "COX1",
    "COX2": "COX2", "CO2": "COX2", "COII": "COX2", "COXII": "COX2",
    "COX3": "COX3", "CO3": "COX3", "COIII": "COX3", "COXIII": "COX3",
    "CYTB": "CYTB", "CYB": "CYTB", "COB": "CYTB",
    "ND4L": "ND4L", "NAD4L": "ND4L", "NADH4L": "ND4L",
}
for _n in "123456":
    for _prefix in ("ND", "NAD", "NADH"):
        _SYNONYMS[_prefix + _n] = "ND" + _n

_ROMAN = {"I": "1", "II": "2", "III": "3"}


def normalize_gene(name: str | None) -> str | None:
    """Map any common spelling of a mitochondrial protein-coding gene to one of
    GENES, or None if it isn't one (rRNA, tRNA, control region, junk)."""
    if not name:
        return None
    clean = re.sub(r"[^A-Z0-9]", "", name.upper())
    if clean.startswith("MT") and clean[2:] in _SYNONYMS:
        clean = clean[2:]
    return _SYNONYMS.get(clean)


def gene_from_product(product: str | None) -> str | None:
    """Gene from an NCBI product description, e.g. 'cytochrome c oxidase subunit II'."""
    if not product:
        return None
    p = product.upper()
    if "CYTOCHROME B" in p:
        return "CYTB"
    m = re.search(r"CYTOCHROME C OXIDASE SUBUNIT\s+(III|II|I|[123])\b", p)
    if m:
        return "COX" + _ROMAN.get(m.group(1), m.group(1))
    m = re.search(r"NADH DEHYDROGENASE SUBUNIT\s+(4L|[1-6])\b", p)
    if m:
        return "ND" + m.group(1)
    m = re.search(r"ATP SYNTHASE (?:F0 |FO )?SUBUNIT\s+([68])\b", p)
    if m:
        return "ATP" + m.group(1)
    return None


def gene_from_header(header: str) -> str | None:
    """Gene from a FASTA header. MitoFinder writes '>Sample@COX1'; other tools
    put the gene in the last word, sometimes after an underscore."""
    header = header.lstrip(">").strip()
    candidates = [header.split("@")[-1]] if "@" in header else []
    words = header.split()
    if words:
        candidates += [words[-1], words[0], words[-1].split("_")[-1]]
    for c in candidates:
        gene = normalize_gene(c)
        if gene:
            return gene
    return gene_from_product(header)


def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    records, header, seq = [], None, []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header, seq = line[1:], []
        else:
            seq.append(line)
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def write_fasta(path: str | Path, records, width: int = 60) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def clean_protein(seq: str) -> str:
    """Drop the trailing stop and turn internal stops/unknowns into X so
    aligners and IQ-TREE accept the sequence."""
    seq = seq.upper().replace("-", "").rstrip("*")
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "X", seq)


def genes_by_taxon(records) -> dict[str, str]:
    """{gene: protein} from one taxon's records, keeping the longest copy when
    a gene appears more than once (split or duplicated annotations)."""
    out: dict[str, str] = {}
    for header, seq in records:
        gene = gene_from_header(header)
        if gene is None:
            continue
        seq = clean_protein(seq)
        if len(seq) > len(out.get(gene, "")):
            out[gene] = seq
    return out


def build_supermatrix(alignments: dict[str, dict[str, str]], taxa: list[str]):
    """Concatenate per-gene alignments ({gene: {taxon: aligned_seq}}) in the
    given gene order. Taxa missing a gene get gaps. Genes with no sequences are
    skipped. Returns ({taxon: row}, [(gene, start, end)]) with 1-based ranges."""
    rows = {t: [] for t in taxa}
    partitions, pos = [], 1
    for gene, aln in alignments.items():
        if not aln:
            continue
        lengths = {len(s) for s in aln.values()}
        if len(lengths) != 1:
            raise ValueError(f"{gene}: sequences have different lengths {sorted(lengths)}; not aligned?")
        length = lengths.pop()
        if length == 0:
            continue
        for t in taxa:
            rows[t].append(aln.get(t, "-" * length))
        partitions.append((gene, pos, pos + length - 1))
        pos += length
    return {t: "".join(parts) for t, parts in rows.items()}, partitions


def partitions_nexus(partitions) -> str:
    lines = ["#nexus", "begin sets;"]
    lines += [f"  charset {gene} = {start}-{end};" for gene, start, end in partitions]
    lines.append("end;")
    return "\n".join(lines) + "\n"


def seq_stats(records) -> dict:
    seqs = [s.upper() for _, s in records]
    total = sum(len(s) for s in seqs)
    gc = sum(s.count("G") + s.count("C") for s in seqs)
    acgt = sum(s.count(b) for s in seqs for b in "ACGT")
    return {
        "contigs": len(seqs),
        "length_bp": total,
        "gc_pct": round(100 * gc / acgt, 1) if acgt else None,
    }
