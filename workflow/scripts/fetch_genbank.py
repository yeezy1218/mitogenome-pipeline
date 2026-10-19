"""Download GenBank records from NCBI Entrez."""

import sys
import time

from Bio import Entrez


def main(snakemake):
    sys.stderr = open(snakemake.log[0], "w")
    if not snakemake.params.email:
        raise SystemExit("Set entrez_email in config/config.yaml (NCBI asks for a contact address).")
    Entrez.email = snakemake.params.email
    accessions = list(snakemake.params.accessions)
    handle = Entrez.efetch(db="nucleotide", id=",".join(accessions), rettype="gbwithparts", retmode="text")
    text = handle.read()
    handle.close()
    found = text.count("\nLOCUS") + text.startswith("LOCUS")
    if found != len(accessions):
        raise SystemExit(f"Asked NCBI for {len(accessions)} records, got {found}: {accessions}")
    with open(snakemake.output[0], "w") as fh:
        fh.write(text)
    print(f"Fetched {found} record(s): {', '.join(accessions)}", file=sys.stderr)
    time.sleep(0.4)  # stay under NCBI's 3 requests/second without an API key


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
