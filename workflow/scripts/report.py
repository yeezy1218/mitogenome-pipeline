"""Self-contained HTML report: tree figure, assembly table, gene presence matrix."""

import base64
import html
import io
import re
import sys
from datetime import date

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from Bio import Phylo  # noqa: E402

CSS = """
body { font-family: system-ui, sans-serif; max-width: 1100px; margin: 2rem auto; padding: 0 1rem; color: #1d1f23; }
h1 { margin-bottom: .2rem; } .sub { color: #5b6270; margin-top: 0; }
table { border-collapse: collapse; font-size: .9rem; margin: 1rem 0; }
th, td { border: 1px solid #d5d9e0; padding: .3rem .55rem; text-align: right; }
th { background: #f1f3f6; } td:first-child, th:first-child { text-align: left; }
td.yes { background: #2f9e66; color: white; text-align: center; }
td.no { background: #f6e1e1; color: #a33; text-align: center; }
img { max-width: 100%; } pre { background: #f6f7f9; padding: .8rem; overflow-x: auto; font-size: .8rem; }
"""


def tree_png(tree_path, outgroup):
    tree = Phylo.read(tree_path, "newick")
    if outgroup:
        tree.root_with_outgroup(outgroup)
    tree.ladderize()
    n = tree.count_terminals()
    fig = plt.figure(figsize=(9, max(3, 0.35 * n + 1)))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        label_func=lambda c: c.name.replace("_", " ") if c.is_terminal() else "",
        branch_labels=lambda c: (
            f"{c.confidence:.0f}" if c.confidence is not None and not c.is_terminal() else ""
        ),
    )
    ax.set_ylabel("")
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel("substitutions per site")
    buf = io.BytesIO()
    fig.tight_layout()
    fig.savefig(buf, format="png", dpi=150)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def best_models(iqtree_text):
    """{gene: model} from IQ-TREE's .iqtree summary. A partitioned run lists
    them as 'model:gene,model:gene,...'."""
    m = re.search(r"Best-fit model according to BIC:\s*(\S+)", iqtree_text)
    if not m:
        return {}
    pairs = [item.rsplit(":", 1) for item in m.group(1).split(",")]
    return {p[1]: p[0] for p in pairs if len(p) == 2} or {"all": m.group(1)}


def models_table(models):
    if not models:
        return "<p>See results/phylo/mito.iqtree.</p>"
    head = "".join(f"<th>{html.escape(g)}</th>" for g in models)
    cells = "".join(f"<td>{html.escape(m)}</td>" for m in models.values())
    return f"<table><tr><th>Gene</th>{head}</tr><tr><td>Model</td>{cells}</tr></table>"


def presence_table(stats, genes):
    head = "".join(f"<th>{g}</th>" for g in genes)
    rows = []
    for _, r in stats.iterrows():
        cells = "".join(
            f'<td class="yes">&#10003;</td>' if bool(r[g]) else '<td class="no">&ndash;</td>' for g in genes
        )
        rows.append(f"<tr><td>{html.escape(r['species'])}</td>{cells}</tr>")
    return f"<table><tr><th>Taxon</th>{head}</tr>{''.join(rows)}</table>"


def main(snakemake):
    sys.stderr = open(snakemake.log[0], "w")
    p = snakemake.params
    stats = pd.read_csv(snakemake.input.stats, sep="\t", keep_default_na=False)
    for g in p.genes:
        stats[g] = stats[g].astype(str) == "True"
    summary = open(snakemake.input.summary).read()
    png = tree_png(snakemake.input.tree, p.outgroup)
    n_sites = re.findall(r"= \d+-(\d+);", open(snakemake.input.partitions).read())

    table = stats[["species", "source", "contigs", "length_bp", "gc_pct", "genes_found"]].rename(
        columns={"species": "Taxon", "source": "Source", "contigs": "Contigs", "length_bp": "Length (bp)",
                 "gc_pct": "GC %", "genes_found": f"Genes (of {len(p.genes)})"}
    )
    qc = '<p><a href="qc/multiqc_report.html">Read quality report (MultiQC)</a></p>' if p.has_qc else ""
    doc = f"""<!doctype html><html><head><meta charset="utf-8">
<title>Mitogenome phylogeny</title><style>{CSS}</style></head><body>
<h1>Mitogenome phylogeny</h1>
<p class="sub">{len(stats)} taxa &middot; {n_sites[-1] if n_sites else '?'} amino-acid sites in
{len(n_sites)} genes &middot; generated {date.today().isoformat()}</p>
<h2>Maximum-likelihood tree</h2>
<p>Partitioned IQ-TREE analysis of concatenated mitochondrial proteins; node labels are ultrafast bootstrap support.</p>
<img alt="Phylogenetic tree" src="data:image/png;base64,{png}">
<h3>Best-fit substitution model per gene (BIC)</h3>{models_table(best_models(summary))}
<h2>Assemblies</h2>{table.to_html(index=False, border=0, na_rep="")}
{qc}
<h2>Protein-coding genes recovered</h2>{presence_table(stats, p.genes)}
</body></html>"""
    open(snakemake.output[0], "w", encoding="utf-8").write(doc)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
