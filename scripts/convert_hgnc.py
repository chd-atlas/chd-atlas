# scripts/convert_hgnc.py
"""Convert the HGNC complete set into `mirrors/genes.tsv`.

Run locally by a curator, never by CI or the build -- design decision D6: CI
validates, builds and deploys; it never authors data. The output is committed,
so a refresh is a reviewable diff rather than a silent change under a build.

    curl -sSL -o /tmp/hgnc.tsv \\
        https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
    .venv/bin/python scripts/convert_hgnc.py /tmp/hgnc.tsv mirrors/genes.tsv

With no third argument the keep set is every gene named by a row of
`mirrors/clingen_gene_validity.tsv` or `mirrors/gencc_submissions.tsv` whose
disease is listed in `curation/chd_scope.yaml` -- 185 genes as committed (154
until the Orphanet scope path widened the keep set on 2026-10-25). The
registry is deliberately wider than what the site publishes (119 genes under
D21), because it is a lookup table rather than published content: widening the
gate later must not require re-mirroring.

`prev_symbols` is a straight copy of HGNC's `prev_symbol`, and it is the
column #33 asked for and this script did not have. HGNC retires a name by
moving it from `symbol` to `prev_symbol`, and essentially never to
`alias_symbol` -- measured over the whole 2026-07-31 complete set, 3 of 15,882
previous symbols are also listed as an alias of the same gene. So a mirror
carrying aliases alone cannot resolve a source still using an old name.
Measured the same day: 56 of the 154
in-scope genes hold a previous symbol -- 70 entries, 69 distinct, the
repeat being `MADH7` which SMAD6 and SMAD7 both retired -- and exactly one of
those (`BVES`, HGNC:1152) is already named by a committed mirror that the
registry could not resolve -- the standing `GEN003`. Copied verbatim, pipes
and all, exactly as `aliases` is: `genes.py` owns the precedence between the
two, and a converter that pre-merged them would destroy the distinction that
precedence turns on.

Two columns need a rule, and both were measured against the live file rather
than assumed:

`mane_select` is published by HGNC as "ENST00000405440.7|NM_181486.4" for all
185 in-scope genes (re-measured 2026-10-25). The mirror carries the RefSeq half, which is what the
committed TBX5 row already holds and what a variant's `mane_transcript` is
written against.

`uniprot_ids` carries more than one accession for exactly one of the 185 --
RBM10, "P0DW28|P98175" (re-measured 2026-10-25 over the 2026-10-16 complete
set) -- and RBM10 is a gene D21 publishes. `UNIPROT_PATTERN` rejects a pipe, so
the mirror needs one, and HGNC's order is alphabetical rather than a ranking.
Until 2026-10-25 the first was taken, which mirrored P0DW28: MINAS-60, a
130-aa microprotein from an alternative reading frame inside RBM10, not
RBM10's own 930-aa P98175. That was tolerable while the column fed only
accession-keyed joins over mirrors that do not exist. Since schema 2.17 it is
the gene's published `identity.uniprot`, and what `build/identity.py` builds
the gene page's UniProt link from, so a wrong accession is a false published
claim -- and a well-formed one fails no check.

So the choice is a human's, never a position's. `_UNIPROT_OVERRIDES` pins
RBM10 to P98175 with the UniProt evidence beside it; any other multi-valued
cell refuses the conversion, naming the gene and its accessions; and a pin
whose accession HGNC no longer lists for that gene refuses too. Null is not a
fallback either: it would read as a gene with no protein.
`tests/unit/test_convert_hgnc.py` pins all three.
"""

from __future__ import annotations

import re
import sys
from collections.abc import Collection
from pathlib import Path

import polars as pl

# The source-to-mirror column mapping, for the reader. `convert` spells the
# same mapping out as expressions; a column HGNC renames raises out of polars
# there rather than mirroring a silently absent field.
_COLUMNS: dict[str, str] = {
    "hgnc_id": "hgnc_id",
    "symbol": "symbol",
    "name": "name",
    "alias_symbol": "aliases",
    "prev_symbol": "prev_symbols",
    "ensembl_gene_id": "ensembl_gene",
    "entrez_id": "ncbi_gene",
    "location": "locus",
    "uniprot_ids": "uniprot",
    "mane_select": "mane_select",
}

# The order `TABLE_SCHEMAS["genes"]` declares. Written out rather than imported
# so this script stays runnable without the package installed.
_ORDER = (
    "hgnc_id",
    "symbol",
    "name",
    "aliases",
    "prev_symbols",
    "ensembl_gene",
    "ncbi_gene",
    "locus",
    "uniprot",
    "mane_select",
)

_SCOPE_ID = re.compile(r"^\s*- id: (MONDO:\d+)", re.M)

# The gene's own protein, for every kept gene whose `uniprot_ids` lists more
# than one accession. HGNC orders the list alphabetically, not by which protein
# is the gene's, so position decides nothing; each entry is a human reading
# UniProt, and carries what they read.
_UNIPROT_OVERRIDES: dict[str, str] = {
    # RBM10. UniProt REST, 2026-10-25: P98175 is RBM10_HUMAN, "RNA-binding
    # protein 10", reviewed, 930 aa, gene RBM10. P0DW28 is MNS60_HUMAN,
    # "Ribosome biogenesis inhibitor MINAS-60", reviewed, 130 aa, also gene
    # RBM10 -- "the product of an alternative open reading frame (alt-ORF) of
    # transcripts coding for the RBM10 (AC P98175) protein". HGNC lists both,
    # "P0DW28|P98175", so first-of-list mirrored the microprotein.
    "HGNC:9896": "P98175",
}


def _uniprot(wanted: pl.DataFrame) -> pl.Expr:
    """The one UniProt accession each kept gene mirrors, refusing any it cannot know.

    A multi-valued cell with no `_UNIPROT_OVERRIDES` entry raises: a new one
    must reach a human, because the only rule a script has is position, and
    position is how RBM10 published the wrong protein. An override whose
    accession HGNC no longer lists for that gene raises too, so a stale pin
    cannot publish a mapping the authority has withdrawn. Checked before
    anything is written, so a refused refresh leaves the committed mirror alone.
    """
    for gene, symbol, cell in wanted.select("hgnc_id", "symbol", "uniprot_ids").rows():
        listed = cell.split("|") if cell is not None else []
        pinned = _UNIPROT_OVERRIDES.get(gene)
        if pinned is not None and pinned not in listed:
            raise ValueError(
                f"{gene} ({symbol}): the pinned UniProt accession {pinned} is not among "
                f"those HGNC lists ({cell}); re-read UniProt and update _UNIPROT_OVERRIDES"
            )
        if pinned is None and len(listed) > 1:
            raise ValueError(
                f"{gene} ({symbol}) lists {len(listed)} UniProt accessions ({cell}); "
                f"pin the gene's own protein in _UNIPROT_OVERRIDES, with the evidence"
            )
    return pl.col("hgnc_id").replace_strict(_UNIPROT_OVERRIDES, default=pl.col("uniprot_ids"))


def _refseq(column: str) -> pl.Expr:
    """The RefSeq half of HGNC's "ENST…|NM_…" MANE Select pair.

    Selected by prefix rather than by position: HGNC publishes the Ensembl
    transcript first today, and a value that is only one of the two must not
    silently publish an Ensembl id in a column every consumer reads as RefSeq.
    A pair with no RefSeq half yields null, which `TABLE_SCHEMAS["genes"]`
    permits.
    """
    return (
        pl.col(column)
        .str.split("|")
        .list.eval(pl.element().filter(pl.element().str.starts_with("NM_")))
        .list.first()
    )


def convert(source: Path, out: Path, keep: Collection[str]) -> int:
    """Write the mirror rows for `keep`. Returns the number of rows written."""
    frame = pl.read_csv(source, separator="\t", infer_schema_length=0)
    wanted = frame.filter(pl.col("hgnc_id").is_in(list(keep)))

    rows = (
        wanted.select(
            pl.col("hgnc_id"),
            pl.col("symbol"),
            pl.col("name"),
            pl.col("alias_symbol").alias("aliases"),
            pl.col("prev_symbol").alias("prev_symbols"),
            pl.col("ensembl_gene_id").alias("ensembl_gene"),
            pl.col("entrez_id").alias("ncbi_gene"),
            pl.col("location").alias("locus"),
            _uniprot(wanted).alias("uniprot"),
            _refseq("mane_select").alias("mane_select"),
        )
        .select(_ORDER)
        # Sorted on the published string, matching `sort_order.py`, so the
        # committed file is already in the order validation requires.
        .sort("hgnc_id")
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    rows.write_csv(out, separator="\t", include_header=True)
    return len(rows)


def _in_scope_genes(root: Path) -> set[str]:
    """Every gene either validity mirror names against an in-scope disease."""
    scope = set(_SCOPE_ID.findall((root / "curation" / "chd_scope.yaml").read_text()))
    genes: set[str] = set()
    for name in ("clingen_gene_validity.tsv", "gencc_submissions.tsv"):
        frame = pl.read_csv(root / "mirrors" / name, separator="\t", infer_schema_length=0)
        genes.update(frame.filter(pl.col("disease").is_in(list(scope)))["gene"].to_list())
    return genes


def main(argv: list[str]) -> int:
    if not 2 <= len(argv) <= 3:
        print(__doc__, file=sys.stderr)
        return 2
    source, out = Path(argv[0]), Path(argv[1])
    root = Path(argv[2]) if len(argv) == 3 else Path.cwd()
    written = convert(source, out, keep=_in_scope_genes(root))
    print(f"wrote {written} gene(s) to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
