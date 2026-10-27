"""Who a gene is, and where else a reader can look it up. Design decision D52.

`mirrors/genes.tsv` has carried every published gene's approved name, locus,
Ensembl and NCBI Gene ids, UniProt accession and MANE Select transcript since
the registry was first mirrored. Until schema 2.17, of those six columns,
`locus`, `ensembl_gene`, `ncbi_gene`, `uniprot` and `mane_select` reached no
published byte, and `name` reached only search's `terms`; a bundle identified
its gene by HGNC id and symbol. Measured 2026-10-25: 0 of the 185 registry
rows lack any of those six columns.

**Two links are required, and a page that cannot carry them is refused rather
than published without them.** A reader leaving a gene page for variant-level
data goes to gnomAD for population frequencies and to ClinVar for clinical
interpretations. gnomAD is keyed on the Ensembl gene id and ClinVar on the NCBI
Gene id -- not on the symbol. HGNC retires symbols, and ClinVar runs a bare
symbol as a search of every field: measured 2026-10-25, `NKX2-5` matched 901
records where `1482[geneid]` matched 891, and the 10 extra are all variants in
other genes (four in TBX5, the rest in TBX20, HAND2, HCN4 and others) whose
records merely mention NKX2-5. `external_links` raises on a gene missing
either id, so a regenerated mirror that drops one fails the build instead of
shipping pages with the link silently gone.

The URL shapes are a published contract, pinned against literals in
`tests/unit/test_build_identity.py`. **A status code checks none of them**:
gnomAD, ClinVar, ClinGen and UniProt serve their pages from applications that
answer HTTP 200 for an id that does not exist (gnomAD does for
ENSG00000000000). Each was checked on 2026-10-25, for NKX2-5, by a request
that does fail on a nonexistent id:

- gnomAD: the GraphQL API's `gene(gene_id: "ENSG00000183072",
  reference_genome: GRCh38)` returns symbol NKX2-5, and "Gene not found" for
  ENSG00000000000. Its `DatasetId` enum lists `gnomad_r4` and no later release.
- ClinVar: the link's own query through E-utilities -- `esearch` with
  `db=clinvar&term=1482[geneid]` -- counts 891 records, and 0 for
  `99999999[geneid]`.
- Ensembl: REST `lookup/id/ENSG00000183072` returns display_name NKX2-5, and
  HTTP 400 "not found" for ENSG00000000000. The Summary URL redirects (308) to
  Ensembl's feature-explorer page, which answers 200 for NKX2-5 and 404 for
  ENSG00000000000.
- UniProt: REST `uniprotkb/P52952.json` returns primaryAccession P52952 with
  gene NKX2-5, and HTTP 404 for Q00000; the entry page answers 200 for both.
- ClinGen: no API was checked. The gene page answers 200 for HGNC:99999999
  too, but is rendered on the server, and its `<title>` reads "NKX2-5 curation
  results" for HGNC:2488 and "ClinGen Web" for the nonexistent id.

OMIM is absent: it needs a MIM number `genes.tsv` does not carry.

`GNOMAD_DATASET` pins the release rather than leaving it to gnomAD's default:
an unpinned link changes what it shows with no change to this atlas, and
gnomAD keeps earlier releases browsable. Re-check and bump it when gnomAD
publishes a new major release -- the GraphQL `DatasetId` enum is where one
appears.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from chd_atlas.build.emit import Json
from chd_atlas.identifiers import ENSEMBL_GENE_PATTERN
from chd_atlas.tables import TABLE_SCHEMAS, cell_text, read_table

GNOMAD_DATASET: Final = "gnomad_r4"

# Ensembl's stable human gene id. The mirror schema constrains nothing about
# this column and its value is interpolated into an `href`, so it is checked
# here, where it becomes a URL. Measured 2026-10-25: all 185 registry rows match.
_ENSEMBL_GENE: Final = re.compile(ENSEMBL_GENE_PATTERN)


@dataclass(frozen=True)
class GeneIdentity:
    """The identity columns of one `mirrors/genes.tsv` row.

    Every field optional: the mirror schema permits a null in every one of these
    columns but `name`, and a gene the registry does not hold has none at all
    (`EMPTY_IDENTITY`). The two that matter are required by `external_links`,
    where a missing one fails the build rather than costing a link.
    `ncbi_gene` is an `int` because the schema reads it as one, and the bundle
    publishes it as one.
    """

    name: str | None = None
    locus: str | None = None
    ensembl_gene: str | None = None
    ncbi_gene: int | None = None
    uniprot: str | None = None
    mane_select: str | None = None


# The identity of a gene the registry does not hold. Its page is refused by
# `external_links`; its bundle publishes six nulls.
EMPTY_IDENTITY: Final = GeneIdentity()


@dataclass(frozen=True)
class ExternalLink:
    """One outbound link: what a reader sees, and where it goes."""

    label: str
    href: str


def read_gene_identities(root: Path) -> dict[str, GeneIdentity] | None:
    """The identity columns of `mirrors/genes.tsv`, keyed on HGNC id. D52, R13.

    **The one reader, for the build and for `chd-atlas validate` alike.**
    GEN004 (`validate/genes.py::validate_identity_ids`) reports, by
    `identity_id_problems`, what `external_links` refuses a page on; until
    review (2026-10-27) each side read the registry through a copy of its own,
    and a mutant to the validator's copy survived every test, so a padded or
    blank cell could be read one way to validate and another to build. Every
    cell is cleaned by `tables.cell_text`, the rule the build's labels use: a
    blank or whitespace-only cell is absent, and a row with no `hgnc_id` has
    nothing to key on and is dropped. `ncbi_gene` arrives as an `int`, as the
    bundle publishes it.

    `None` where the file is absent or unreadable, or lacks `hgnc_id`,
    `ensembl_gene` or `ncbi_gene` -- each reported by a check of its own
    (TBL008 for the first three, TBL001 for a missing column) -- so GEN004
    does not report every published gene for the file's one fault. A registry
    that reads and holds no row is `{}`, not `None`: GEN004 reports that once.
    """
    path = root / "mirrors" / "genes.tsv"
    if not path.is_file():
        return None
    frame, _ = read_table(path, TABLE_SCHEMAS["genes"])
    if frame is None or any(
        column not in frame.columns for column in ("hgnc_id", "ensembl_gene", "ncbi_gene")
    ):
        return None

    identities: dict[str, GeneIdentity] = {}
    for row in frame.to_dicts():
        gene = cell_text(row.get("hgnc_id"))
        if gene is None:
            continue
        ncbi = row.get("ncbi_gene")
        identities[gene] = GeneIdentity(
            name=cell_text(row.get("name")),
            locus=cell_text(row.get("locus")),
            ensembl_gene=cell_text(row.get("ensembl_gene")),
            ncbi_gene=int(ncbi) if ncbi is not None else None,
            uniprot=cell_text(row.get("uniprot")),
            mane_select=cell_text(row.get("mane_select")),
        )
    return identities


def identity_payload(identity: GeneIdentity) -> dict[str, Json]:
    """The bundle's `identity` object: every key on every gene, `null` where absent.

    Ids only, no URLs. The ids are the facts; a URL template is this site's
    presentation choice, and publishing it would make a third party's URL
    scheme part of this atlas's schema.
    """
    return {
        "name": identity.name,
        "locus": identity.locus,
        "ensembl_gene": identity.ensembl_gene,
        "ncbi_gene": identity.ncbi_gene,
        "uniprot": identity.uniprot,
        "mane_select": identity.mane_select,
    }


def identity_id_problems(identity: GeneIdentity) -> tuple[list[str], list[str]]:
    """The two ids a gene page links by, checked: `(missing, malformed)`. D52, R13.

    **The one rule, read in two places** (owner decision 2026-10-26, R13):
    `external_links` refuses a published gene's page on it, and GEN004
    (`validate/genes.py`) reports it for every published gene, so `chd-atlas
    validate` and the build cannot disagree about an identity row. Until then
    a malformed Ensembl id passed validate with 0 errors and the build refused
    the repository after writing had begun -- measured 2026-10-26 on a copy
    with TBX5's id cut to `ENSG0089225`, and likewise for a missing NCBI Gene
    id. `missing` names the absent ids; `malformed` says what is wrong with a
    present one: an Ensembl id that is not `ENSG` + 11 digits, or an NCBI Gene
    id below 1, which a ClinVar link would query as `0[geneid]` (the build
    accepted that until this rule was shared).
    """
    missing = [
        label
        for label, value in (
            ("Ensembl gene id", identity.ensembl_gene),
            ("NCBI Gene id", identity.ncbi_gene),
        )
        if value is None
    ]
    malformed: list[str] = []
    if identity.ensembl_gene is not None and not _ENSEMBL_GENE.fullmatch(identity.ensembl_gene):
        malformed.append(f"{identity.ensembl_gene!r} is not an Ensembl gene id (ENSG + 11 digits)")
    if identity.ncbi_gene is not None and identity.ncbi_gene < 1:
        malformed.append(f"{identity.ncbi_gene} is not an NCBI Gene id (a positive integer)")
    return missing, malformed


def external_links(gene: str, identity: GeneIdentity) -> tuple[ExternalLink, ...]:
    """The gene page's link row: gnomAD and ClinVar first, both required. D52.

    Raises `ValueError` naming the gene and the missing id rather than leaving
    a link out -- the failure this project calls its characteristic one is a
    green build that quietly publishes less than it promised. `raise`, never
    `assert`: `-O` strips `assert`. The rule is `identity_id_problems`, which
    GEN004 reports from `chd-atlas validate` before a build reaches here.
    """
    missing, malformed = identity_id_problems(identity)
    ensembl, ncbi = identity.ensembl_gene, identity.ncbi_gene
    if missing or ensembl is None or ncbi is None:
        raise ValueError(
            f"{gene} has no {' or '.join(missing)} in mirrors/genes.tsv; D52 requires every "
            f"published gene page to link gnomAD and ClinVar, which are keyed on them"
        )
    if malformed:
        raise ValueError(f"{gene}: {'; '.join(malformed)}")
    links = [
        ExternalLink(
            "gnomAD", f"https://gnomad.broadinstitute.org/gene/{ensembl}?dataset={GNOMAD_DATASET}"
        ),
        ExternalLink("ClinVar", f"https://www.ncbi.nlm.nih.gov/clinvar/?term={ncbi}%5Bgeneid%5D"),
        ExternalLink("ClinGen", f"https://search.clinicalgenome.org/kb/genes/{gene}"),
        ExternalLink("Ensembl", f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={ensembl}"),
    ]
    if identity.uniprot is not None:
        links.append(
            ExternalLink("UniProt", f"https://www.uniprot.org/uniprotkb/{identity.uniprot}/entry")
        )
    return tuple(links)
