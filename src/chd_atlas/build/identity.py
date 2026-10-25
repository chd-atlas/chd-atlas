"""Who a gene is, and where else a reader can look it up. Design decision D52.

`mirrors/genes.tsv` has carried every published gene's approved name, locus,
Ensembl and NCBI Gene ids, UniProt accession and MANE Select transcript since
the registry was first mirrored. Until schema 2.17, of those six columns,
`locus`, `ensembl_gene`, `ncbi_gene`, `uniprot` and `mane_select` reached no
published byte, and `name` reached only search's `terms`; a bundle identified
its gene by HGNC id and symbol. Measured 2026-10-25: 0 of the 185 registry
rows lack any of those six columns.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

from chd_atlas.build.emit import Json


@dataclass(frozen=True)
class GeneIdentity:
    """The identity columns of one `mirrors/genes.tsv` row.

    Every field optional: the mirror schema permits a null in every one of these
    columns but `name`, and a gene the registry does not hold has none at all
    (`EMPTY_IDENTITY`).
    `ncbi_gene` is an `int` because the schema reads it as one, and the bundle
    publishes it as one.
    """

    name: str | None = None
    locus: str | None = None
    ensembl_gene: str | None = None
    ncbi_gene: int | None = None
    uniprot: str | None = None
    mane_select: str | None = None


# The identity of a gene the registry does not hold. Its bundle publishes six
# nulls.
EMPTY_IDENTITY: Final = GeneIdentity()


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
