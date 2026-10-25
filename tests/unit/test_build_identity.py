"""A gene's identity, and the external records its page links to. D52."""

from __future__ import annotations

from typing import Final

from chd_atlas.build.identity import GeneIdentity, identity_payload

# NKX2-5's row of `mirrors/genes.tsv`, measured 2026-10-25.
_IDENTITY: Final = GeneIdentity(
    name="NK2 homeobox 5",
    locus="5q35.1",
    ensembl_gene="ENSG00000183072",
    ncbi_gene=1482,
    uniprot="P52952",
    mane_select="NM_004387.4",
)


def test_the_identity_payload_carries_every_key_and_an_integer_gene_id() -> None:
    """Every key on every gene, `null` where the registry has nothing.

    `type(...) is int` beside the equality, because `{"ncbi_gene": 1482} ==
    {"ncbi_gene": 1482.0}` is true in Python and the published contract is an
    integer (CLAUDE.md section 5).
    """
    payload = identity_payload(_IDENTITY)
    assert payload == {
        "name": "NK2 homeobox 5",
        "locus": "5q35.1",
        "ensembl_gene": "ENSG00000183072",
        "ncbi_gene": 1482,
        "uniprot": "P52952",
        "mane_select": "NM_004387.4",
    }
    assert type(payload["ncbi_gene"]) is int
    assert identity_payload(GeneIdentity()) == dict.fromkeys(payload)
