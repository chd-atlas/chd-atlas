"""A gene's identity, and the external records its page links to. D52."""

from __future__ import annotations

from dataclasses import replace
from typing import Final

import pytest

from chd_atlas.build.identity import GeneIdentity, external_links, identity_payload

# NKX2-5's row of `mirrors/genes.tsv`, measured 2026-10-25. Real ids, so each URL
# below names a record checked to exist that day -- by gnomAD's GraphQL API,
# ClinVar's E-utilities, Ensembl's and UniProt's REST services and ClinGen's page
# title, each of which fails on a nonexistent id (`identity.py`'s module
# docstring has the requests and numbers). Not by status code: gnomAD, ClinVar,
# ClinGen and UniProt answer HTTP 200 for any id.
_NKX2_5: Final = "HGNC:2488"
_IDENTITY: Final = GeneIdentity(
    name="NK2 homeobox 5",
    locus="5q35.1",
    ensembl_gene="ENSG00000183072",
    ncbi_gene=1482,
    uniprot="P52952",
    mane_select="NM_004387.4",
)


def test_the_link_row_is_pinned_against_literals_with_the_required_links_first() -> None:
    """The URL shapes are a published contract, so literals, never the templates.

    gnomAD by Ensembl gene id and ClinVar by NCBI Gene id are the two D52
    requires, and they come first. UniProt is the one optional link: a registry
    row without an accession gets no UniProt link rather than a broken one.
    """
    assert [(link.label, link.href) for link in external_links(_NKX2_5, _IDENTITY)] == [
        ("gnomAD", "https://gnomad.broadinstitute.org/gene/ENSG00000183072?dataset=gnomad_r4"),
        ("ClinVar", "https://www.ncbi.nlm.nih.gov/clinvar/?term=1482%5Bgeneid%5D"),
        ("ClinGen", "https://search.clinicalgenome.org/kb/genes/HGNC:2488"),
        ("Ensembl", "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000183072"),
        ("UniProt", "https://www.uniprot.org/uniprotkb/P52952/entry"),
    ]
    without = external_links(_NKX2_5, replace(_IDENTITY, uniprot=None))
    assert [link.label for link in without] == ["gnomAD", "ClinVar", "ClinGen", "Ensembl"]


@pytest.mark.parametrize(
    ("missing", "phrase"),
    [
        pytest.param({"ensembl_gene": None}, "has no Ensembl gene id in", id="ensembl"),
        pytest.param({"ncbi_gene": None}, "has no NCBI Gene id in", id="ncbi"),
        pytest.param(
            {"ensembl_gene": None, "ncbi_gene": None},
            "has no Ensembl gene id or NCBI Gene id in",
            id="both",
        ),
    ],
)
def test_a_gene_missing_a_required_id_is_refused_with_the_id_it_lacks(
    missing: dict[str, None], phrase: str
) -> None:
    """D52's refusal, one case per half of the condition and one for both.

    Each case matches the *specific* phrase, so a message naming the wrong id --
    the curator then looks in the wrong column -- fails as surely as a guard that
    stopped firing. Before this test, deleting the condition or either half of
    it survived all 1,206 tests (code review, 2026-10-25).
    """
    with pytest.raises(ValueError, match=f"^{_NKX2_5} {phrase} mirrors/genes.tsv; D52 requires"):
        external_links(_NKX2_5, replace(_IDENTITY, **missing))


@pytest.mark.parametrize("ensembl", ['ENSG00000183072" onmouseover="x', "ENSG0000018307", "5q35.1"])
def test_a_value_that_is_not_an_ensembl_gene_id_never_reaches_an_href(ensembl: str) -> None:
    """`ensembl_gene` has no pattern in the mirror schema and is interpolated into
    an `href`; escaping makes it inert, not correct (measured 2026-10-25: all 185
    registry rows match `ENSG` + 11 digits, so this fires only on a bad mirror)."""
    with pytest.raises(ValueError, match="not an Ensembl gene id"):
        external_links(_NKX2_5, replace(_IDENTITY, ensembl_gene=ensembl))


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
