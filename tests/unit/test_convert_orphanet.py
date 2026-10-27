# tests/unit/test_convert_orphanet.py
"""The Orphanet/MONDO converter's contracts, driven through its row builders.

The mirror-reading test at the end catches a *stale* mirror; the builder tests
catch a converter defect, which a committed mirror cannot show (CLAUDE.md
section 4.23).
"""

from __future__ import annotations

import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Final

import pytest

from chd_atlas.orphanet_scope import load_orphanet_tables, reach
from chd_atlas.tables import MONDO_XREFS, ORPHANET_OMIM, ORPHANET_PHENOTYPES
from chd_atlas.validate.sources import load_sources

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from convert_orphanet import (  # noqa: E402
    HPO_RELEASE,
    OMIM_ORDER,
    PHENOTYPE_ORDER,
    XREF_ORDER,
    _omim_rows,
    _phenotype_rows,
    _xref_rows,
    jdbor_release,
    record_release,
)

_PRODUCT4_DISORDER = """
<HPODisorderSetStatus id="1">
  <Disorder id="2">
    <OrphaCode>648</OrphaCode>
    <Name lang="en">Noonan syndrome</Name>
    <HPODisorderAssociationList count="3">
      <HPODisorderAssociation id="10">
        <HPO id="4"><HPOId>HP:0001642</HPOId><HPOTerm>Pulmonic stenosis</HPOTerm></HPO>
        <HPOFrequency id="5"><Name lang="en">Frequent (79-30%)</Name></HPOFrequency>
      </HPODisorderAssociation>
      <HPODisorderAssociation id="11">
        <HPO id="6"><HPOId>HP:0000256</HPOId><HPOTerm>Macrocephaly</HPOTerm></HPO>
        <HPOFrequency id="7"><Name lang="en">Very frequent (99-80%)</Name></HPOFrequency>
      </HPODisorderAssociation>
      <HPODisorderAssociation id="12">
        <HPO id="8"><HPOId>HP:0001631</HPOId><HPOTerm>Atrial septal defect</HPOTerm></HPO>
        <HPOFrequency id="9"/>
      </HPODisorderAssociation>
    </HPODisorderAssociationList>
  </Disorder>
</HPODisorderSetStatus>
"""

_PRODUCT1_DISORDER = """
<Disorder id="17601">
  <OrphaCode>648</OrphaCode>
  <Name lang="en">Noonan syndrome</Name>
  <ExternalReferenceList count="2">
    <ExternalReference id="1">
      <Source>ICD-10</Source>
      <Reference>Q87.1</Reference>
      <DisorderMappingRelation id="2">
        <Name lang="en">NTBT (narrower)</Name>
      </DisorderMappingRelation>
      <DisorderMappingValidationStatus id="3">
        <Name lang="en">Validated</Name>
      </DisorderMappingValidationStatus>
    </ExternalReference>
    <ExternalReference id="4">
      <Source>OMIM</Source>
      <Reference>610733</Reference>
      <DisorderMappingRelation id="5">
        <Name lang="en">BTNT (ORPHAcode is broader than the targeted code ...)</Name>
      </DisorderMappingRelation>
      <DisorderMappingValidationStatus id="6">
        <Name lang="en">Validated</Name>
      </DisorderMappingValidationStatus>
    </ExternalReference>
  </ExternalReferenceList>
</Disorder>
"""


def test_a_phenotype_row_keeps_orphanets_band_and_only_the_cardiovascular_annotations() -> None:
    """Pulmonic stenosis kept with its band; macrocephaly dropped; an empty band kept as ''."""
    element = ET.fromstring(_PRODUCT4_DISORDER)
    rows = _phenotype_rows(element, cardiovascular=frozenset({"HP:0001642", "HP:0001631"}))
    assert rows == [
        {
            "orpha": "ORPHA:648",
            "disorder_name": "Noonan syndrome",
            "hpo": "HP:0001642",
            "hpo_label": "Pulmonic stenosis",
            "frequency": "Frequent (79-30%)",
            "source": "orphanet_product4",
        },
        {
            "orpha": "ORPHA:648",
            "disorder_name": "Noonan syndrome",
            "hpo": "HP:0001631",
            "hpo_label": "Atrial septal defect",
            "frequency": "",
            "source": "orphanet_product4",
        },
    ]
    assert tuple(rows[0]) == PHENOTYPE_ORDER == ORPHANET_PHENOTYPES.column_names


def test_an_omim_row_takes_the_relation_code_and_drops_every_other_source() -> None:
    element = ET.fromstring(_PRODUCT1_DISORDER)
    rows = _omim_rows(element)
    assert rows == [
        {
            "orpha": "ORPHA:648",
            "omim": "OMIM:610733",
            "relation": "BTNT",
            "validation": "Validated",
            "source": "orphanet_product1",
        }
    ]
    assert tuple(rows[0]) == OMIM_ORDER == ORPHANET_OMIM.column_names


def test_an_xref_row_is_an_exact_match_to_orphanet_or_omim_for_a_mirrored_term() -> None:
    """MONDO spells Orphanet `Orphanet:`; the mirror spells it `ORPHA:`, as Orphanet does."""
    triples = [
        ("MONDO:0018997", "skos:exactMatch", "Orphanet:648"),
        ("MONDO:0018997", "skos:exactMatch", "OMIM:163950"),
        ("MONDO:0018997", "skos:closeMatch", "Orphanet:98733"),  # not exact
        ("MONDO:0018997", "skos:exactMatch", "DOID:3490"),  # not a target the rule follows
        ("MONDO:0018997", "skos:exactMatch", "OMIMPS:163950"),  # a phenotypic series
        ("MONDO:9999999", "skos:exactMatch", "Orphanet:1"),  # not a mirrored term
    ]
    rows = _xref_rows(triples, terms={"MONDO:0018997"})
    assert rows == [
        {"mondo": "MONDO:0018997", "target": "OMIM:163950", "source": "mondo_sssom"},
        {"mondo": "MONDO:0018997", "target": "ORPHA:648", "source": "mondo_sssom"},
    ]
    assert tuple(rows[0]) == XREF_ORDER == MONDO_XREFS.column_names


def test_the_committed_mirrors_carry_noonan_4s_route_to_its_parent() -> None:
    """MONDO:0012547 -> OMIM:610733 -> ORPHA:648 (BTNT, Validated), measured 2026-10-24."""
    tables = load_orphanet_tables(REPO_ROOT)
    assert tables is not None
    assert reach("MONDO:0012547", "ORPHA:648", tables) == ("subtype", "OMIM:610733")
    assert reach("MONDO:0018997", "ORPHA:648", tables) == ("direct", None)


def test_the_converter_filters_on_the_hpo_release_sources_yaml_pins() -> None:
    """One HPO release for the mirror and for the rule that reads it.

    The converter keeps an Orphanet annotation only if its HPO id falls under
    HP:0001626 in `HPO_RELEASE`; SCP006 and the build decide what is structural
    over the release `mirrors/sources.yaml` pins (`pinned_structural_phenotypes`).
    The two are separate literals -- the converter runs without the package
    installed, so it does not import the registry loader -- and re-pinning HPO
    in one alone would apply the rule over one release to a mirror filtered
    through another: an annotation whose term the new release places under the
    heart would never have been mirrored, and the syndrome it qualifies would
    fail the rule with nothing saying why.
    """
    registry, _ = load_sources(REPO_ROOT)
    assert registry.ontology_files()["HP"] == HPO_RELEASE


_JDBOR = (
    '<?xml version="1.0" encoding="UTF-8"?>\n'
    '<JDBOR date="{date}" version="{version}" copyright="Orphanet (c) 2026">\n'
    "  <Availability/>\n</JDBOR>\n"
)


def test_the_converter_writes_each_files_own_header_beside_its_checksum(tmp_path: Path) -> None:
    """R11 (owner decision 2026-10-26): the data version cannot drift from the file.

    Orphadata's citation asks for a "Data version", and both files state one on
    their root `<JDBOR>` element. `jdbor_release` reads it; `record_release`
    writes it, with the sha256, into the file's own `sources.yaml` entry on
    every converter run, by line edit so the file's comments survive.
    Re-recording the committed values leaves the committed file byte-identical
    -- measured 2026-10-26 with the real headers of both local downloads, which
    is the check this cannot run in CI (`local/` is not committed) -- and a new
    release changes exactly its own three lines. A header with no version is
    refused rather than recorded as nothing.
    """
    product = tmp_path / "en_product4.xml"
    product.write_text(
        _JDBOR.format(date="2027-01-02 03:04:05", version="1.3.50 / 4.1.9 [2026-12-01]")
    )
    assert jdbor_release(product) == ("1.3.50 / 4.1.9 [2026-12-01]", "2027-01-02 03:04:05")

    committed = (REPO_ROOT / "mirrors" / "sources.yaml").read_text(encoding="utf-8")
    sources = tmp_path / "sources.yaml"
    sources.write_text(committed, encoding="utf-8")
    registry, _ = load_sources(REPO_ROOT)
    for source_id in ("orphanet_product4", "orphanet_product1"):
        entry = registry.by_id[source_id]
        assert entry.file_version is not None and entry.file_date is not None
        record_release(sources, source_id, entry.version, entry.file_version, entry.file_date)
    assert sources.read_text(encoding="utf-8") == committed

    record_release(
        sources,
        "orphanet_product4",
        "sha256:" + "a" * 64,
        "1.3.50 / 4.1.9 [2026-12-01]",
        "2027-01-02 03:04:05",
    )
    changed = [
        (old, new)
        for old, new in zip(committed.splitlines(), sources.read_text().splitlines(), strict=True)
        if old != new
    ]
    assert [new for _, new in changed] == [
        '    version: "sha256:' + "a" * 64 + '"',
        '    file_version: "1.3.50 / 4.1.9 [2026-12-01]"',
        '    file_date: "2027-01-02 03:04:05"',
    ]

    product.write_text('<?xml version="1.0"?>\n<JDBOR date="2027-01-02 03:04:05"/>\n')
    with pytest.raises(ValueError, match="no version"):
        jdbor_release(product)

    # Two more refusals, each of which survived every test until review
    # (2026-10-27): a file whose root is not `<JDBOR>` -- some other XML, whose
    # attributes would otherwise be cited as Orphadata's release -- and an
    # entry that does not carry each recorded key exactly once, which a line
    # edit would otherwise rewrite twice or not at all.
    product.write_text(
        _JDBOR.replace("JDBOR", "HPODisorderSetStatusList").format(
            date="2027-01-02 03:04:05", version="1.3.50"
        )
    )
    with pytest.raises(ValueError, match="not <JDBOR>"):
        jdbor_release(product)
    sources.write_text(
        committed.replace("    file_date: ", "    file_date: x\n    file_date: ", 1),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="expected one `file_date:` line, found 2"):
        record_release(sources, "orphanet_product4", "sha256:" + "b" * 64, "v", "d")


_RECORDED_RELEASE: Final = {
    "orphanet_product4": ("1.3.42 / 4.1.8 [2025-03-03]", "2026-06-23 07:57:18"),
    "orphanet_product1": ("1.3.42 / 4.1.8 [2025-03-03]", "2026-06-23 07:53:50"),
}


def test_the_recorded_orphadata_release_is_the_header_of_the_mirrored_files() -> None:
    """The values `sources.json` publishes and the landing page cites, pinned (R11).

    Transcribed from the `<JDBOR>` headers of the two downloads whose sha256s
    the registry records (2026-06-23 release, read 2026-10-26). The downloads
    themselves are read by the test below, where they exist.
    """
    registry, _ = load_sources(REPO_ROOT)
    recorded = {
        source_id: (registry.by_id[source_id].file_version, registry.by_id[source_id].file_date)
        for source_id in ("orphanet_product4", "orphanet_product1")
    }
    assert recorded == _RECORDED_RELEASE


@pytest.mark.parametrize(
    ("source_id", "name"),
    [("orphanet_product4", "en_product4.xml"), ("orphanet_product1", "en_product1.xml")],
)
def test_a_local_orphadata_download_states_the_recorded_release(source_id: str, name: str) -> None:
    """Where the download is present -- a curator's machine, never CI -- it must agree.

    So a re-download recorded by hand rather than by the converter fails here.
    `local/` is gitignored, so in CI the file is absent and the test is a
    visible skip naming it; until review (2026-10-27) it was an `if` that
    passed silently, which read the same as a check that ran.
    """
    path = REPO_ROOT / "local" / "data" / "orphanet" / name
    if not path.is_file():
        pytest.skip(f"local/data/orphanet/{name} is not present (gitignored; a curator's copy)")
    assert jdbor_release(path) == _RECORDED_RELEASE[source_id]
