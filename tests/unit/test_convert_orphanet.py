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

from chd_atlas.orphanet_scope import load_orphanet_tables, reach
from chd_atlas.tables import MONDO_XREFS, ORPHANET_OMIM, ORPHANET_PHENOTYPES

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from convert_orphanet import (  # noqa: E402
    OMIM_ORDER,
    PHENOTYPE_ORDER,
    XREF_ORDER,
    _omim_rows,
    _phenotype_rows,
    _xref_rows,
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
