# tests/unit/test_build_scope.py
"""`scope.json`: the registry a bundle's bare disease ids resolve against."""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import pytest

from chd_atlas.build.emit import Emitter
from chd_atlas.build.scope import build_scope, scope_warrants
from chd_atlas.models.scope import ScopeEntry
from chd_atlas.orphanet_scope import Annotation, OrphanetTables

_PULMONIC = Annotation("HP:0001642", "Pulmonic stenosis", "Frequent (79-30%)")
_TABLES = OrphanetTables(
    annotations={"ORPHA:648": (_PULMONIC,)},
    names={"ORPHA:648": "Noonan syndrome"},
    omim_links={("ORPHA:648", "OMIM:610733"): ("BTNT", "Validated")},
    xrefs={
        "MONDO:0018997": frozenset({"ORPHA:648"}),
        "MONDO:0012547": frozenset({"OMIM:610733"}),
    },
)
_STRUCTURAL = frozenset({"HP:0001642"})


def _entry(id_: str, label: str, **fields: object) -> ScopeEntry:
    payload: dict[str, object] = {
        "id": id_,
        "label": label,
        "reason": "transcribed",
        "admitted_by": "orphanet_phenotype_frequency",
        "attributed_to": "Orphanet",
        "orphanet": "ORPHA:648",
        "admitted_on": date(2026, 10, 25),
    }
    payload.update(fields)
    return ScopeEntry.model_validate(payload)


_HOLT_ORAM = _entry(
    "MONDO:0007732",
    "Holt-Oram syndrome",
    admitted_by="clingen_expert_panel",
    attributed_to="Syndromic Disorders Gene Curation Expert Panel",
    orphanet=None,
)
_NOONAN = _entry("MONDO:0018997", "Noonan syndrome")
_NOONAN_4 = _entry("MONDO:0012547", "Noonan syndrome 4")


def test_scope_json_publishes_every_term_sorted_with_its_warrant_derived_from_the_mirrors(
    tmp_path: Path,
) -> None:
    """Sorted by id, every key on every object, and the warrant read from data.

    The phenotypes are not stored in `chd_scope.yaml`; they are re-derived from
    the Orphanet mirror on every build, so `scope.json` cannot drift from what
    Orphanet records. A ClinGen-path term carries `orphanet: null`, never an
    absent key (CLAUDE.md section 5: check with `is`, not only equality).
    """
    entries = [_NOONAN_4, _HOLT_ORAM, _NOONAN]
    warrants = scope_warrants(entries, _TABLES, _STRUCTURAL)
    emitter = Emitter(root=tmp_path)
    build_scope(entries, warrants, emitter)
    terms = json.loads((tmp_path / "scope.json").read_text())["terms"]

    assert [term["id"] for term in terms] == ["MONDO:0007732", "MONDO:0012547", "MONDO:0018997"]
    assert terms[0]["orphanet"] is None
    assert set(terms[0]) == {"id", "label", "admitted_by", "attributed_to", "orphanet"}
    assert terms[1] == {
        "id": "MONDO:0012547",
        "label": "Noonan syndrome 4",
        "admitted_by": "orphanet_phenotype_frequency",
        "attributed_to": "Orphanet",
        "orphanet": {
            "code": "ORPHA:648",
            "name": "Noonan syndrome",
            "via": "subtype",
            "omim": "OMIM:610733",
            "phenotypes": [
                {
                    "hpo": "HP:0001642",
                    "label": "Pulmonic stenosis",
                    "frequency": "Frequent (79-30%)",
                }
            ],
        },
    }
    assert terms[2]["orphanet"]["via"] == "direct"
    assert terms[2]["orphanet"]["omim"] is None


@pytest.mark.parametrize(
    ("entry", "tables", "structural", "message"),
    [
        pytest.param(
            _entry("MONDO:0018997", "Noonan syndrome", orphanet="ORPHA:3071"),
            _TABLES,
            _STRUCTURAL,
            "MONDO:0018997",
            id="unwarranted-term",
        ),
        pytest.param(_NOONAN, None, _STRUCTURAL, "no Orphanet mirrors or HPO", id="no-mirrors"),
        pytest.param(_NOONAN, _TABLES, None, "no Orphanet mirrors or HPO", id="no-hpo-release"),
    ],
)
def test_a_build_that_reaches_an_unwarranted_orphanet_term_refuses_rather_than_publishing(
    entry: ScopeEntry,
    tables: OrphanetTables | None,
    structural: frozenset[str] | None,
    message: str,
) -> None:
    """A bypassed gate must fail, never publish (CLAUDE.md section 9): SCP006 refuses first.

    And a term with nothing to derive its warrant from is the same bypass, not
    a term with no warrant to state: returning `{}` there would publish every
    gene page the term admits with no Orphanet note and `scope.json` with
    `orphanet: null` for it, on a green build. Measured 2026-10-26 at c034dd6:
    a `return {}` before this raise survived all 634 tests of the release
    review's mutation set, which had no case reaching it.
    """
    with pytest.raises(ValueError, match=message):
        scope_warrants([entry], tables, structural)
