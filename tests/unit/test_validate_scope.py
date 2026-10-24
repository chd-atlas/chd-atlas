# tests/unit/test_validate_scope.py
"""MONDO is not pinned as an ontology, so scope terms are checked against the
mirror instead of against a resolver: see `validate/scope.py` for why.
"""

from __future__ import annotations

from datetime import date

import pytest

from chd_atlas.issues import Severity
from chd_atlas.models.scope import ScopeEntry
from chd_atlas.orphanet_scope import Annotation, OrphanetTables
from chd_atlas.validate.scope import (
    CHD_EXPERT_PANEL,
    orphanet_candidates,
    scope_candidates,
    validate_curation_is_in_scope,
    validate_orphanet_warrants,
    validate_scope_attribution,
    validate_scope_terms,
)

_LOCATION = "curation/chd_scope.yaml"


def _entry(**overrides: object) -> ScopeEntry:
    payload: dict[str, object] = {
        "id": "MONDO:0007732",
        "label": "Holt-Oram syndrome",
        "reason": "cardiac septation defects are a defining feature",
        "admitted_by": "clingen_chd_panel",
        "attributed_to": "Congenital Heart Disease Gene Curation Expert Panel",
        "admitted_on": date(2026, 8, 3),
    }
    payload.update(overrides)
    return ScopeEntry.model_validate(payload)


def test_a_term_absent_from_the_mirror_is_reported() -> None:
    issues = validate_scope_terms([_entry()], mirrored={}, location=_LOCATION)

    assert [i.code for i in issues] == ["SCP001"]
    assert issues[0].severity is Severity.ERROR
    assert "MONDO:0007732" in issues[0].message


def test_a_disagreeing_label_is_reported_naming_both_labels() -> None:
    mirrored = {"MONDO:0007732": "Holt-Oram syndrome, ClinGen wording"}

    issues = validate_scope_terms([_entry()], mirrored=mirrored, location=_LOCATION)

    assert [i.code for i in issues] == ["SCP002"]
    assert issues[0].severity is Severity.ERROR
    assert "Holt-Oram syndrome" in issues[0].message
    assert "Holt-Oram syndrome, ClinGen wording" in issues[0].message


def test_no_mirror_reports_exactly_one_skip_warning_and_checks_nothing_else() -> None:
    """`mirrored=None` means the validity mirror could not be read at all.

    Checking terms against `{}` would report SCP001 for every one of them,
    burying the single cause -- the same cascade REF000/SRC000/ONT000 exist to
    avoid. So a bad label on an unreadable-mirror run must not also surface as
    SCP002; only the one skip warning is reported.
    """
    entries = [_entry(), _entry(id="MONDO:0005453", label="congenital heart disease")]

    issues = validate_scope_terms(entries, mirrored=None, location=_LOCATION)

    assert len(issues) == 1
    assert issues[0].code == "SCP000"
    assert issues[0].severity is Severity.WARNING


def test_a_matching_term_is_accepted() -> None:
    mirrored = {"MONDO:0007732": "Holt-Oram syndrome"}

    assert validate_scope_terms([_entry()], mirrored=mirrored, location=_LOCATION) == []


# --- scope_candidates -----------------------------------------------------


def test_net1_finds_a_cardiac_label_and_excludes_admitted_and_non_cardiac_terms() -> None:
    """Net 1 (label keyword) must neither re-report an already-admitted term
    nor flag a term whose label carries no cardiac word.
    """
    mirrored = {
        "MONDO:0001": "atrial septal defect",  # candidate: cardiac, not in scope
        "MONDO:0005453": "congenital heart disease",  # already admitted: excluded
        "MONDO:0002": "diabetes mellitus",  # non-cardiac: excluded
    }
    scope_terms = {"MONDO:0005453"}

    issues = scope_candidates(
        scope_terms=scope_terms,
        mirrored=mirrored,
        diseases_by_gene={},
        in_scope_genes=set(),
        location=_LOCATION,
    )

    assert [i.code for i in issues] == ["SCP003"]
    assert issues[0].severity is Severity.WARNING
    assert "MONDO:0001" in issues[0].message
    assert "MONDO:0005453" not in issues[0].message
    assert "MONDO:0002" not in issues[0].message


def test_net2_finds_a_gene_bridge_candidate_that_net1_misses() -> None:
    """Holt-Oram syndrome's label carries no cardiac keyword, so net 1 is
    blind to it -- this is the measured defect the two-net design fixes. Net
    2 still finds it because TBX5 is already in scope for a different term.
    """
    mirrored = {
        "MONDO:0005453": "congenital heart disease",  # TBX5's in-scope term
        "MONDO:0007732": "Holt-Oram syndrome",  # TBX5's other curation
    }
    scope_terms = {"MONDO:0005453"}
    diseases_by_gene = {"HGNC:11604": {"MONDO:0005453", "MONDO:0007732"}}
    in_scope_genes = {"HGNC:11604"}

    issues = scope_candidates(
        scope_terms=scope_terms,
        mirrored=mirrored,
        diseases_by_gene=diseases_by_gene,
        in_scope_genes=in_scope_genes,
        location=_LOCATION,
    )

    assert [i.code for i in issues] == ["SCP003"]
    assert "net 2" in issues[0].message
    assert "MONDO:0007732" in issues[0].message


def test_at_most_two_issues_with_candidates_sorted() -> None:
    mirrored = {
        "MONDO:0009": "vascular ring",
        "MONDO:0001": "atrial septal defect",
        "MONDO:0005": "aortic coarctation",
    }
    scope_terms: set[str] = set()
    diseases_by_gene = {
        "HGNC:2": {"MONDO:0099", "MONDO:0011"},
        "HGNC:1": {"MONDO:0011"},
    }
    in_scope_genes = {"HGNC:1", "HGNC:2"}

    issues = scope_candidates(
        scope_terms=scope_terms,
        mirrored=mirrored,
        diseases_by_gene=diseases_by_gene,
        in_scope_genes=in_scope_genes,
        location=_LOCATION,
    )

    assert len(issues) == 2
    net1 = next(i for i in issues if "net 1" in i.message)
    net2 = next(i for i in issues if "net 2" in i.message)
    # MONDO:0001 sorts before MONDO:0005, which sorts before MONDO:0009.
    assert net1.message.index("MONDO:0001") < net1.message.index("MONDO:0005")
    assert net1.message.index("MONDO:0005") < net1.message.index("MONDO:0009")
    assert "MONDO:0011" in net2.message
    assert "MONDO:0099" in net2.message
    assert net2.message.index("MONDO:0011") < net2.message.index("MONDO:0099")


def test_nothing_to_report_is_an_empty_list() -> None:
    mirrored = {"MONDO:0005453": "congenital heart disease"}
    scope_terms = {"MONDO:0005453"}
    diseases_by_gene = {"HGNC:11604": {"MONDO:0005453"}}
    in_scope_genes = {"HGNC:11604"}

    issues = scope_candidates(
        scope_terms=scope_terms,
        mirrored=mirrored,
        diseases_by_gene=diseases_by_gene,
        in_scope_genes=in_scope_genes,
        location=_LOCATION,
    )

    assert issues == []


# --- validate_curation_is_in_scope ----------------------------------------


def test_an_out_of_scope_curated_gene_is_reported_at_its_file() -> None:
    curated_genes = {"HGNC:11604": "curation/assertions/TBX5.yaml"}

    issues = validate_curation_is_in_scope(curated_genes, in_scope_genes=set())

    assert [i.code for i in issues] == ["SCP004"]
    assert issues[0].severity is Severity.ERROR
    assert issues[0].location == "curation/assertions/TBX5.yaml"
    assert "HGNC:11604" in issues[0].message


def test_an_in_scope_curated_gene_is_accepted() -> None:
    curated_genes = {"HGNC:11604": "curation/assertions/TBX5.yaml"}

    issues = validate_curation_is_in_scope(curated_genes, in_scope_genes={"HGNC:11604"})

    assert issues == []


# --- validate_scope_attribution (SCP005) -----------------------------------


def test_an_attribution_the_mirror_does_not_support_is_reported() -> None:
    """The check that turns scope provenance from a claim into a fact.

    Until 2026-08-06 every entry in `curation/chd_scope.yaml` named the project
    owner as the admitting authority, which made a non-clinician the authority
    on what counts as congenital heart disease. `ScopeAuthority` now makes a
    person's name unrepresentable -- but an enum constrains only the *shape* of
    the claim. A curator could still write `clingen_chd_panel` beside a term
    that panel has never curated, and without this check nothing would notice,
    while every published page asserts the term came from an external authority.

    ERROR, not WARNING: a false attribution is a provenance claim the atlas
    cannot support, and `_SCOPE_RULE` publishes it on all 24 pages.

    The truthful entry is asserted in the same test rather than a separate one:
    a check that reports everything passes this file's negative case too, and
    only the pair distinguishes a working check from a broken one.
    """
    panels = {"MONDO:0005453": frozenset({"Congenital Heart Disease Gene Curation Expert Panel"})}
    submitters = {"MONDO:0005453": frozenset({"Ambry Genetics"})}
    truthful = _entry(id="MONDO:0005453")
    invented = _entry(id="MONDO:0005453", attributed_to="A Panel That Never Curated This")

    assert validate_scope_attribution([truthful], panels, submitters, _LOCATION) == []

    issues = validate_scope_attribution([invented], panels, submitters, _LOCATION)
    assert [i.code for i in issues] == ["SCP005"]
    assert issues[0].severity is Severity.ERROR
    assert "A Panel That Never Curated This" in issues[0].message
    # The message names what the mirror *does* say, so a curator can fix it
    # without going and reading the mirror themselves.
    assert "Congenital Heart Disease Gene Curation Expert Panel" in issues[0].message


def test_an_unreadable_mirror_skips_attribution_rather_than_failing_every_term() -> None:
    """One skip warning, not 59 errors naming the symptom and none the cause.

    The same reasoning as SCP000, REF000, SRC000 and ONT000. A term already
    reported absent by SCP001 is also skipped here, so one bad id costs one
    issue rather than two.
    """
    entry = _entry(id="MONDO:0005453")

    skipped = validate_scope_attribution([entry], None, None, _LOCATION)
    assert [(i.code, i.severity) for i in skipped] == [("SCP005", Severity.WARNING)]

    # Absent from the mirror entirely: SCP001's business, not this check's.
    assert validate_scope_attribution([entry], {}, {}, _LOCATION) == []


def test_a_commercial_submitters_term_cannot_be_relabelled_as_clingens_chd_panel() -> None:
    """The forgery the first version of SCP005 accepted, in the same commit that
    claimed to prevent it.

    `admitted_by` selects WHICH mirror the name must appear in. The original
    check merged ClinGen expert-panel names and GenCC submitter names into one
    set and asked only "did some authority use this term" -- so a curator could
    take a term whose sole warrant is a commercial laboratory's GenCC
    submission, write `admitted_by: clingen_chd_panel`, leave the laboratory in
    `attributed_to`, and pass. `_SCOPE_RULE` then tells every reader on 24 pages
    that a ClinGen expert panel scoped the disease. Found by adversarial review
    of #30.

    The third case is the narrower half of the same rule: `CLINGEN_CHD_PANEL`
    asserts the *chartered CHD panel* specifically, so another ClinGen panel
    claiming that member is also refused -- otherwise the strongest tier in the
    vocabulary would mean nothing more than "some ClinGen panel".
    """
    panels = {
        "MONDO:0005453": frozenset({CHD_EXPERT_PANEL}),
        "MONDO:0007732": frozenset({"Hearing Loss Gene Curation Expert Panel"}),
    }
    submitters = {"MONDO:0000119": frozenset({"Ambry Genetics"})}

    def check(entry: ScopeEntry) -> list[str]:
        return [i.code for i in validate_scope_attribution([entry], panels, submitters, _LOCATION)]

    # Honest, all three classes.
    assert (
        check(
            _entry(
                id="MONDO:0000119", admitted_by="gencc_submitter", attributed_to="Ambry Genetics"
            )
        )
        == []
    )
    assert (
        check(
            _entry(
                id="MONDO:0005453", admitted_by="clingen_chd_panel", attributed_to=CHD_EXPERT_PANEL
            )
        )
        == []
    )
    assert (
        check(
            _entry(
                id="MONDO:0007732",
                admitted_by="clingen_expert_panel",
                attributed_to="Hearing Loss Gene Curation Expert Panel",
            )
        )
        == []
    )

    # The forgery: a laboratory's term wearing the chartered panel's tier.
    assert check(
        _entry(id="MONDO:0000119", admitted_by="clingen_chd_panel", attributed_to="Ambry Genetics")
    ) == ["SCP005"]
    # And a real ClinGen panel claiming to be the CHD panel.
    assert check(
        _entry(
            id="MONDO:0007732",
            admitted_by="clingen_chd_panel",
            attributed_to="Hearing Loss Gene Curation Expert Panel",
        )
    ) == ["SCP005"]


# --- The Orphanet frequency warrant (issue #8 part 1, 2026-10-25) ---------------

_STRUCTURAL = frozenset({"HP:0001642", "HP:0001631"})
_PULMONIC = Annotation("HP:0001642", "Pulmonic stenosis", "Frequent (79-30%)")
_HCM = Annotation("HP:0001638", "Cardiomyopathy", "Very frequent (99-80%)")
_TABLES = OrphanetTables(
    annotations={"ORPHA:648": (_PULMONIC,), "ORPHA:79": (_HCM,), "ORPHA:3071": (_PULMONIC,)},
    names={"ORPHA:648": "Noonan syndrome", "ORPHA:79": "Fabry disease", "ORPHA:3071": "Costello"},
    omim_links={("ORPHA:648", "OMIM:610733"): ("BTNT", "Validated")},
    xrefs={
        "MONDO:0018997": frozenset({"ORPHA:648"}),
        "MONDO:0012547": frozenset({"OMIM:610733"}),
        "MONDO:0010526": frozenset({"ORPHA:79"}),
        "MONDO:0009026": frozenset({"ORPHA:3071"}),
    },
)


def _orphanet_entry(**overrides: object) -> ScopeEntry:
    fields: dict[str, object] = {
        "id": "MONDO:0018997",
        "label": "Noonan syndrome",
        "admitted_by": "orphanet_phenotype_frequency",
        "attributed_to": "Orphanet",
        "orphanet": "ORPHA:648",
    }
    return _entry(**{**fields, **overrides})


def test_attribution_leaves_an_orphanet_warrant_to_its_own_check() -> None:
    """SCP005 would otherwise look "Orphanet" up among ClinGen panels and refuse it."""
    panels = {"MONDO:0018997": frozenset({"RASopathy Gene Curation Expert Panel"})}
    issues = validate_scope_attribution([_orphanet_entry()], panels, {}, _LOCATION)
    assert issues == []


def test_a_nominated_syndrome_and_its_subtype_with_a_warrant_pass() -> None:
    subtype = _orphanet_entry(id="MONDO:0012547", label="Noonan syndrome 4")
    assert (
        validate_orphanet_warrants([_orphanet_entry(), subtype], _TABLES, _STRUCTURAL, _LOCATION)
        == []
    )


@pytest.mark.parametrize(
    ("entry", "reason"),
    [
        # The term does not reach the code it names: a wrong `orphanet` field.
        pytest.param(_orphanet_entry(orphanet="ORPHA:3071"), "does not reach", id="unreachable"),
        # It reaches its code, and Orphanet records only cardiomyopathy there.
        pytest.param(
            _orphanet_entry(id="MONDO:0010526", label="Fabry disease", orphanet="ORPHA:79"),
            "no structural heart phenotype",
            id="not-qualifying",
        ),
    ],
)
def test_a_nomination_orphanet_does_not_support_is_refused_saying_why(
    entry: ScopeEntry, reason: str
) -> None:
    """SCP006, one case per half of the rule, because the two are fixed differently."""
    issues = validate_orphanet_warrants([entry], _TABLES, _STRUCTURAL, _LOCATION)
    assert [issue.code for issue in issues] == ["SCP006"]
    assert issues[0].severity is Severity.ERROR
    assert reason in issues[0].message


@pytest.mark.parametrize(
    ("tables", "structural"),
    [pytest.param(None, _STRUCTURAL, id="mirrors"), pytest.param(_TABLES, None, id="hpo")],
)
def test_an_unreadable_input_is_one_error_and_one_skip_never_a_refusal_per_term(
    tables: OrphanetTables | None, structural: frozenset[str] | None
) -> None:
    """SCP008 is a skip warning, so it must arrive with the error that caused it (§4.41)."""
    entries = [_orphanet_entry(), _orphanet_entry(id="MONDO:0012547", label="Noonan syndrome 4")]
    issues = validate_orphanet_warrants(entries, tables, structural, _LOCATION)
    assert [(issue.code, issue.severity) for issue in issues] == [
        ("SCP006", Severity.ERROR),
        ("SCP008", Severity.WARNING),
    ]
    assert "2 Orphanet-scoped term(s)" in issues[0].message
    # And with nothing to check, nothing at all: a skip with no work is a lone warning.
    assert validate_orphanet_warrants([_entry()], tables, structural, _LOCATION) == []


def test_an_unnominated_term_that_would_qualify_is_listed_once_and_a_nominated_one_is_not() -> None:
    """SCP007: the review signal that makes a non-nomination a visible choice."""
    labels = {
        "MONDO:0018997": "Noonan syndrome",
        "MONDO:0009026": "Costello syndrome",
        "MONDO:0012547": "Noonan syndrome 4",
        "MONDO:0010526": "Fabry disease",
    }
    issues = orphanet_candidates({"MONDO:0018997"}, labels, _TABLES, _STRUCTURAL, _LOCATION)
    assert [(issue.code, issue.severity) for issue in issues] == [("SCP007", Severity.WARNING)]
    message = issues[0].message
    assert message.startswith("2 disease term(s) not in scope")
    # Sorted by term id; Fabry (cardiomyopathy only) and Noonan (nominated) absent.
    assert message.index("MONDO:0009026") < message.index("MONDO:0012547")
    assert "MONDO:0010526" not in message and "MONDO:0018997 " not in message
