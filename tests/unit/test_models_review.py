# tests/unit/test_models_review.py
"""A curated record's open question: a status and its note, together or not at all."""

from __future__ import annotations

import pytest
from pydantic import BaseModel, ValidationError

from chd_atlas.models.cohort import CohortFile
from chd_atlas.models.scope import ChdScopeFile

_COHORT = {
    "id": "pmid22969434_cases",
    "name": "CHD case series of PMID:22969434",
    "description": "178 individuals with congenital heart disease.",
}
_TERM = {
    "id": "MONDO:1010098",
    "label": "TFAP2B-related congenital heart disease spectrum disorder",
    "reason": "Patent ductus arteriosus is definitional.",
    "admitted_by": "clingen_expert_panel",
    "attributed_to": "Syndromic Disorders Gene Curation Expert Panel",
    "admitted_on": "2026-08-03",
}

# Both record kinds that carry one, each through its own file model: each
# model calls the shared check from its own validator, so a call dropped from
# either is a refusal one kind of record no longer makes.
_KINDS = [
    pytest.param(CohortFile, "cohorts", _COHORT, id="cohort"),
    pytest.param(ChdScopeFile, "diseases", _TERM, id="scope-term"),
]


@pytest.mark.parametrize(("model", "key", "entry"), _KINDS)
@pytest.mark.parametrize(
    ("fields", "match"),
    [
        pytest.param({"review_status": "pending_review"}, "needs a review_note", id="no-note"),
        pytest.param({"review_note": "S2 and the abstract differ."}, "needs a review_status",
                     id="no-status"),
        pytest.param({"review_status": "pending_review", "review_note": "  \n"},
                     "at least 1 character", id="blank-note"),
        pytest.param({"review_status": "disputed", "review_note": "x"}, "Input should be",
                     id="unknown-status"),
    ],
)  # fmt: skip
def test_a_review_status_is_refused_at_load_without_its_note_and_a_note_without_its_status(
    model: type[BaseModel],
    key: str,
    entry: dict[str, str],
    fields: dict[str, str],
    match: str,
) -> None:
    """The owner's standing rule (2026-10-27): the status says which open question, the
    note what is unresolved, and neither means anything alone.

    A status with no note tells a reader there is a question and not what it
    is, which invites the speculation the rule forbids; a note with no status
    reads as part of the record. A note of whitespace renders as nothing and
    looks answered. The vocabulary is closed, so a curator cannot invent a
    third status that no consumer reads. Both together load, and neither is
    the default: absent means no question recorded.
    """
    with pytest.raises(ValidationError, match=match):
        model.model_validate({key: [{**entry, **fields}]})

    together = {"review_status": "not_accessible", "review_note": "The source could not be read."}
    loaded = getattr(model.model_validate({key: [{**entry, **together}]}), key)[0]
    assert (loaded.review_status, loaded.review_note) == (
        "not_accessible",
        "The source could not be read.",
    )
    absent = getattr(model.model_validate({key: [entry]}), key)[0]
    assert absent.review_status is None and absent.review_note is None
