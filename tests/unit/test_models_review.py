# tests/unit/test_models_review.py
"""A curated record's open question: a status and its note, together or not at all."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pytest
from pydantic import BaseModel, ValidationError

from chd_atlas.build.burden import build_cohorts
from chd_atlas.build.emit import Emitter
from chd_atlas.build.paths import COHORTS, SCOPE
from chd_atlas.build.scope import build_scope
from chd_atlas.models.cohort import Cohort, CohortFile
from chd_atlas.models.review import ReviewStatus
from chd_atlas.models.scope import ChdScopeFile, ScopeEntry

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


def _emit_cohort(record: BaseModel, emitter: Emitter) -> None:
    assert isinstance(record, Cohort)
    build_cohorts([record], emitter)


def _emit_term(record: BaseModel, emitter: Emitter) -> None:
    assert isinstance(record, ScopeEntry)
    build_scope([record], {}, emitter)


@pytest.mark.parametrize(
    ("model", "entry", "emit", "published"),
    [
        pytest.param(Cohort, _COHORT, _emit_cohort, COHORTS, id="cohorts.json"),
        pytest.param(ScopeEntry, _TERM, _emit_term, SCOPE, id="scope.json"),
    ],
)
@pytest.mark.parametrize(
    ("fields", "match"),
    [
        pytest.param(
            {"review_status": ReviewStatus.PENDING_REVIEW, "review_note": None},
            "needs a review_note",
            id="no-note",
        ),
        pytest.param(
            {"review_status": None, "review_note": "S2 and the abstract differ."},
            "needs a review_status",
            id="no-status",
        ),
    ],
)
def test_each_emitter_refuses_a_record_that_bypassed_the_check_before_writing_it(
    tmp_path: Path,
    model: type[BaseModel],
    entry: dict[str, str],
    emit: Callable[[BaseModel, Emitter], None],
    published: str,
    fields: dict[str, object],
    match: str,
) -> None:
    """Both or neither, again where the JSON is written (review, 2026-10-27, N2).

    The test above holds the refusal at load, and `render.review_status_note`
    raises on a page; but a record built without validation reached
    `build_cohorts` and `build_scope` and published `"pending_review"` beside
    a `null` note. `model_copy(update=...)` skips validation, so the record
    here is a valid one with the pair broken afterwards -- every other field
    is what the emitter expects, and nothing but the check can refuse it.
    Each emitter raises before writing, so no file is left behind either.
    """
    record = model.model_validate(entry).model_copy(update=fields)
    with pytest.raises(ValueError, match=match):
        emit(record, Emitter(root=tmp_path))
    assert not (tmp_path / published).exists()
