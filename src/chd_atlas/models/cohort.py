# src/chd_atlas/models/cohort.py
"""The sample collections a published burden analysis drew on.

This registry exists for one reason: **sample reuse between studies is invisible
in a table of per-study rows, and a reader will meta-analyse them anyway.**

The atlas deliberately publishes no meta-statistic (D33) precisely because the
CHD literature reuses cohorts -- DDD samples appear in more than one of the
papers this atlas cites, so a pooled p-value would count the same children
twice. But declining to compute one does not stop a reader eyeballing four rows
and forming the same wrong impression by hand.

Naming the cohorts turns that from an invisible trap into a displayed fact: two
rows that share `ddd` can be shown to share samples. That is the only place this
registry is used, and it is why `role` is not a field here -- a collection is
cases in one study and could be controls in another, so which it was is a
property of the burden row, not of the cohort.
"""

from __future__ import annotations

from typing import Annotated

from pydantic import BaseModel, ConfigDict, Field, StringConstraints, model_validator

from chd_atlas.duplicates import duplicates
from chd_atlas.identifiers import CohortId
from chd_atlas.models.review import ReviewStatus, check_review

# Same rule and rationale as `models/scope.py`: whitespace passes a bare
# `min_length` and renders as nothing, which looks answered rather than missing.
_NonBlankStr = Annotated[str, StringConstraints(strip_whitespace=True, min_length=1)]


class Cohort(BaseModel):
    """One named sample collection, as a burden row's `case_cohorts` cites it.

    `description` is where the caveats that qualify every number drawn from this
    collection belong -- an adult control cohort compared against children
    enrolled in infancy carries a survivorship bias that no column can express,
    and the paper's own limitations section is the right source for it.
    """

    model_config = ConfigDict(extra="forbid")

    id: CohortId
    name: _NonBlankStr
    description: _NonBlankStr
    url: str | None = None
    # Collections this one shares people with though no row names both
    # (owner decision 2026-10-27, I3). `shared_cohorts` sees an overlap only
    # where two studies cite one id; DDD and DECIPHER are two ids for
    # partly the same children, so two studies citing one each were counted
    # as independent datasets. Declared on both entries (`CohortFile`
    # checks it), so neither record reads as independent on its own.
    overlaps_with: list[CohortId] = Field(default_factory=list)
    # A question about this record that the atlas leaves open (owner's
    # standing rule, 2026-10-27; `models/review.py`): both or neither, the
    # note stating only what is unresolved. Published in `cohorts.json` and
    # rendered after the description wherever it renders.
    review_status: ReviewStatus | None = None
    review_note: _NonBlankStr | None = None

    @model_validator(mode="after")
    def a_review_status_carries_its_note(self) -> Cohort:
        check_review(self.id, self.review_status, self.review_note)
        return self


class CohortFile(BaseModel):
    """Top level of ``curation/cohorts.yaml``."""

    model_config = ConfigDict(extra="forbid")

    cohorts: list[Cohort] = Field(min_length=1)

    @model_validator(mode="after")
    def ids_are_unique(self) -> CohortFile:
        found = duplicates(cohort.id for cohort in self.cohorts)
        if found:
            raise ValueError(f"duplicate cohort ids: {found}")
        return self

    @model_validator(mode="after")
    def overlaps_are_known_symmetric_and_not_self(self) -> CohortFile:
        """`overlaps_with` names registered cohorts, never itself, and both ends declare it.

        Symmetric rather than one-sided, so the record a reader opens for
        either collection says so; an overlap declared on one end only is the
        curator's half-finished edit, refused rather than completed silently.
        """
        known = {cohort.id: set(cohort.overlaps_with) for cohort in self.cohorts}
        for cohort in self.cohorts:
            repeated = duplicates(cohort.overlaps_with)
            if repeated:
                raise ValueError(f"{cohort.id} lists {repeated} twice in overlaps_with")
            for other in cohort.overlaps_with:
                if other == cohort.id:
                    raise ValueError(f"{cohort.id} declares an overlap with itself")
                if other not in known:
                    raise ValueError(f"{cohort.id} overlaps_with {other}, which is not a cohort")
                if cohort.id not in known[other]:
                    raise ValueError(
                        f"{cohort.id} overlaps_with {other}, but {other} does not list "
                        f"{cohort.id}; declare the overlap on both"
                    )
        return self
