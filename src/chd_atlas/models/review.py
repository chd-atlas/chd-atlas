# src/chd_atlas/models/review.py
"""An open question a curated record carries, stated rather than settled.

The owner's standing rule (answer of 2026-10-27, "from now on"): where a curated
record rests on a source that disagrees with the paper it cites, or on a source
that could not be read, the record says so. It carries a status naming which,
and a short note stating exactly what is unresolved. It never asserts either
side, never speculates about which is right, and is never dropped: dropping
the record would hide the question along with it.

**One vocabulary for every record kind that carries one**, so "pending review"
means the same thing on a cohort as on a scope term, and a consumer reads one
token set. Today that is `Cohort` and `ScopeEntry`; each declares the two
fields itself and calls `check_review` from its own validator.

Absent means no open question has been recorded, not that the record was
checked and found clean.
"""

from __future__ import annotations

from enum import StrEnum


class ReviewStatus(StrEnum):
    """Why a curated record is marked as carrying an open question.

    `PENDING_REVIEW`: a question awaits an expert -- the sources a record rests
    on disagree, or, as on the two scope terms the owner marked on 2026-10-27,
    a question about the record is left open for clinical review.
    `NOT_ACCESSIBLE`: a source it rests on could not be read, so what it says
    is unchecked. Neither is a verdict on the record.
    """

    PENDING_REVIEW = "pending_review"
    NOT_ACCESSIBLE = "not_accessible"


def check_review(record: str, status: ReviewStatus | None, note: str | None) -> None:
    """Refuse a status without its note, or a note without its status.

    A status alone tells a reader there is a question and not what it is, which
    invites the speculation the rule forbids; a note alone is an open question
    nothing marks as one, so it reads as part of the record. Raised from each
    model's validator, so the refusal arrives where the file loads (the build
    then refuses the corpus). A blank note is refused before this, by the
    fields' own non-blank type.
    """
    if status is not None and note is None:
        raise ValueError(
            f"{record}: review_status {status.value!r} needs a review_note stating what is "
            f"unresolved"
        )
    if status is None and note is not None:
        raise ValueError(
            f"{record}: a review_note needs a review_status saying which open question it records"
        )
