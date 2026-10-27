# tests/unit/test_convert_audain2021.py
"""The Audain 2021 converter's row note and denominators, on the row it builds and in the mirror.

The note renders unfolded under "N carriers / 6,746" on every gene page carrying
this study, so it has to be true of that row. Until 2026-10-26 it read "777
sporadic thoracic aortic aneurysm probands (of 7,958 contributing individuals)",
which put all 777 in the case count printed above it -- and no test read the
note at all, because nothing had ever driven this converter (CLAUDE.md §4.23).

The denominators are the 6,746 cases and 14,024 controls the article's methods
say the permutation test was run on, after filtering. The rows carried S1
Table's deletion column, 7,107 and 14,040, until 2026-10-26; each test below
fails on that pair, one on the converter and one on a mirror not re-run.
"""

from __future__ import annotations

import sys
from pathlib import Path

import polars as pl
import pytest

from chd_atlas.tables import BURDEN, read_table

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import convert_audain2021  # noqa: E402

STUDY = "PMID:34324492"

# Written out rather than imported, so a change to `METHOD_NOTE` has to be made
# here too, deliberately. It names no case count of its own: the converter's
# comment on `METHOD_NOTE` records why.
_METHOD_NOTE = (
    "Cases combine congenital heart disease with a series of 777 sporadic "
    "thoracic aortic aneurysm probands among the study's 7,958 cases; the study "
    "does not report how many of the 777 are in this row's case set. The study "
    "reports the CHD-only contribution separately for significant genes. A "
    "carrier is an individual with a rare deletion (MAF below 0.01) overlapping "
    "the gene."
)


def test_the_row_the_converter_builds_carries_the_current_note_and_denominators(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Driven through `convert` on two synthetic S6 records: one registered gene, one not."""
    record = {"symbol": "TBX1", "case": "120", "control": "3", "p_value": "1E-4"}
    records = [{**record, "p_value_adj": "2E-4"}, {**record, "symbol": "NOT-REGISTERED"}]
    sheets: list[str] = []

    def read_sheet(path: Path, name: str) -> list[dict[str, str]]:
        sheets.append(name)
        return records

    monkeypatch.setattr(convert_audain2021, "read_sheet", read_sheet)
    rows, _ = convert_audain2021.convert(Path("s011.xlsx"), {"TBX1": "HGNC:11592"})

    assert sheets == ["Sheet A"]
    assert [row["gene"] for row in rows] == ["HGNC:11592"]
    (row,) = rows
    assert row["method_note"] == _METHOD_NOTE
    # Every case source registered for this study, as a literal: the four S2
    # collections, then the seven S2 series identified by PMID with 50 or more
    # cases (owner decision R12, 2026-10-27). The committed mirror reads the
    # same, but cannot see a converter that dropped one and was not re-run.
    assert row["case_cohorts"].split(";") == [
        "decipher",
        "isca",
        "pmid18535492_22349727_cases",
        "pmid22318994_cases",
        "pmid22912587_cases",
        "pmid22929023_cases",
        "pmid22939634_cases",
        "pmid22969434_cases",
        "pmid26655555_cases",
        "signature",
        "taa_cases",
    ]
    # Literals, not `N_CASES`/`N_CONTROLS`: a test reading the constants
    # compares the converter to itself.
    assert (row["n_cases"], row["n_controls"]) == ("6746", "14024")


def test_every_committed_row_carries_the_current_note_and_denominators() -> None:
    """Reads the committed mirror, so a converter edited without being re-run fails here."""
    frame, issues = read_table(REPO_ROOT / "mirrors" / "burden.tsv", BURDEN)
    assert frame is not None, issues
    mine = frame.filter(pl.col("study") == STUDY)
    assert mine.height, f"no rows for {STUDY}; did the merge drop them?"
    assert mine["method_note"].unique().to_list() == [_METHOD_NOTE]
    assert mine["n_cases"].unique().to_list() == [6746]
    assert mine["n_controls"].unique().to_list() == [14024]
