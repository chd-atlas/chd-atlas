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

`n_chd_case_carriers` is read from S8 Table (owner decision R2, 2026-10-27),
from the column its header rows locate. No two cells of the fixture's split
gene's S8 row are alike, so a converter reading any wrong column writes a
number this fixture does not expect. The fixture's genes are synthetic, and so
are their numbers: none is a fact about a real gene.
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
from burden_mirror import ORDER  # noqa: E402

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


# S8 Table's two header rows as the workbook has them, written out rather than
# built from the converter's constants, then one gene per line: symbol, chr,
# start, end, and each band's case, control, p, corrected p. Note the trailing
# space in "significant CHD ", as the sheet spells it.
_S8_GROUPS = ["", "", "", "", "CHD+TAA", "", "", "", "only CHD", "", "", "", "only TAA"]
_S8_HEADER = [
    "symbol", "chr", "start", "end",
    "case", "control", "p_value", "p_value_adj",
    "case", "control", "p_value", "p_value_adj",
    "case", "control", "p_value", "p_value_adj",
    "significant CHD ", "significant TAA", "LOEUF", "DNM", "*Expression RNA-Seq", "Info",
]  # fmt: skip
# Synthetic genes, synthetic numbers. SPLIT's combined `case` is its S6 `case`
# (120), split 117 CHD-only and 3 TAA-only, and no two of its cells are alike.
_SPLIT_COMBINED = ["120", "21", "1E-4", "2E-4"]
_SPLIT_CHD_ONLY = ["117", "17", "3E-4", "4E-4"]
_SPLIT_TAA_ONLY = ["3", "13", "5E-4", "6E-4"]
_S8_SPLIT = ["SPLIT", "22", "11", "12", *_SPLIT_COMBINED, *_SPLIT_CHD_ONLY, *_SPLIT_TAA_ONLY]
# S8 leaves a band blank where it reports nothing for it.
_S8_BLANK = ["BLANK", "8", "31", "32", "29", "25", "7E-4", "8E-4", "", "", "", "", "29"]
_SYMBOLS = {"SPLIT": "HGNC:900001", "BLANK": "HGNC:900002", "UNLISTED": "HGNC:900003"}


def _drive(
    monkeypatch: pytest.MonkeyPatch, s8: list[list[str]], symbols: dict[str, str]
) -> list[dict[str, str]]:
    """`convert` on synthetic S6 records: SPLIT and BLANK, which S8 lists,
    UNLISTED, which it does not, and one gene no registry holds."""
    record = {"symbol": "SPLIT", "case": "120", "control": "21", "p_value": "1E-4"}
    records = [
        {**record, "p_value_adj": "2E-4"},
        {**record, "symbol": "BLANK", "case": "29", "p_value_adj": "2E-4"},
        {**record, "symbol": "UNLISTED", "case": "30", "p_value_adj": "2E-4"},
        {**record, "symbol": "UNREGISTERED", "p_value_adj": "2E-4"},
    ]
    sheets: list[str] = []

    def read_sheet(path: Path, name: str) -> list[dict[str, str]]:
        sheets.append(name)
        return records

    def read_sheet_rows(path: Path, name: str) -> list[list[str]]:
        sheets.append(name)
        return s8

    monkeypatch.setattr(convert_audain2021, "read_sheet", read_sheet)
    monkeypatch.setattr(convert_audain2021, "read_sheet_rows", read_sheet_rows)
    rows, _ = convert_audain2021.convert(Path("s011.xlsx"), Path("s013.xlsx"), symbols)
    assert sheets == ["Sheet A", "S8 Table"]
    return rows


def test_the_row_the_converter_builds_carries_the_current_note_and_denominators(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Driven through `convert` on synthetic S6 and S8 records.

    `n_chd_case_carriers` is S8's "only CHD" `case` (117 for SPLIT, not the
    combined 120 or the TAA-only 3), empty where S8 leaves that cell blank
    (BLANK) and where S8 does not list the gene (UNLISTED). Every row carries
    exactly the mirror's columns, in the schema's order.
    """
    rows = _drive(monkeypatch, [_S8_GROUPS, _S8_HEADER, _S8_SPLIT, _S8_BLANK], _SYMBOLS)

    assert [(row["gene"], row["n_chd_case_carriers"]) for row in rows] == [
        ("HGNC:900001", "117"),
        ("HGNC:900002", ""),
        ("HGNC:900003", ""),
    ]
    assert all(tuple(row) == ORDER == BURDEN.column_names for row in rows)
    row = rows[0]
    assert row["method_note"] == _METHOD_NOTE
    # Every case source registered for this study, as a literal: the four
    # registered with the study, then the seven CHD series S2 identifies by
    # PMID with 50 or more cases (owner decision R12, 2026-10-26). The
    # committed mirror reads the same, but cannot see a converter that dropped
    # one and was not re-run.
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


@pytest.mark.parametrize(
    ("groups", "header", "split"),
    [
        # The TAA-only band before the CHD-only one: an index fixed at the
        # workbook's column 8 would publish the TAA-only 3.
        pytest.param(
            [*_S8_GROUPS[:8], "only TAA", "", "", "", "only CHD"],
            _S8_HEADER,
            ["SPLIT", "22", "11", "12", *_SPLIT_COMBINED, *_SPLIT_TAA_ONLY, *_SPLIT_CHD_ONLY],
            id="bands-reordered",
        ),
        # `case` second in every band: a column taken from the band's label
        # alone, without reading the second header row, is each band's
        # `control`.
        pytest.param(
            _S8_GROUPS,
            ["symbol", "chr", "start", "end", *("control", "case", "p_value", "p_value_adj") * 3],
            [
                "SPLIT", "22", "11", "12",
                "21", "120", "1E-4", "2E-4",
                "17", "117", "3E-4", "4E-4",
                "13", "3", "5E-4", "6E-4",
            ],
            id="case-second-in-each-band",
        ),
    ],
)  # fmt: skip
def test_the_chd_only_count_is_read_from_the_column_its_headers_name(
    monkeypatch: pytest.MonkeyPatch, groups: list[str], header: list[str], split: list[str]
) -> None:
    """Located from S8's header rows, so a moved column is still the right one.

    BUR003 cannot catch a misread band: every band's `case` is at most the
    row's case carriers (see the converter's CHD-ONLY CARRIERS), so the wrong
    one validates clean. The combined band is located the same way and still
    checked against S6, which these layouts pass only if it is found too.
    """
    rows = _drive(monkeypatch, [groups, header, split], {"SPLIT": "HGNC:900001"})
    assert [row["n_chd_case_carriers"] for row in rows] == ["117"]


@pytest.mark.parametrize(
    ("s8", "match"),
    [
        pytest.param(
            [[*_S8_GROUPS[:8], "CHD only", "", "", "", "only TAA"], _S8_HEADER, _S8_SPLIT],
            "first header row reads 'only CHD' 0 times",
            id="band-unlabelled",
        ),
        pytest.param(
            [[*_S8_GROUPS[:8], "only CHD", "", "", "", "only CHD"], _S8_HEADER, _S8_SPLIT],
            "first header row reads 'only CHD' 2 times",
            id="band-labelled-twice",
        ),
        # The band's columns renamed: nothing under it reads `case`.
        pytest.param(
            [_S8_GROUPS, [*_S8_HEADER[:8], "cases", *_S8_HEADER[9:]], _S8_SPLIT],
            "second header row under 'only CHD' reads 'case' 0 times",
            id="case-unnamed",
        ),
        # With the next band unlabelled, the CHD-only band runs to the end of
        # the row and holds two `case` columns, and neither can be chosen.
        pytest.param(
            [_S8_GROUPS[:12], _S8_HEADER, _S8_SPLIT],
            "second header row under 'only CHD' reads 'case' 2 times",
            id="next-band-unlabelled",
        ),
        # S8's combined count is not S6's: the two tables are not one analysis.
        pytest.param(
            [_S8_GROUPS, _S8_HEADER, [*_S8_SPLIT[:4], "121", *_S8_SPLIT[5:]]],
            "differs from S6",
            id="combined-is-not-s6",
        ),
        # A gene S6 does not list cannot be checked against it.
        pytest.param(
            [_S8_GROUPS, _S8_HEADER, ["ELSEWHERE", *_S8_SPLIT[1:]]],
            "S6 lists no such gene, for 1 gene\\(s\\): ELSEWHERE",
            id="gene-not-in-s6",
        ),
        pytest.param(
            [_S8_GROUPS, _S8_HEADER, _S8_SPLIT, _S8_SPLIT],
            "S8 Table lists SPLIT twice",
            id="gene-listed-twice",
        ),
    ],
)
def test_the_converter_stops_on_an_s8_table_it_cannot_read_as_s6s_split(
    monkeypatch: pytest.MonkeyPatch, s8: list[list[str]], match: str
) -> None:
    """Refusals rather than reports: each would otherwise write a wrong count, or none."""
    with pytest.raises(ValueError, match=match):
        _drive(monkeypatch, s8, {"SPLIT": "HGNC:900001"})


def test_every_committed_row_carries_the_current_note_and_denominators() -> None:
    """Reads the committed mirror, so a converter edited without being re-run fails here.

    And the S8 split (owner decision R2, 2026-10-27), on the 12 registered
    genes S8 lists, as the converter transcribed them on 2026-10-27; the other
    109 rows carry none. Only this study's rows carry the column, which is what
    lets `pages._chd_only` name "S8 Table": a second study's split fails here
    first, and that sentence must then name its own study's table.
    """
    frame, issues = read_table(REPO_ROOT / "mirrors" / "burden.tsv", BURDEN)
    assert frame is not None, issues
    mine = frame.filter(pl.col("study") == STUDY)
    assert mine.height, f"no rows for {STUDY}; did the merge drop them?"
    assert mine["method_note"].unique().to_list() == [_METHOD_NOTE]
    assert mine["n_cases"].unique().to_list() == [6746]
    assert mine["n_controls"].unique().to_list() == [14024]

    split = frame.filter(pl.col("n_chd_case_carriers").is_not_null())
    assert split["study"].unique().to_list() == [STUDY]
    assert {
        gene: (carriers, chd)
        for gene, carriers, chd in split.select(
            "gene", "n_case_carriers", "n_chd_case_carriers"
        ).iter_rows()
    } == {
        "HGNC:11592": (120, 120),  # TBX1
        "HGNC:17075": (12, 12),  # TAB2
        "HGNC:17985": (10, 10),  # ROBO4
        "HGNC:18196": (21, 21),  # SOX7
        "HGNC:21316": (13, 12),  # ANKRD11, unpublished
        "HGNC:2488": (10, 10),  # NKX2-5
        "HGNC:3227": (12, 12),  # EFNB2
        "HGNC:3327": (30, 30),  # ELN
        "HGNC:3488": (16, 16),  # ETS1
        "HGNC:4173": (29, 29),  # GATA4
        "HGNC:7976": (10, 10),  # NR2F2
        "HGNC:8569": (10, 10),  # PCSK6, unpublished
    }
