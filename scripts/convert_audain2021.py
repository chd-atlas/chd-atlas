# scripts/convert_audain2021.py
"""Convert the Audain 2021 CNV-deletion burden supplement into `mirrors/burden.tsv`.

Run locally by a curator, never by CI or the build -- design decision D6.

    .venv/bin/python scripts/convert_audain2021.py \\
        local/data/audain2021/s011.xlsx local/data/audain2021/s013.xlsx \\
        mirrors/burden.tsv

The first workbook is S6 Table, which every row's statistic is transcribed
from; the second is S8 Table, which gives one more number for the genes it
lists (see CHD-ONLY CARRIERS below).

Merges beside the rows other studies already contributed rather than replacing
the file; see `burden_mirror.merge_study`.

    Audain E, Wilsdon A, Breckpot J, et al. (2021) Integrative analysis of
    genomic variants reveals new associations of candidate haploinsufficient
    genes with congenital heart disease. PLoS Genet 17(7): e1009679.

Licence CC0 -- read off the article, which states it is "free of all copyright
... made available under the Creative Commons CC0 public domain dedication".
Not assumed from the publisher: PLOS uses CC BY on other titles.

THE ERRATUM WAS READ FIRST
--------------------------
PMID:34547032 (PLoS Genet 17(9): e1009809) corrects this paper. It reads, in
full: "The thirteenth author's name is spelled incorrectly. The correct name is:
D. Woodrow Benson." No data is corrected, so every number below stands. Checked
before a single figure was transcribed, because an erratum on a burden paper
could as easily have corrected the counts.

WHICH TABLE, AND WHY IT IS NOT THE OBVIOUS ONE
----------------------------------------------
S6 Table Sheet A: all 11,515 genes PLINK tested, with `case`, `control`, EMP1
and EMP2. S8 Table has something better in one respect -- it splits the case set
into CHD-only and TAA-only -- and something disqualifying in another: **it
contains only the 528 genes that reached significance.**

Measured 2026-08-05 against the registry: S6 covers 103 of the 154 in-scope
genes and 13 of the 23 published ones; S8 covers 12 and 4. Every one of S8's
four (GATA4, NKX2-5, NR2F2, TAB2) is significant by construction, while S6 also
carries the nulls that qualify them -- CHD7 at EMP1 0.0068 and EMP2 0.99, EVC2
at 0.22, ISL1 at 0.34. Mirroring S8 alone would have put a significant CNV
result on every gene page that showed one, which is publication bias built into
the atlas rather than merely reported by it.

So S6 it is, and with it the cost below.

THE CASE SET IS NOT PURELY CHD
------------------------------
S6's `case` column counts CHD **and** sporadic thoracic aortic aneurysm
probands together. The study says so in its abstract -- "probands with CHD
including cases with sporadic thoracic aortic aneurysm" -- and S2 Table puts
the TAA series (PMID 21092924) at 777 of the 7,958 contributing individuals,
about 10%. Thoracic aortic aneurysm is outside this atlas's scope:
`curation/chd_scope.yaml` lists no aortic aneurysm term.

Every row therefore names `taa_cases` in `case_cohorts`, whose registry entry
states the exclusion in as many words, and carries a `method_note` saying it on
the row itself. Nothing here silently relabels a CHD+TAA analysis as CHD.

WHY THE DE NOVO TABLE IS NOT MIRRORED
-------------------------------------
S10 Table (gene-based DNV analysis, 3,258 genes, 2,489 trios) fits this schema
technically -- it even carries `n_syn`, so the synonymous negative control would
work. It is deliberately left out.

**Its de novo variants are Jin 2017's and Sifrim 2016's.** The paper describes
itself as "a genome-wide integrative meta-analysis of published and publicly
available datasets of CNVs and DNVs", says "the data used in this study have
already been published", and cites Sifrim 2016 as reference [6] and Jin 2017 as
reference [7]. Those are the next two studies on this atlas's roadmap. Mirroring
S10 now and their own tables later would publish the same de novo variants twice
under two study ids, with the cohort columns unable to show it -- exactly the
double counting design decision D33 forbids, arriving through the one route the
overlap machinery cannot see.

The CNV part does not have that problem, and the difference is the point. S2
Table lists 189 case sources for 7,958 cases (measured 2026-10-26; this
paragraph said "~190 sources are mostly single-case reports plus four
collections" until then). The 11 registered as this study's case cohorts hold
6,211 of them (owner decision R12, 2026-10-26: the four it named before, 4,242,
and the seven CHD case sources S2 identifies by PMID with 50 or more cases,
1,969); of the other 178, 105 single-case sources hold 105, 69 sources of 2 to
49 cases hold 471, and 4 sources of 50 or more hold 1,171. None is a CNV study
this atlas carries, so aggregating them puts no deletion in the mirror twice.
The de novo part aggregates two studies this atlas intends to carry
individually. What the unregistered sources cost is visibility, not a double
count; see the comment on `CASE_COHORTS`.

When Jin 2017 and Sifrim 2016 land, their own de novo tables are the right
representation, and their cohorts will name the overlap between them properly.

CHD-ONLY CARRIERS (owner decision R2, 2026-10-27)
-------------------------------------------------
S8 Table repeats S6's test for the 528 genes that reached significance and
splits its case carriers by the case series they came from. Its first header
row names three bands, "CHD+TAA", "only CHD" and "only TAA" (columns 4, 8 and
12 of the workbook read on 2026-10-27), and its second reads `case`,
`control`, `p_value`, `p_value_adj` under each. The sheet's "Info" tab
defines `case` as "number of CNV cases", and its caption says the bands are
the 528 genes "combining both CHD and TAA cases (CHD+TAA), the contribution of
only CHD cases (only CHD) and the contribution of only TAA cases (only TAA)".
It is a gene-level table, one row per gene, like S6's Sheet A.

`n_chd_case_carriers` is the "only CHD" band's `case`: of the row's case
carriers, how many are CHD cases rather than aneurysm probands. It is
transcribed only where S8 lists the gene and is empty elsewhere; S8 is not a
source of rows. Before it is written, S8's "CHD+TAA" `case` must equal S6's
`case` for every gene S8 lists, or the two tables are not the same analysis
and the converter stops. Measured 2026-10-27: equal for all 528, and so are
the controls and the corrected p. S8 leaves the CHD-only `case` blank for 7
genes, none registered, and a blank is left blank rather than read as 0. The
workbook's sha256 is recorded as `audain2021_s8` in `mirrors/sources.yaml`.

**Both columns are located from the header rows, never from a fixed index**
(`_band_case`): the band's label once on the first row, then `case` once on
the second, between that label and the next. A sheet where either is not
found exactly once is refused. Nothing downstream could catch a misread band:
BUR003 holds `n_chd_case_carriers` to at most `n_case_carriers`, and every
band's `case` is within that bound -- measured 2026-10-27, "only CHD" plus
"only TAA" equals "CHD+TAA" for all 528 genes, a blank read as none -- so the
TAA-only or the combined count read in its place validated clean. The review
re-ran the converter both ways (2026-10-27): 0 errors each time, and read from
the TAA-only band, four cards said S8 counts 0 of the carriers among its CHD
cases.

DENOMINATORS
------------
Every row carries 6,746 cases and 14,024 controls, the set the article's
Methods put the permutation test on: "After filtering, 6,746 cases (3,929
harbouring CNV deletions) and 14,024 controls (12,585 harbouring CNV
deletions) remained for further analysis. A region-based permutation test
(using PLINK version 1.07 ...) was used on the filtered set to perform a
case-control association analysis." That is the rule applied to PMID:42230622,
whose abstract likewise quotes a larger set than its statistics used: a row's
denominators are the ones its statistic was computed on.

S1 Table ("Number of probands in the CNV cohort", by CNV type) is not that set.
The article cites it for the cohort as assembled, 7,958 cases and 14,082
controls, before describing the quality-control filters, and its deletion
column gives 7,107 cases and 14,040 controls. The mirror carried that column
until 2026-10-26, when this was resolved: the statistic was computed on the
filtered set the Methods state, so the filtered set is what the rows publish.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path
from typing import Final

from burden_mirror import merge_study
from convert_burden import gene_symbols, read_sheet, read_sheet_rows

STUDY: Final = "PMID:34324492"
SOURCE: Final = "audain2021_s6"

# S8 Table's header cells, quoted: two band labels from its first header row,
# and the two column names this reads from its second. The columns are found
# from these and never from a fixed index; see CHD-ONLY CARRIERS above.
S8_COMBINED: Final = "CHD+TAA"
S8_CHD_ONLY: Final = "only CHD"
S8_CASE: Final = "case"
S8_SYMBOL: Final = "symbol"

# S2/S3 Tables. 11 of S2's 189 case sources are registered, and they hold 6,211
# of its 7,958 cases (78%). Four were registered with the study on 2026-08-05,
# holding 4,242 (53%); this comment called the rest "~190 single-case PMID
# sources" until S2 was measured on 2026-10-26. Owner decision R12
# (2026-10-26) added the seven CHD case sources S2 identifies by PMID with 50
# or more cases: 742 + 432 + 340 + 178 + 123 + 91 + 63 = 1,969, read from S2's
# own counts on 2026-10-27 and recorded in each `curation/cohorts.yaml`
# description. (S2 identifies an eighth source of that size by PMID, its
# aneurysm series, which is `taa_cases`, one of the four.) Of the 178 still
# unregistered, measured from S2 that day, 105 are single cases (105
# individuals), 69 hold 2 to 49 (471 in all) and 4 hold 50 or more (1,171 in
# all). S2 labels 6 of the 178 by a collection name rather than a
# publication: those 4, which the owner decided are not registered, and 2 of
# the 69. The unregistered sources are named on no row, so `shared_cohorts`
# cannot report an overlap between one of them and a study curated later:
# check S2 before curating another CHD cohort. Alphabetical, as before.
CASE_COHORTS: Final = (
    "decipher;isca;pmid18535492_22349727_cases;pmid22318994_cases;pmid22912587_cases;"
    "pmid22929023_cases;pmid22939634_cases;pmid22969434_cases;pmid26655555_cases;"
    "signature;taa_cases"
)
CONTROL_COHORTS: Final = "gain_controls;ottawa_controls;wtccc_controls"

# The Methods' filtered set, on which the permutation test was run; see
# DENOMINATORS above. Not S1 Table's deletion column (7,107 and 14,040), which
# counts the cohort as assembled and which the mirror carried until 2026-10-26.
N_CASES: Final = 6746
N_CONTROLS: Final = 14024

# Renders unfolded directly under "N carriers / 6,746" on every gene page that
# carries this study. It read "777 sporadic thoracic aortic aneurysm probands
# (of 7,958 contributing individuals)" until 2026-10-26, which put all 777 in
# the case count printed above it. Measured 2026-10-26: S2 Table gives the
# series as 777 of 7,958; S8 Table's CHD+TAA case count equals S6's for all 528
# of its genes, and 144 of them have a TAA-only carrier, so the series is in
# the analysis. How many of the 777 are is stated nowhere read -- S1, S2, S3,
# S6, S8, S10, S11 and the article's text and table captions.
#
# The note deliberately names no case count for the row. `n_cases` prints
# directly above it, and the 7,958 the note does name is the assembled cohort
# S2 Table counts the 777 against, not the 6,746 the test ran on. How many of
# the 777 survived the filtering is stated nowhere read (above), which is why
# the note says so rather than putting a count on it.
METHOD_NOTE: Final = (
    "Cases combine congenital heart disease with a series of 777 sporadic "
    "thoracic aortic aneurysm probands among the study's 7,958 cases; the study "
    "does not report how many of the 777 are in this row's case set. The study "
    "reports the CHD-only contribution separately for significant genes. A "
    "carrier is an individual with a rare deletion (MAF below 0.01) overlapping "
    "the gene."
)


def _only(path: Path, cells: list[str], value: str, start: int, end: int, where: str) -> int:
    """The index of the one cell in `cells[start:end]` reading `value`; refuses none or several."""
    found = [index for index in range(start, min(end, len(cells))) if cells[index] == value]
    if len(found) != 1:
        raise ValueError(
            f"{path}: S8 Table's {where} reads {value!r} {len(found)} times, not once, so "
            f"the column it names cannot be located"
        )
    return found[0]


def _band_case(path: Path, groups: list[str], header: list[str], band: str) -> int:
    """The index of `case` under `band`: its label once on the first header row,
    then `case` once on the second, from that label up to the next one.

    A band label sits in the first cell of the band it heads and the cells
    after it are blank, so the band ends where the next label begins, or with
    the row.
    """
    start = _only(path, groups, band, 0, len(groups), "first header row")
    end = next((index for index in range(start + 1, len(groups)) if groups[index]), len(header))
    return _only(path, header, S8_CASE, start, end, f"second header row under {band!r}")


def read_chd_only(path: Path) -> dict[str, tuple[str, str]]:
    """S8 Table as symbol -> (its "CHD+TAA" `case`, its "only CHD" `case`).

    Every column is located from the two header rows (`_band_case`); refuses
    a sheet where one cannot be, and a gene listed twice.
    """
    rows = read_sheet_rows(path, "S8 Table")
    groups, header = ([cell.strip() for cell in row] for row in [*rows, [], []][:2])
    symbol_at = _only(path, header, S8_SYMBOL, 0, len(header), "second header row")
    combined_at = _band_case(path, groups, header, S8_COMBINED)
    chd_only_at = _band_case(path, groups, header, S8_CHD_ONLY)
    split: dict[str, tuple[str, str]] = {}
    for row in rows[2:]:
        # Padded: a row read from the sheet ends at its last cell, which may
        # fall short of the columns this reads.
        cells = row + [""] * (len(header) - len(row))
        symbol = cells[symbol_at].strip()
        if not symbol:
            continue
        if symbol in split:
            raise ValueError(f"{path}: S8 Table lists {symbol} twice")
        split[symbol] = (cells[combined_at].strip(), cells[chd_only_at].strip())
    return split


def convert(
    source: Path, chd_only: Path, symbols: dict[str, str]
) -> tuple[list[dict[str, str]], list[str]]:
    """Build the mirror rows, and the report lines describing what was dropped.

    `source` is S6 Table and `chd_only` S8 Table; see CHD-ONLY CARRIERS above.
    """
    records = read_sheet(source, "Sheet A")
    split = read_chd_only(chd_only)
    cases = {record["symbol"]: record["case"] for record in records}
    disagree = sorted(
        symbol for symbol, (combined, _) in split.items() if cases.get(symbol) != combined
    )
    if disagree:
        raise ValueError(
            f"S8 Table's CHD+TAA case count differs from S6 Table's, or S6 lists no such "
            f"gene, for {len(disagree)} gene(s): {', '.join(disagree)}"
        )
    rows: list[dict[str, str]] = []
    seen: set[str] = set()

    for record in records:
        symbol = record["symbol"]
        seen.add(symbol)
        hgnc = symbols.get(symbol)
        if hgnc is None:
            continue
        rows.append(
            {
                "study": STUDY,
                "gene": hgnc,
                # Not `all`: `CohortStratum.ALL` means a study's undivided case
                # set, and this one is undivided in a way that includes a
                # disease the atlas excludes. The stratum stays `all` because
                # that is what the study analysed; `taa_cases` and the method
                # note are what say the set is not purely CHD.
                "cohort_stratum": "all",
                "lesion_group": "",
                "variant_class": "cnv_deletion",
                # A deletion removes the gene; the study frames it as
                # haploinsufficiency throughout, and PLINK counts individuals
                # with an overlapping deletion rather than a coding consequence.
                "consequence_class": "lof",
                "origin": "any",
                # "Rare CNV deletions (MAF 0.01)" -- S15 Table's own title.
                "maf_max": "0.01",
                # PLINK's CNV association counts individuals carrying a deletion
                # that overlaps the gene, not alleles and not deletions: a
                # proband with two overlapping deletions is one case.
                "count_unit": "individuals",
                "n_case_carriers": record["case"],
                "n_cases": str(N_CASES),
                # S8 Table's "only CHD" `case`, for a gene S8 lists; empty
                # otherwise, and where S8 leaves the cell blank.
                "n_chd_case_carriers": split.get(symbol, ("", ""))[1],
                "comparator": "control_cohort",
                "n_control_carriers": record["control"],
                "n_controls": str(N_CONTROLS),
                "expected_count": "",
                # PLINK's CNV association test reports no effect size and no
                # confidence interval, so none is invented. `_effect` renders an
                # em dash, and the counts beside it carry the finding.
                "effect": "",
                "effect_measure": "",
                "effect_bound": "",
                "ci_low": "",
                "ci_high": "",
                "pvalue": _number(record["p_value"]),
                "pvalue_test": "permutation",
                # EMP2, PLINK's max(T) family-wise corrected p over the same
                # permutations that produced EMP1. This study publishes its own
                # correction, which PMID:42230622 does not.
                "pvalue_adjusted": _number(record["p_value_adj"]),
                "pvalue_adjustment": "familywise_permutation",
                "case_cohorts": CASE_COHORTS,
                "control_cohorts": CONTROL_COHORTS,
                "method_note": METHOD_NOTE,
                "source": SOURCE,
            }
        )

    kept = {row["gene"] for row in rows}
    missing = sorted(symbol for symbol, hgnc in symbols.items() if hgnc not in kept)
    split_rows = [row for row in rows if row["n_chd_case_carriers"]]
    smaller = sum(
        1 for row in split_rows if int(row["n_chd_case_carriers"]) < int(row["n_case_carriers"])
    )
    blank = sorted(symbol for symbol, (_, chd) in split.items() if not chd)
    report = [
        f"{len(rows)} row(s) for {len(kept)} of the {len(symbols)} registered genes",
        f"{len(seen)} genes tested in the supplement; {len(seen) - len(kept)} out of scope",
        f"{len(split_rows)} row(s) carry S8 Table's CHD-only case carriers, of {len(split)} "
        f"genes S8 lists; on {smaller} of them it is below the row's case carriers",
        f"S8 leaves the CHD-only case count blank for {len(blank)} gene(s): {', '.join(blank)}",
    ]
    if missing:
        report.append(
            f"WARNING: {len(missing)} registered gene(s) were not tested by this study, "
            f"so they will show no CNV burden from it: {', '.join(missing)}"
        )
    return rows, report


def _number(value: str) -> str:
    """A published float at its shortest exact spelling; see `convert_burden`."""
    return repr(float(value))


def main(argv: list[str]) -> int:
    if not 3 <= len(argv) <= 4:
        print(__doc__, file=sys.stderr)
        return 2
    source, chd_only, out = Path(argv[0]), Path(argv[1]), Path(argv[2])
    root = Path(argv[3]) if len(argv) == 4 else Path.cwd()

    rows, report = convert(source, chd_only, gene_symbols(root))
    kept, written = merge_study(out, STUDY, rows)

    for line in report:
        print(line)
    print(f"wrote {written} row(s) for {STUDY}, keeping {kept} row(s) from other studies")
    for path in (source, chd_only):
        print(f"sha256 of {path.name}: {hashlib.sha256(path.read_bytes()).hexdigest()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
