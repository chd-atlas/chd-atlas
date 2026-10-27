# scripts/convert_audain2021.py
"""Convert the Audain 2021 CNV-deletion burden supplement into `mirrors/burden.tsv`.

Run locally by a curator, never by CI or the build -- design decision D6.

    .venv/bin/python scripts/convert_audain2021.py \\
        local/data/audain2021/s011.xlsx mirrors/burden.tsv

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
6,211 of them (2026-10-27, owner decision R12: the four it named before, 4,242,
and the seven S2 identifies by PMID with 50 or more cases, 1,969); of the
other 178, 105 single-case sources hold 105, 69 sources of 2 to 49 cases hold
471, and 4 sources of 50 or more hold 1,171. None is a CNV study this atlas
carries, so aggregating them puts no deletion in the mirror twice. The de novo
part aggregates two studies this atlas intends to carry individually. What the
unregistered sources cost is visibility, not a double count; see the comment
on `CASE_COHORTS`.

When Jin 2017 and Sifrim 2016 land, their own de novo tables are the right
representation, and their cohorts will name the overlap between them properly.

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
from convert_burden import gene_symbols, read_sheet

STUDY: Final = "PMID:34324492"
SOURCE: Final = "audain2021_s6"

# S2/S3 Tables. 11 of S2's 189 case sources are registered, and they hold 6,211
# of its 7,958 cases (78%). This comment called the rest "~190 single-case PMID
# sources" until 2026-10-26, when four were registered, holding 4,242 (53%).
# Owner decision R12 (2026-10-27) added the seven S2 identifies by PMID with 50
# or more cases: 742 + 432 + 340 + 178 + 123 + 91 + 63 = 1,969, read from S2's
# own counts and recorded in each `curation/cohorts.yaml` description. Of the
# 178 still unregistered, measured from S2 that day, 105 are single cases (105
# individuals), 69 hold 2 to 49 (471 in all) and 4 hold 50 or more (1,171 in
# all); those 4 are the S2 sources labelled by a collection name rather than a
# publication, which the owner decided are not registered. The unregistered
# sources are named on no row, so `shared_cohorts` cannot report an overlap
# between one of them and a study curated later: check S2 before curating
# another CHD cohort. Alphabetical, as before.
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


def convert(source: Path, symbols: dict[str, str]) -> tuple[list[dict[str, str]], list[str]]:
    """Build the mirror rows, and the report lines describing what was dropped."""
    records = read_sheet(source, "Sheet A")
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
    report = [
        f"{len(rows)} row(s) for {len(kept)} of the {len(symbols)} registered genes",
        f"{len(seen)} genes tested in the supplement; {len(seen) - len(kept)} out of scope",
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
    if not 2 <= len(argv) <= 3:
        print(__doc__, file=sys.stderr)
        return 2
    source, out = Path(argv[0]), Path(argv[1])
    root = Path(argv[2]) if len(argv) == 3 else Path.cwd()

    rows, report = convert(source, gene_symbols(root))
    kept, written = merge_study(out, STUDY, rows)

    for line in report:
        print(line)
    print(f"wrote {written} row(s) for {STUDY}, keeping {kept} row(s) from other studies")
    print(f"sha256 of {source.name}: {hashlib.sha256(source.read_bytes()).hexdigest()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
