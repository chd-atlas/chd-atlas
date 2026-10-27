# src/chd_atlas/build/concordance.py
"""Whether independent datasets agree about a gene.

**This module computes no statistic.** It counts datasets, and reads the
direction of statistics the studies themselves published. That is arithmetic
over study identities rather than over their results, so design decision D33 --
no pooled statistic, ever -- stands: these cohorts overlap, and combining their
p-values would count the same children twice.

The distinction is the entire reason `cohort_families` exists rather than a
count of studies. A reader handed "3 of 3 datasets agree" will treat it as
replication, and it is only replication if the three describe different people.

Split from `burden.py` rather than added to it: that module loads and shapes
rows, this one derives a view over them, and the two change for different
reasons.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from enum import StrEnum
from typing import Final

from chd_atlas.build.burden import BurdenRow, shared_cohorts
from chd_atlas.build.emit import Json
from chd_atlas.models.cohort import Cohort

# The threshold both states are read against, and the one number in the rule
# `methods.html#rare-variants` states as the atlas's own (owner decision
# 2026-10-26, W2), which renders it from here. This said 0.05 was "the
# studies' convention and not the atlas's" until review (2026-10-27), which is
# not true of a raw p: PMID:42230622's block tells a reader to judge its raw p
# against its whole scan of 138,609 tests. Named rather than inlined so the two
# call sites below, and the page that states the rule, cannot drift apart.
ALPHA: Final = 0.05


class FamilyState(StrEnum):
    """What one cohort family says about one gene.

    `NOT_TESTED` is a distinct member and must never be folded into
    `NO_ENRICHMENT`; `family_state` records why.
    """

    CORRECTED = "corrected"
    NOMINAL = "nominal"
    NO_ENRICHMENT = "no_enrichment"
    NOT_TESTED = "not_tested"


def cohort_families(
    rows: Iterable[BurdenRow], cohorts: Mapping[str, Cohort]
) -> tuple[frozenset[str], ...]:
    """Studies grouped into sets that do not describe independent people.

    A family is a connected component of the graph whose nodes are studies and
    whose edges are shared sample collections, and collections the registry
    declares overlapping (`cohorts`, `Cohort.overlaps_with`). `shared_cohorts`
    computes the edges; this walks them.

    **Transitivity is why this is a graph walk and not a pairwise read of
    `shared_cohorts`.** If A shares a collection with B and B with C, all three
    describe overlapping people even where A and C share nothing directly. A
    reader told "two independent results" would be wrong, and the pairwise
    mapping alone cannot say so.

    Measured 2026-08-05 over the committed corpus: three studies, three
    singleton families, because they cite disjoint collections. That read
    "draw on disjoint collections" until owner decision 2026-10-27 (I3), and
    the collections were not disjoint: PMID:34324492's `decipher` cases may
    include DDD participants, and PMID:42230622 names `ddd`. Declared in the
    registry, the overlap merges the two into one family -- two families now,
    measured on a real build -- with no change here, as this was written to
    allow, and nobody is told that a reused cohort is independent evidence.

    Sorted on the sorted tuple of study ids. Components come out of a set walk,
    and `encode_json`'s `sort_keys` orders dict keys only, so an array derived
    from a set must be sorted by its caller or the build stops being
    reproducible.
    """
    materialised = list(rows)
    parent: dict[str, str] = {row.study: row.study for row in materialised}

    def find(study: str) -> str:
        # Path halving: the chains here are two or three long, so this is for
        # the property rather than the speed -- `find` must be idempotent for
        # the grouping pass below to be stable.
        while parent[study] != study:
            parent[study] = parent[parent[study]]
            study = parent[study]
        return study

    for left, right in shared_cohorts(materialised, cohorts):
        parent[find(left)] = find(right)

    grouped: dict[str, set[str]] = {}
    for study in parent:
        grouped.setdefault(find(study), set()).add(study)

    return tuple(
        frozenset(members) for members in sorted(grouped.values(), key=lambda m: sorted(m))
    )


def is_upward(row: BurdenRow) -> bool:
    """Whether the row's own statistic points to enrichment in cases.

    Read from what the study published, in three steps, because the three
    studies do not publish the same columns:

    1. `effect_bound == "unbounded_above"` -- an infinite odds ratio, which
       `allow_nan=False` refuses to publish as a number, so the flag *is* the
       result. 19 published rows.
    2. An effect against 1. Every member of `EFFECT_MEASURES` is a ratio, so one
       comparison covers odds ratios, rate ratios and enrichment ratios alike.
    3. Only where a study published no effect at all -- PMID:34324492's PLINK
       CNV test reports none -- the rates the study *did* publish.

    Step 3 compares two published numbers and authors nothing, but it is the
    closest this module comes to D12's line, so it is documented rather than
    left for a reviewer to infer. Note what it is not: it does not compute a
    ratio, a confidence interval or a p-value, and nothing it returns is
    published as a quantity -- only as a direction that decides a glyph.

    A row with neither an effect nor a usable comparator returns `False`. That
    is the safe default: it means "this does not count as agreement", and a row
    that cannot state its direction should not be counted as supporting.
    """
    if row.effect_bound == "unbounded_above":
        return True
    if row.effect is not None:
        return row.effect > 1
    if row.comparator == "mutation_model" and row.expected_count is not None:
        return row.n_case_carriers > row.expected_count
    if row.n_control_carriers is not None and row.n_controls:
        return (row.n_case_carriers / row.n_cases) > (row.n_control_carriers / row.n_controls)
    return False


def _supports(row: BurdenRow) -> bool:
    """Whether one row is evidence *for* the gene at `ALPHA`.

    `synonymous` is excluded because it is the study's negative control rather
    than a result. A significant synonymous row says that gene's comparison is
    poorly calibrated -- it is a warning about the two rows above it -- so
    counting it as support would invert the meaning of the only calibration
    check the table offers.
    """
    if row.consequence_class == "synonymous":
        return False
    if row.pvalue is None or row.pvalue >= ALPHA:
        return False
    return is_upward(row)


def survives_own_correction(row: BurdenRow) -> bool:
    """Whether a row is support for the gene that clears its own study's correction.

    The `CORRECTED` tier's rule, as one predicate, because three renderers
    decide by it -- the evidence matrix (its tier and the row its cell names),
    the gene page's card (D56) and the forest panels' fill for every row but a
    synonymous one (`pages._survived_correction`) -- and two copies of a rule
    are two rules (CLAUDE.md section 4.28). The forest drew its fill from the
    study's corrected p alone until 2026-10-26, so a depleted row with q below
    0.05 would have been drawn filled beside a card saying none survives. No
    row was in that state, measured that day over the 1,805 mirrored rows: the
    89 non-synonymous rows with a corrected p below 0.05 are the 89 this
    predicate passes. Never a correction this atlas computed: `pvalue_adjusted`
    is the study's own, and a row whose study published none does not survive
    one.
    """
    return _supports(row) and row.pvalue_adjusted is not None and row.pvalue_adjusted < ALPHA


def deciding_rows(rows: Iterable[BurdenRow]) -> list[BurdenRow]:
    """The rows that give a family its `family_state`, in the order given.

    Under `CORRECTED`, the rows that survive their own study's correction;
    under `NOMINAL`, the supporting rows, none of which does; empty under
    `NO_ENRICHMENT` and `NOT_TESTED`, where no row supports the gene.
    `family_state` reads its tier from this list, so a renderer naming one row
    for a state names a row that state is true of.

    **Why it exists.** The evidence matrix headlined the most significant row
    in a cell whatever that row's role, and GJA1's de novo cell in
    PMID:40127276 was drawn `nominal` by its missense row (p 0.0468) while
    naming the `damaging` union (p 0.0613, q 0.200), a row that is not
    nominally significant at all.
    """
    supporting = [row for row in rows if _supports(row)]
    return [row for row in supporting if survives_own_correction(row)] or supporting


def family_state(rows: Iterable[BurdenRow]) -> FamilyState:
    """The strongest thing one family says about one gene.

    **`NOT_TESTED` is not `NO_ENRICHMENT`, and the distinction prevents real
    harm.** KDM6A is ClinGen `definitive` -- it causes Kabuki syndrome -- and
    shows nothing in either of the two datasets that tested it. Collapsing the
    two states would render that as "0 of 3" beside a green definitive chip, and
    a clinical geneticist reads that as the data contradicting the
    classification. It is not: burden tests at these cohort sizes routinely
    detect nothing for genes with overwhelming family and functional evidence.

    So the tally counts *tested* families, the renderer gives the two states
    different glyphs, and the page carries a fixed caption saying that no
    enrichment here is not evidence against a gene.

    `CORRECTED` outranks `NOMINAL` because fill encodes the correction: one
    glyph carries both statistics, which is what lets a reader see that a result
    is nominally significant and does not survive its own study's correction.
    """
    materialised = list(rows)
    if not materialised:
        return FamilyState.NOT_TESTED

    deciding = deciding_rows(materialised)
    if not deciding:
        return FamilyState.NO_ENRICHMENT
    if survives_own_correction(deciding[0]):
        return FamilyState.CORRECTED
    return FamilyState.NOMINAL


def gene_concordance(
    rows: Iterable[BurdenRow], families: tuple[frozenset[str], ...]
) -> dict[str, Json]:
    """One gene's `independent_datasets` object, JSON-ready.

    `families` is passed in rather than derived here because it is a property of
    the *corpus*, not of the gene: a gene no study tested must still render the
    same family slots as every other gene, or a hole stops looking like a hole.

    Named `independent_datasets` in the payload rather than `concordance`. D12
    says the atlas authors no validity classification of its own, and a
    verdict-shaped key sitting beside a mirrored ClinGen `definitive` invites
    being read as a competing call. The name states what is counted.

    Every key is present on every gene, `families` empty where no study reported
    one -- the rule `burden_payload` follows, for the same reason: an object
    whose shape varies is a trap for a consumer that reads a field off one gene
    and expects it on the next.
    """
    by_study: dict[str, list[BurdenRow]] = {}
    for row in rows:
        by_study.setdefault(row.study, []).append(row)

    # The states are kept as enum members alongside the JSON entries rather than
    # read back out of them, so the tally is counted from the same values the
    # entries were built from and the two cannot disagree about one gene.
    states: list[FamilyState] = []
    entries: list[Json] = []
    for family in families:
        member_rows = [row for study in sorted(family) for row in by_study.get(study, ())]
        state = family_state(member_rows)
        states.append(state)
        entries.append({"studies": sorted(family), "state": state.value})

    return {
        # The denominator, and the reason the tally cannot read as a verdict.
        # A consumer must use this and never `len(families)`.
        "tested": sum(1 for state in states if state is not FamilyState.NOT_TESTED),
        "enriched": sum(
            1 for state in states if state in {FamilyState.CORRECTED, FamilyState.NOMINAL}
        ),
        "corrected": sum(1 for state in states if state is FamilyState.CORRECTED),
        "families": entries,
    }


def evidence_axes(rows: Iterable[BurdenRow]) -> tuple[tuple[str, str], ...]:
    """The distinct `(variant_class, comparator)` designs in the corpus.

    Passed the *whole* corpus, never one gene's rows. A per-gene derivation
    would give each page a different matrix, and a design nobody ran for that
    gene would vanish rather than render as the hole it is -- which is the same
    confusion between "not tested" and "tested and found nothing" that
    `FamilyState` has four members for.

    Measured 2026-08-05 over the published corpus: three axes -- SNV/indel
    case-control, SNV/indel de novo, CNV deletion case-control.

    **What this does NOT do, and an earlier version of this docstring claimed it
    did:** it does not render a column for a design *no* study ran. CNV de novo
    is such a design, and it produces no column at all rather than an empty one,
    because the axes come from the rows that exist. So the matrix shows holes at
    the level of "this family did not run a design *another* family ran", and is
    silent about a design nobody ran. That is a weaker guarantee than the
    original sentence promised, and it is the true one.

    Sorted: derived from a set, and it orders the columns of published HTML.
    """
    return tuple(sorted({(row.variant_class, row.comparator) for row in rows}))
