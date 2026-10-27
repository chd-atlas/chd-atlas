# tests/unit/test_build_concordance.py
"""Whether independent datasets agree -- and what makes two datasets independent.

The committed corpus exercises few of these guards: its three studies cite
disjoint collection ids and form two families, PMID:34324492 and PMID:42230622
merged by the overlap the registry declares between DECIPHER and DDD (owner
decision 2026-10-27, I3), and PMID:40127276 alone. This read "three singleton
families" until then. The guards are driven by fixtures rather than left until
the day a study reuses a cohort, because a feature nobody exercised is a feature
that breaks the first time it matters -- the same reason `shared_cohorts` was
tested before a second study existed.
"""

from __future__ import annotations

import pytest

from chd_atlas.build.burden import BurdenRow, Overlap, shared_cohorts
from chd_atlas.build.concordance import (
    FamilyState,
    cohort_families,
    evidence_axes,
    family_state,
    gene_concordance,
    is_upward,
    survives_own_correction,
)
from chd_atlas.models.cohort import Cohort


def _row(**overrides: object) -> BurdenRow:
    """TAB2's all-cases loss-of-function row from PMID:42230622, as a base."""
    payload: dict[str, object] = {
        "study": "PMID:1",
        "gene": "HGNC:17075",
        "cohort_stratum": "all",
        "lesion_group": None,
        "variant_class": "snv_indel",
        "consequence_class": "lof",
        "origin": "any",
        "maf_max": 0.001,
        "count_unit": "individuals",
        "n_case_carriers": 6,
        "n_cases": 3876,
        "comparator": "control_cohort",
        "n_control_carriers": 0,
        "n_controls": 45082,
        "expected_count": None,
        "effect": 12.0,
        "effect_measure": "odds_ratio",
        "effect_bound": None,
        "ci_low": None,
        "ci_high": None,
        "pvalue": 2.45e-07,
        "pvalue_test": "fisher_exact",
        "pvalue_adjusted": None,
        "pvalue_adjustment": None,
        "case_cohorts": ("cnchd",),
        "control_cohorts": ("ukbb",),
        "method_note": None,
        "source": "audain2026_sd3",
    }
    payload.update(overrides)
    return BurdenRow(**payload)  # type: ignore[arg-type]


def test_two_studies_sharing_a_cohort_are_one_family() -> None:
    """The whole point of counting families rather than studies.

    If the strip counted studies, a reused collection would render as
    independent confirmation -- which is the inference `curation/cohorts.yaml`
    and `shared_cohorts` exist to prevent, arriving through the front door.
    """
    rows = [
        _row(study="PMID:1", case_cohorts=("cnchd", "ddd"), control_cohorts=("ukbb",)),
        _row(study="PMID:2", case_cohorts=("ddd", "pcgc"), control_cohorts=()),
        _row(study="PMID:3", case_cohorts=("decipher",), control_cohorts=("gain",)),
    ]

    assert cohort_families(rows, {}) == (
        frozenset({"PMID:1", "PMID:2"}),
        frozenset({"PMID:3"}),
    )


def test_a_declared_overlap_is_the_only_link_two_studies_need() -> None:
    """Owner decision 2026-10-27 (I3): two ids for partly the same people.

    PMID:34324492 cites `decipher` and PMID:42230622 `ddd`, and no id is
    shared, so they were two families on every page -- while DDD links its
    participants into DECIPHER, and an author of the first study answered that
    its cases do or likely do include DDD participants. `cohorts.yaml` now
    declares the two overlapping, and here the declaration is the **only**
    link: without the registry the same rows are three families, with it two,
    and `shared_cohorts` names the declared pair rather than a shared id. The
    third study, linked to nothing, stays apart.
    """
    rows = [
        _row(study="PMID:1", case_cohorts=("decipher",), control_cohorts=("gain",)),
        _row(study="PMID:2", case_cohorts=("ddd",), control_cohorts=("ukbb",)),
        _row(study="PMID:3", case_cohorts=("pcgc",), control_cohorts=()),
    ]
    registry = {
        "ddd": Cohort(id="ddd", name="DDD", description="d", overlaps_with=["decipher"]),
        "decipher": Cohort(id="decipher", name="DECIPHER", description="d", overlaps_with=["ddd"]),
    }

    assert cohort_families(rows, {}) == (
        frozenset({"PMID:1"}),
        frozenset({"PMID:2"}),
        frozenset({"PMID:3"}),
    )
    assert cohort_families(rows, registry) == (
        frozenset({"PMID:1", "PMID:2"}),
        frozenset({"PMID:3"}),
    )
    assert shared_cohorts(rows, registry) == {
        ("PMID:1", "PMID:2"): Overlap(shared=(), declared=(("decipher", "ddd"),))
    }


def test_a_transitive_chain_collapses_to_one_family() -> None:
    """A shares with B, B shares with C, A shares nothing with C.

    All three describe overlapping people, so all three are one dataset for the
    purpose of counting agreement. **This is what a pairwise read of
    `shared_cohorts` cannot say**: it reports two overlaps, and a reader
    counting distinct studies concludes there is some independent replication
    here. There is none.
    """
    rows = [
        _row(study="PMID:1", case_cohorts=("a",), control_cohorts=()),
        _row(study="PMID:2", case_cohorts=("a", "b"), control_cohorts=()),
        _row(study="PMID:3", case_cohorts=("b",), control_cohorts=()),
    ]

    assert cohort_families(rows, {}) == (frozenset({"PMID:1", "PMID:2", "PMID:3"}),)


def test_families_sort_deterministically_because_they_reach_published_bytes() -> None:
    """Derived from sets, so iteration order varies with `PYTHONHASHSEED`.

    `encode_json`'s `sort_keys` orders dict keys only; an array built from a set
    must be sorted by its caller. CLAUDE.md names this the likeliest way this
    build loses determinism, and the concordance object is published on 23
    browse rows and 23 bundles.
    """
    rows = [
        _row(study="PMID:9", case_cohorts=("z",), control_cohorts=()),
        _row(study="PMID:2", case_cohorts=("y",), control_cohorts=()),
        _row(study="PMID:5", case_cohorts=("x",), control_cohorts=()),
    ]

    assert cohort_families(rows, {}) == (
        frozenset({"PMID:2"}),
        frozenset({"PMID:5"}),
        frozenset({"PMID:9"}),
    )


@pytest.mark.parametrize(
    ("overrides", "upward"),
    [
        # An unbounded odds ratio is the strongest upward result in the data:
        # `allow_nan=False` refuses to publish the infinity, so the flag is it.
        pytest.param({"effect": None, "effect_bound": "unbounded_above"}, True, id="unbounded"),
        pytest.param({"effect": 6.53}, True, id="above-one"),
        # Below 1 is depletion. Every measure in `EFFECT_MEASURES` is a ratio
        # against 1, so one comparison covers all three.
        pytest.param({"effect": 0.47}, False, id="below-one"),
        pytest.param({"effect": 1.0}, False, id="exactly-one"),
        # PMID:34324492 publishes no effect at all -- PLINK's CNV association
        # test reports none -- so direction falls back to the published rates.
        pytest.param(
            {
                "effect": None,
                "effect_measure": None,
                "n_case_carriers": 12,
                "n_cases": 6746,
                "n_control_carriers": 0,
                "n_controls": 14024,
            },
            True,
            id="rates-enriched",
        ),
        pytest.param(
            {
                "effect": None,
                "effect_measure": None,
                "n_case_carriers": 1,
                "n_cases": 6746,
                "n_control_carriers": 40,
                "n_controls": 14024,
            },
            False,
            id="rates-depleted",
        ),
        # A mutation-model row compares its count against a modelled expectation.
        pytest.param(
            {
                "effect": None,
                "effect_measure": None,
                "comparator": "mutation_model",
                "n_control_carriers": None,
                "n_controls": None,
                "control_cohorts": (),
                "expected_count": 0.0276,
                "n_case_carriers": 1,
            },
            True,
            id="above-expectation",
        ),
        pytest.param(
            {
                "effect": None,
                "effect_measure": None,
                "comparator": "mutation_model",
                "n_control_carriers": None,
                "n_controls": None,
                "control_cohorts": (),
                "expected_count": 4.0,
                "n_case_carriers": 1,
            },
            False,
            id="below-expectation",
        ),
        # Neither an effect nor a usable comparator: not agreement. A row that
        # cannot state its direction must not be counted as supporting one.
        pytest.param(
            {
                "effect": None,
                "effect_measure": None,
                "n_control_carriers": None,
                "n_controls": None,
            },
            False,
            id="direction-unknowable",
        ),
    ],
)
def test_direction_is_read_from_what_the_study_published(
    overrides: dict[str, object], upward: bool
) -> None:
    """A significantly *depleted* row is not agreement.

    Three sources in order, because the studies do not publish the same columns:
    the unbounded flag, the effect against 1, then -- only where a study
    published no effect at all -- the rates it did publish. That last step
    compares two published numbers rather than computing one, and it is the
    closest this module comes to D12's line.
    """
    assert is_upward(_row(**overrides)) is upward


def test_a_family_that_did_not_test_the_gene_is_not_one_that_found_nothing() -> None:
    """**The state that exists to prevent harm.**

    KDM6A is ClinGen `definitive` -- it causes Kabuki syndrome -- and shows
    nothing in either dataset that tested it. Merging `NOT_TESTED` into
    `NO_ENRICHMENT` renders that as "0 of 3" beside a green definitive chip,
    which a clinical geneticist reads as the data contradicting the
    classification. Burden tests at these cohort sizes routinely detect nothing
    for genes with overwhelming family evidence.
    """
    assert family_state([]) is FamilyState.NOT_TESTED
    assert family_state([_row(pvalue=0.9)]) is FamilyState.NO_ENRICHMENT


def test_a_corrected_result_outranks_a_nominal_one() -> None:
    """Fill encodes the correction, so the state has to distinguish them.

    One glyph carries both statistics, which is what lets a reader see that a
    result is nominally significant *and* does not survive its own study's
    correction -- TAB2 in PMID:40127276, p 0.024 and q 0.222.
    """
    nominal = _row(pvalue=0.02, pvalue_adjusted=None)
    corrected = _row(pvalue=0.02, pvalue_adjusted=0.01, pvalue_adjustment="bonferroni")

    assert family_state([nominal]) is FamilyState.NOMINAL
    assert family_state([corrected]) is FamilyState.CORRECTED
    # The strongest row in the family decides it, whatever order they arrive in.
    assert family_state([nominal, corrected]) is FamilyState.CORRECTED
    assert family_state([corrected, nominal]) is FamilyState.CORRECTED


def test_a_corrected_p_above_alpha_does_not_promote_the_family() -> None:
    """TAB2's real shape in PMID:40127276: p 0.0242, q 0.222.

    A corrected p that is *present* is not a corrected p that is *significant*,
    and the whole value of the fill distinction is that it separates them.
    """
    row = _row(pvalue=0.0242, pvalue_adjusted=0.222, pvalue_adjustment="benjamini_hochberg")
    assert family_state([row]) is FamilyState.NOMINAL


@pytest.mark.parametrize(
    ("overrides", "survives"),
    [
        ({"pvalue": 0.001, "pvalue_adjusted": 0.01}, True),
        ({"pvalue": 0.001, "pvalue_adjusted": None}, False),
        ({"pvalue": 0.001, "pvalue_adjusted": 0.2}, False),
        ({"pvalue": 0.001, "pvalue_adjusted": 0.001, "consequence_class": "synonymous"}, False),
        ({"pvalue": 0.001, "pvalue_adjusted": 0.001, "effect": 0.5}, False),
        ({"pvalue": 0.05, "pvalue_adjusted": 0.01}, False),
    ],
    ids=[
        "supporting, corrected 0.01",
        "supporting, no correction published",
        "supporting, corrected 0.2",
        "synonymous negative control",
        "depleted, effect 0.5",
        "p at alpha, corrected 0.01",
    ],
)
def test_a_row_survives_its_own_correction_only_as_support_below_alpha_on_both(
    overrides: dict[str, object], survives: bool
) -> None:
    """The one rule the evidence matrix and the gene page's card both read (D56).

    Two renderers decide "survives its own study's correction", and two copies
    of a rule are two rules (CLAUDE.md section 4.28): this is the predicate, so
    a card that disagreed with its own matrix would have to stop calling it.
    Each case fails one half -- support (`_supports`: not synonymous, p below
    `ALPHA`, upward) or the study's own corrected p below `ALPHA` -- and the
    last pins the boundary: a p exactly at 0.05 is not below it, whatever the
    corrected value says.
    """
    row = _row(**overrides)
    assert survives_own_correction(row) is survives
    # And the matrix tier agrees row for row, which is the point of one predicate.
    assert (family_state([row]) is FamilyState.CORRECTED) is survives


def test_a_significant_synonymous_row_is_never_agreement() -> None:
    """It is the negative control.

    A significant synonymous row is a *warning* about that gene's comparison --
    the row that says whether to believe the two above it -- so counting it as
    support would invert the meaning of the only calibration check the table
    offers.
    """
    row = _row(consequence_class="synonymous", pvalue=1e-9)
    assert family_state([row]) is FamilyState.NO_ENRICHMENT


def test_a_significantly_depleted_row_is_not_agreement() -> None:
    """Significance is not support; direction is half the claim."""
    assert family_state([_row(pvalue=1e-9, effect=0.2)]) is FamilyState.NO_ENRICHMENT


def test_the_tally_counts_tested_families_and_reports_both_statistics() -> None:
    """ "3 of 3 tested, 1 corrected" -- both, because fill encodes the correction.

    TAB2's real shape. A tally reporting only `enriched` hides that two of the
    three do not survive correction; one reporting only `corrected` hides that
    PMID:42230622 publishes no correction on any row and so can never earn one.
    """
    families = (frozenset({"PMID:1"}), frozenset({"PMID:2"}), frozenset({"PMID:3"}))
    rows = [
        _row(study="PMID:1", pvalue=2.45e-07, pvalue_adjusted=None),
        _row(
            study="PMID:2",
            pvalue=0.0242,
            pvalue_adjusted=0.222,
            pvalue_adjustment="benjamini_hochberg",
        ),
        _row(
            study="PMID:3",
            pvalue=0.0001,
            pvalue_adjusted=0.0019,
            pvalue_adjustment="familywise_permutation",
        ),
    ]

    assert gene_concordance(rows, families) == {
        "tested": 3,
        "enriched": 3,
        "corrected": 1,
        "families": [
            {
                "studies": ["PMID:1"],
                "state": "nominal",
                "members": [{"study": "PMID:1", "state": "nominal"}],
            },
            {
                "studies": ["PMID:2"],
                "state": "nominal",
                "members": [{"study": "PMID:2", "state": "nominal"}],
            },
            {
                "studies": ["PMID:3"],
                "state": "corrected",
                "members": [{"study": "PMID:3", "state": "corrected"}],
            },
        ],
    }


def test_an_untested_family_is_published_and_excluded_from_the_denominator() -> None:
    """KDM6A's shape: two families tested, neither enriched, one not on the panel.

    The untested family is still published -- a consumer must be able to tell
    "not on that study's panel" from "that study looked and found nothing" --
    and it is kept out of `tested`, because that denominator is what stops the
    tally reading as a verdict on the gene.
    """
    families = (frozenset({"PMID:1"}), frozenset({"PMID:2"}), frozenset({"PMID:3"}))
    rows = [_row(study="PMID:1", pvalue=0.9), _row(study="PMID:2", pvalue=0.6)]

    result = gene_concordance(rows, families)

    assert result["tested"] == 2
    assert result["enriched"] == 0
    assert result["corrected"] == 0
    assert [entry["state"] for entry in result["families"]] == [  # type: ignore[union-attr]
        "no_enrichment",
        "no_enrichment",
        "not_tested",
    ]


def test_a_gene_no_study_reported_gets_a_not_tested_entry_per_family() -> None:
    """An object whose shape varies is a trap for a consumer reading a field off
    one gene and expecting it on the next -- the rule `burden_payload` follows.

    **This asserted the wrong condition under the right name.** It passed
    `gene_concordance([], ())`, which is an empty *corpus*, and concluded that a
    gene no study reported gets `families: []`. It does not: the corpus-wide
    family list is iterated whatever the gene has, so an unreported gene gets one
    `not_tested` entry per family. That is the whole point -- every gene shows the
    same slots, so a dataset that did not test this gene reads as an absence
    rather than as a shorter list. The old assertion let a false sentence stand in
    `docs/data-api.md` and in `manifest.py`'s changelog, both of which said
    `families` is `[]` for such a gene.
    """
    families = (frozenset({"PMID:1"}), frozenset({"PMID:2"}))

    assert gene_concordance([], families) == {
        "tested": 0,
        "enriched": 0,
        "corrected": 0,
        "families": [
            {
                "studies": ["PMID:1"],
                "state": "not_tested",
                "members": [{"study": "PMID:1", "state": "not_tested"}],
            },
            {
                "studies": ["PMID:2"],
                "state": "not_tested",
                "members": [{"study": "PMID:2", "state": "not_tested"}],
            },
        ],
    }

    # `[]` happens only when the corpus itself has no burden rows at all.
    assert gene_concordance([], ()) == {
        "tested": 0,
        "enriched": 0,
        "corrected": 0,
        "families": [],
    }


def test_a_family_of_two_studies_is_evaluated_over_both_their_rows() -> None:
    """One dataset, several papers: the family speaks once, with everything it has.

    Otherwise the merge in `cohort_families` would be cosmetic -- two studies
    would collapse into one row of the strip while only one of them was read.

    **And each member keeps its own state** (`members`, review 2026-10-27),
    because the strip's hover names every study, and it attributed the
    family's state to each: on a build of 5f7161b, 10 browse rows said
    PMID:42230622 survives its own correction, a study that publishes no
    corrected p, and 46 gave a state to a study with no row for the gene. Here
    the family is `corrected` on PMID:2's row alone; PMID:1 looked and found
    nothing, and PMID:3 has no row and is `not_tested`, the state a study that
    did not report on the gene keeps.
    """
    families = (frozenset({"PMID:1", "PMID:2", "PMID:3"}),)
    rows = [
        _row(study="PMID:1", pvalue=0.9),
        _row(study="PMID:2", pvalue=0.001, pvalue_adjusted=0.01, pvalue_adjustment="bonferroni"),
    ]

    result = gene_concordance(rows, families)

    assert result == {
        "tested": 1,
        "enriched": 1,
        "corrected": 1,
        "families": [
            {
                "studies": ["PMID:1", "PMID:2", "PMID:3"],
                "state": "corrected",
                "members": [
                    {"study": "PMID:1", "state": "no_enrichment"},
                    {"study": "PMID:2", "state": "corrected"},
                    {"study": "PMID:3", "state": "not_tested"},
                ],
            }
        ],
    }


def test_axes_come_from_the_whole_corpus_so_a_hole_looks_like_a_hole() -> None:
    """Derived corpus-wide, never per gene.

    Per gene, a gene tested only by case-control would render a one-column
    matrix, and the absence of a de novo analysis would be invisible -- exactly
    the "not tested is not a null result" confusion the four states exist for.
    Corpus-wide, every page has the same columns and an untested design renders
    as the empty cell it is.

    Sorted, because it comes from a set and orders the columns of published HTML.
    """
    rows = [
        _row(variant_class="snv_indel", comparator="control_cohort"),
        _row(
            variant_class="snv_indel",
            comparator="mutation_model",
            n_control_carriers=None,
            n_controls=None,
            control_cohorts=(),
            expected_count=0.1,
        ),
        _row(variant_class="cnv_deletion", comparator="control_cohort"),
        _row(variant_class="snv_indel", comparator="control_cohort"),
    ]

    assert evidence_axes(rows) == (
        ("cnv_deletion", "control_cohort"),
        ("snv_indel", "control_cohort"),
        ("snv_indel", "mutation_model"),
    )
