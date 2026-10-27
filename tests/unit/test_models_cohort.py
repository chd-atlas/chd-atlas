# tests/unit/test_models_cohort.py
"""The cohort registry, and the committed file it validates."""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from pydantic import ValidationError
from ruamel.yaml import YAML

from chd_atlas.build.burden import load_burden, shared_cohorts
from chd_atlas.corpus import Corpus, load_curation
from chd_atlas.models.cohort import CohortFile

REPO_ROOT = Path(__file__).resolve().parents[2]

_ENTRY = {
    "id": "ddd",
    "name": "Deciphering Developmental Disorders (DDD)",
    "description": "UK-wide study of children with severe undiagnosed developmental disorders.",
}


@pytest.mark.parametrize(
    ("mutation", "match"),
    [
        # The registry's own rule, not pydantic's: two entries under one id make
        # a `;`-joined burden cell ambiguous about which description it means.
        pytest.param({}, "duplicate cohort ids", id="duplicate-id"),
        # `_NonBlankStr` rather than a bare `min_length=1`. Whitespace renders as
        # nothing on a page while passing every presence check, so the entry
        # looks answered instead of missing -- the defect `models/scope.py` and
        # `validate/burden.py::_blank` both exist to catch on their own side.
        pytest.param({"description": "   "}, "at least 1 character", id="blank-description"),
        pytest.param({"name": "\t\n"}, "at least 1 character", id="blank-name"),
        # The grammar is shared with the `;`-joined mirror columns, so an id
        # carrying a separator would split into two ids nothing defines.
        pytest.param({"id": "ddd;ukbb"}, "should match pattern", id="separator-in-id"),
        pytest.param({"id": "DDD"}, "should match pattern", id="capital-in-id"),
    ],
)
def test_a_registry_entry_that_could_not_be_resolved_is_refused(
    mutation: dict[str, str], match: str
) -> None:
    """Each case makes some burden row's `case_cohorts` cell unresolvable."""
    second = {**_ENTRY, **mutation} if mutation else _ENTRY
    with pytest.raises(ValidationError, match=match):
        CohortFile.model_validate({"cohorts": [_ENTRY, second] if not mutation else [second]})


def test_an_unknown_field_is_refused_rather_than_silently_dropped() -> None:
    """`extra="forbid"`. A curator writing `descripton:` would otherwise get an
    entry with no description and no error, which is this project's
    characteristic failure in miniature.
    """
    with pytest.raises(ValidationError, match="Extra inputs"):
        CohortFile.model_validate({"cohorts": [{**_ENTRY, "role": "case"}]})


def test_the_committed_registry_parses_and_covers_every_cited_cohort() -> None:
    """A contract pin: `mirrors/burden.tsv` cites these ids by name.

    Named rather than merely counted, because the point of the file is *which*
    collections a row draws on -- a registry that still had the right number of
    entries after `ddd` was renamed would satisfy a count and break every row
    citing it.

    `taa_cases` is in this set and is deliberately not a CHD cohort: 777 of the
    7,958 cases in PMID:34324492 have sporadic thoracic aortic aneurysm, which
    `curation/chd_scope.yaml` puts outside this atlas. It is registered so that
    every row drawing on it says so, rather than being quietly folded into a
    case count labelled CHD.
    """
    path = REPO_ROOT / "curation" / "cohorts.yaml"
    parsed = CohortFile.model_validate(YAML(typ="safe").load(path.read_text(encoding="utf-8")))

    assert {cohort.id for cohort in parsed.cohorts} == {
        # PMID:42230622
        "cnchd",
        "ddd",
        "nottingham",
        "ukbb",
        # PMID:34324492
        "signature",
        "decipher",
        "isca",
        "taa_cases",
        # PMID:34324492's S2 Table sources identified by PMID with 50 or more
        # cases (owner decision R12, 2026-10-27).
        "pmid22318994_cases",
        "pmid22939634_cases",
        "pmid22912587_cases",
        "pmid22969434_cases",
        "pmid22929023_cases",
        "pmid26655555_cases",
        "pmid18535492_22349727_cases",
        "gain_controls",
        "wtccc_controls",
        "ottawa_controls",
        # PMID:40127276. `pcgc` carries the constraint that makes this registry
        # do its job prospectively rather than only descriptively: its entry
        # records that it subsumes Jin et al. 2017's trios, so a curator who
        # later adds that paper meets the double-count before writing the row.
        # `shared_cohorts` could never surface it, because the overlap is inside
        # one study's cohort list rather than between two.
        "pcgc",
        "gnomad_controls",
    }


@pytest.fixture(scope="module")
def corpus() -> Corpus:
    """The committed curation, loaded the way `runner.build_site` loads it."""
    loaded, issues = load_curation(REPO_ROOT)
    assert issues == [], issues
    return loaded


def test_the_gnomad_caveat_is_held_to_the_denominators_it_describes(corpus: Corpus) -> None:
    """`gnomad_controls` says its allele denominators sometimes exceed twice its samples.

    It read "run as high as 271,354, which is above twice that figure on 22 of
    the 90 rows mirrored here" until 2026-10-26: true when it was written and
    false once the mirror widened with the gene registry. Measured that day,
    123 PMID:40127276 rows cite this collection, all counting alleles, 36 of
    them above twice 133,743, the largest 271,396 -- which RAF1's page printed
    beside the sentence giving the maximum as 271,354. Reworded so no count in
    it can expire; this holds what is left to the mirror.

    The figure is read from the description rather than written here, so a
    reworded figure is checked against the data too. A mirror in which no
    denominator exceeds it fails, and that is the point: "the two are not
    reconcilable" then has lost its reason and the caveat must be rewritten.
    """
    description = {cohort.id: cohort for cohort in corpus.cohorts}["gnomad_controls"].description
    stated = re.search(r"states this control set as ([\d,]+) samples", description)
    assert stated, "the description no longer states a sample count to check against"
    assert "denominators it publishes are ALLELE counts, and some exceed twice that" in description

    rows = [
        row
        for rows in load_burden(REPO_ROOT).values()
        for row in rows
        if row.study == "PMID:40127276" and "gnomad_controls" in row.control_cohorts
    ]
    assert rows, "no row cites gnomad_controls; this guard would pass vacuously"
    assert {row.count_unit for row in rows} == {"alleles"}
    samples = int(stated.group(1).replace(",", ""))
    assert max(row.n_controls or 0 for row in rows) > 2 * samples


def test_taa_cases_says_no_scope_term_is_an_aortic_aneurysm_and_none_is(corpus: Corpus) -> None:
    """`taa_cases` tells a reader this atlas's scope "lists no aortic aneurysm".

    Nothing tied that to `curation/chd_scope.yaml` (review, 2026-10-26), and the
    list is one nomination away from making it false: SCP007 reports
    aneurysm-osteoarthritis syndrome as a term that would qualify on Orphanet's
    frequencies if a curator nominated it, beside Loeys-Dietz and Marfan
    syndromes. Measured that day: no scope term's label contains "aneurysm".

    Checked on the label, which is what the sentence is about: a Loeys-Dietz
    term, whose label does not say "aneurysm", would not make it false as
    written, though a curator admitting one should reread it. The phrase is
    asserted present rather than tested for, so rewording the description fails
    here until this check is re-pointed at whatever the new wording claims.
    """
    description = {cohort.id: cohort for cohort in corpus.cohorts}["taa_cases"].description
    assert "outside this atlas's scope, which lists no aortic aneurysm" in description
    assert corpus.chd_scope, "no scope terms loaded; this guard would pass vacuously"
    assert [term.label for term in corpus.chd_scope if "aneurysm" in term.label.lower()] == []


def test_ddd_and_decipher_are_declared_one_overlap_and_say_so(corpus: Corpus) -> None:
    """Two claims in curated text, each held to the data (owner decisions W4, W5, I3).

    `ddd` says "only one study cited today names it", a count that expires the
    day a second study's rows cite the collection: measured 2026-10-26, every
    mirrored row naming `ddd` is PMID:42230622's. It fails here then, and the
    sentence is reworded rather than left false.

    `ddd` and `decipher` both said the registry "cannot see" that DECIPHER's
    cases may include DDD participants, true while the two were separate ids
    and `shared_cohorts` compared ids. Owner decision 2026-10-27 (I3) declared
    them overlapping (`overlaps_with`), which made that clause false. It is
    dropped, not reworded: a page carrying both studies states the overlap in
    its own notice, and the descriptions keep the reason -- DDD links its
    participants into DECIPHER. This test held the old clause and is
    re-pointed at the reason and the declaration. The declaration is what
    links the real rows: a row citing `ddd` and a row citing `decipher`, from
    two studies, overlap through it and through nothing else.
    """
    described = {cohort.id: cohort.description for cohort in corpus.cohorts}
    assert "only one study cited today names it" in described["ddd"]
    assert "plans to curate" not in described["ddd"]
    rows = [row for rows in load_burden(REPO_ROOT).values() for row in rows]
    naming = {row.study for row in rows if "ddd" in (*row.case_cohorts, *row.control_cohorts)}
    assert len(naming) == 1, f"{sorted(naming)} name ddd; 'only one study' is false"

    for cohort in ("ddd", "decipher"):
        text = " ".join(described[cohort].split())
        assert "may include DDD participants." in text, cohort
        assert "cannot see that overlap" not in text, cohort
    registry = {str(cohort.id): cohort for cohort in corpus.cohorts}
    assert registry["ddd"].overlaps_with == ["decipher"]
    assert registry["decipher"].overlaps_with == ["ddd"]
    ddd = next(row for row in rows if "ddd" in row.case_cohorts)
    decipher = next(row for row in rows if "decipher" in row.case_cohorts)
    assert ddd.study != decipher.study
    pair = tuple(sorted((ddd.study, decipher.study)))
    overlap = shared_cohorts([ddd, decipher], registry)[pair]  # type: ignore[index]
    assert overlap.shared == ()
    assert shared_cohorts([ddd, decipher], {}) == {}


@pytest.mark.parametrize(
    ("ddd", "decipher", "match"),
    [
        pytest.param(["decipher"], [], "does not list ddd", id="one-sided"),
        pytest.param(["ddd"], [], "overlap with itself", id="self"),
        pytest.param(["gel"], [], "gel, which is not a cohort", id="unknown"),
        pytest.param(["decipher", "decipher"], ["ddd"], "twice", id="repeated"),
    ],
)
def test_a_declared_overlap_must_name_a_cohort_other_than_itself_on_both_ends(
    ddd: list[str], decipher: list[str], match: str
) -> None:
    """`overlaps_with` is refused at load unless it is known, not self, and symmetric (I3).

    It merges two studies into one cohort family on every page, so a typo'd
    id, a cohort overlapping itself, or a declaration on one end only is a
    malformed registry, reported when the file loads rather than published.
    """
    entries = [
        {**_ENTRY, "overlaps_with": ddd},
        {**_ENTRY, "id": "decipher", "name": "DECIPHER", "overlaps_with": decipher},
    ]
    with pytest.raises(ValidationError, match=match):
        CohortFile.model_validate({"cohorts": entries})
    CohortFile.model_validate(
        {
            "cohorts": [
                {**_ENTRY, "overlaps_with": ["decipher"]},
                {**_ENTRY, "id": "decipher", "name": "DECIPHER", "overlaps_with": ["ddd"]},
            ]
        }
    )
