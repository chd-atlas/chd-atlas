# tests/unit/test_build_methods.py
"""`methods.html`: what every gene page assumes, written once (D50).

Every sentence on this page was moved off a gene page, and the move is where it
could go false: a sentence's truth conditions travel with its position
(CLAUDE.md section 4.27). So each one is pinned here against a **literal**,
sliced to its own section, never against the constant that produced it
(section 4.38) -- and the phrases the move retired are asserted absent from the
whole page, because a pasted old constant is how they would come back.
"""

from __future__ import annotations

import html
import re

import pytest

from chd_atlas.build.emit import Emitter
from chd_atlas.build.methods import build_methods
from chd_atlas.build.paths import METHODS_PAGE


@pytest.fixture(scope="module")
def page(tmp_path_factory: pytest.TempPathFactory) -> str:
    """The page as built, with no cohorts: `#scope` reads none of them."""
    root = tmp_path_factory.mktemp("methods")
    build_methods(Emitter(root=root), cohorts={})
    return (root / METHODS_PAGE).read_text(encoding="utf-8")


def _section(page: str, anchor: str) -> str:
    """One section: from its `<h2 id=...>` to the next `<h2`, or the end of the page."""
    start = page.index(f'<h2 id="{anchor}">')
    end = page.find("<h2", start + 1)
    return page[start:] if end == -1 else page[start:end]


def test_the_page_sits_at_the_site_root_under_its_own_title(page: str) -> None:
    """P1: one page at the root, named for what a reader comes to it for.

    `root=""` is what makes its own nav resolve from the site root; a page
    rendered with `"../"` would point every visitor's "Genes" link above the
    deploy root.
    """
    assert METHODS_PAGE == "methods.html"
    assert "<title>How to read this atlas — CHD Atlas</title>" in page
    assert "<h1>How to read this atlas</h1>" in page
    assert '<a href="genes/index.html">Genes</a>' in page
    assert '<a href="methods.html">How to read this atlas</a>' in page
    assert '<h2 id="scope">' in page


def test_scope_states_the_rule_that_admits_a_gene_in_full(page: str) -> None:
    """The rule both page kinds stated in full until 2026-10-26, and now link to (D50).

    Moved from `pages._SCOPE_RULE`, which was pinned by
    `test_both_page_kinds_state_the_rule_that_admits_a_gene_to_this_atlas` and,
    for its Orphanet clause, by
    `test_an_orphanet_scoped_syndrome_says_why_it_counts_as_chd_above_both_groups`.
    Those tests' recorded decisions travel with it:

    * **Both warrants, or the sentence describes a gate the site does not
      have.** Stating only the ClinGen half would tell a reader that every gene
      here was graded by an expert panel; 21 of the 119 were not (measured
      2026-10-26), and they are the ones whose provenance most needs stating.
    * **The negative is the load-bearing half.** The rule used to say a gene is
      definitive "for a disease in this atlas's CHD scope", which told 24
      published pages that this atlas decides what counts as congenital heart
      disease. It does not; scope is a selection among terms external
      authorities already use. Without the negative a future edit could
      restore the old sentence beside the new one and satisfy every positive.
    * **Sliced to the section, because page-wide the fragments could not
      fail.** The grade key's `ungraded` gloss carries "two or more Gene
      Curation Coalition submitters independently assert it" in the same
      section's `<dl>`, so the slice is the rule's own `<p>`: rewriting the
      rule to state a **one**-submitter gate survived all 796 tests when the
      check was page-wide (measured 2026-08-06, CLAUDE.md section 4.19).

    Its last clause said the disease "is shown beside the gene", true on a gene
    page and on the gene index and unanchored here, so it now names both
    places -- and only for a gene a ClinGen panel graded: the 21 genes admitted
    on submitter agreement have no panel-named disease to show (2026-10-26).
    The full literal pins the wording; the fragments beside it are the
    claims any rewording must keep, so a deliberate edit to the literal is not
    also a silent edit to the gate it describes.
    """
    section = _section(page, "scope")
    rule = re.search(r'<p class="scope-rule">.*?</p>', section, re.S)
    assert rule, "#scope carries no scope-rule element"
    sentence = rule.group(0)
    assert sentence == (
        '<p class="scope-rule">A gene is published here when an external authority already '
        "treats it as a congenital heart disease gene: either a <strong>ClinGen expert "
        "panel</strong> classifies it <strong>Limited or better</strong> for a disease that "
        "an external authority treats as congenital heart disease (or in which "
        "<strong>Orphanet records a structural heart malformation in at least 30% of "
        "cases</strong>), or <strong>two or more "
        "Gene Curation Coalition submitters</strong> independently assert it and no ClinGen "
        "panel disputes it. "
        "<strong>No disease is in scope on this atlas's own judgement</strong>, and neither "
        "is any gene: every scope term names the authority it came from "
        "(<code>scope.json</code>, from <code>curation/chd_scope.yaml</code>), and every gene "
        "page names the authority that admitted the gene. A classification is not the same "
        "as definitive for congenital heart disease: for many of these genes the disease "
        "graded is a syndrome of which a cardiac lesion is one feature, and for a gene a "
        "ClinGen panel graded, the disease the panel named is shown beside the gene, on "
        "its page and in the gene index.</p>"
    )
    for claim in (
        "Limited or better",
        "two or more Gene Curation Coalition submitters",
        "no ClinGen panel disputes it",
        "an external authority treats as congenital heart disease",
        "No disease is in scope on this atlas's own judgement",
        "not the same as definitive for congenital heart disease",
        "for a disease that an external authority treats as congenital heart disease (or in "
        "which <strong>Orphanet records a structural heart malformation in at least 30% of "
        "cases</strong>)",
        "(<code>scope.json</code>, from <code>curation/chd_scope.yaml</code>)",
    ):
        assert claim in sentence, f"the scope rule no longer says {claim!r}"


def test_scope_glosses_every_rung_as_one_label_chip_and_gloss(page: str) -> None:
    """The grade key, open, under the rule (P3). Moved from both page kinds (D50).

    Moved from `test_both_page_kinds_gloss_every_grade_they_can_render`, with its
    recorded decisions: a colour ramp with no key is a ranking the reader has to
    guess at, and since 2026-08-06 the site renders four rungs and an ungraded
    state -- `limited`, a panel saying the case is not yet made, is the
    plurality. Each rung's gloss is checked by a distinctive phrase rather than
    in full, so this fails when a rung stops being explained and not when a
    comma moves.

    **The label, the chip class and the gloss are asserted as one triple, and
    that is the whole point of this test.** It checked only that seven phrases
    appeared *somewhere on the page* until 2026-08-06, which says nothing about
    which rung each is attached to. Measured that day: painting `moderate` with
    `strong`'s chip -- inverting the colour ramp mid-ladder -- **survived all 796
    tests**, as did swapping two glosses. Scoped to the key's `<dl>`, because
    these words also occur in the scope rule above it (section 4.19).

    **Rendered open, not folded** (P3). On a gene page the key was a collapsed
    `<details>` because it was not what the page was for; here it is.
    """
    section = _section(page, "scope")
    assert "<details" not in section, "the grade key is folded on the page that exists for it"
    key = re.search(r'<div class="grade-key">.*?</div>', section, re.S)
    assert key, "#scope carries no grade key"
    legend = re.search(r"<dl>(.*?)</dl>", key.group(0), re.S)
    assert legend, "the grade key carries no rows"

    # (label a reader sees, chip class that colours it, phrase unique to its gloss)
    rungs = (
        ("definitive", "chip-definitive", "upheld over time"),
        ("strong", "chip-strong", "without the years of replication"),
        ("moderate", "chip-moderate", "several unrelated probands"),
        ("limited", "chip-limited", "not yet made"),
        ("no known association", "chip-warn", "found no evidence of a relationship"),
        (
            "not classified",
            "chip-ungraded",
            "two or more Gene Curation Coalition submitters independently assert it",
        ),
    )
    pairs = re.findall(r"<dt>(.*?)</dt><dd>(.*?)</dd>", legend.group(1), re.S)
    assert len(pairs) == len(rungs), f"the key glosses {len(pairs)} rungs, expected {len(rungs)}"
    for (label, chip, gloss), (dt, dd) in zip(rungs, pairs, strict=True):
        assert chip in dt, f"{label} is not chipped {chip} -- got {dt}"
        assert f">{label}<" in dt, f"expected the label {label}, got {dt}"
        assert gloss in dd, f"{label} is not glossed by {gloss!r}"


def test_scope_states_where_a_dispute_is_listed_without_pointing_below(page: str) -> None:
    """The dispute paragraph, reworded for a page with no validity table on it (P3).

    Moved from `test_the_grade_key_no_longer_denies_an_in_scope_dispute_it_now_displays`.
    The sentence went false once already, when Costello and CFC syndrome entered
    scope on 2026-10-25: it read "No gene here was disputed or refuted by an
    expert panel for a disease an external authority treats as congenital heart
    disease", and ClinGen's RASopathy panel disputes six published genes for
    those two syndromes -- PTPN11's page carried the denial beside two
    `Disputed` rows.

    **Reworded again on the move (2026-10-26, D50).** It said "the record is in
    the table below", which was false on the gene index (no validity table) and
    is false here; it now names the gene's validity table. "This page is not a
    complete account" became "a gene page", for the same reason.

    **And its condition named one path of two.** "Where that disease is one an
    external authority treats as congenital heart disease" excluded the
    Orphanet path -- and Costello and cardiofaciocutaneous syndrome, the
    diseases all six disputed genes (BRAF, KRAS, MAP2K1, PTPN11, RAF1, SOS1) are
    disputed for, are in scope only through it (`admitted_by:
    orphanet_phenotype_frequency`). It now reads "in scope here, by either path
    above", which `test_the_key_follows_the_rule_its_dispute_sentence_points_back_to`
    keeps true.
    """
    section = _section(page, "scope")
    key = re.search(r'<div class="grade-key"><p>(.*?)</p>', section, re.S)
    assert key, "#scope carries no grade-key paragraph"
    assert key.group(1) == (
        "Every classification on this site is mirrored from the authority that made it. "
        "The rung names are "
        '<a href="https://clinicalgenome.org/curation-activities/gene-disease-validity/">'
        "ClinGen&#x27;s</a>, glossed here from their Gene-Disease Validity SOP; Gene "
        "Curation Coalition submitters grade on the same ladder, and also use "
        "<em>Supportive</em>, which asserts an association without grading its evidence "
        "and which this atlas maps to no rung. "
        "<strong>No gene here is published on a disputed or refuted grade.</strong> A panel "
        "may still have disputed the same gene for a different disease. Where that disease "
        "is in scope here, by either path above, the record is in that gene&#x27;s validity "
        "table, marked <em>Disputed</em> or <em>Refuted</em>, and the gene carries a "
        "<em>conflicting evidence</em> flag. Where it is not, the record "
        "is not shown: a gene page lists only records for diseases in this atlas&#x27;s "
        "remit, so it is not a complete account of what ClinGen has said about the gene."
    )


@pytest.mark.parametrize(
    "retired",
    [
        # Told 24 pages this atlas decides what counts as CHD (2026-08-06).
        "this atlas's CHD scope",
        # False since Costello and CFC entered scope (2026-10-25).
        "No gene here was disputed or refuted",
        # False for 5 genes with an out-of-scope `Disputed` (2026-08-06).
        "not published here at all",
        # There is no table on this page (D50).
        "the table below",
        # A one-submitter gate is not the gate (2026-08-06).
        "one or more Gene Curation",
    ],
)
def test_no_sentence_the_move_retired_comes_back_here(page: str, retired: str) -> None:
    """Whole page, both spellings of an apostrophe (section 4.35).

    A positive literal fails when a sentence is deleted; these fail when a
    deleted one comes back, and coming back -- an old constant pasted onto the
    page it was moved to -- is the failure nobody watches for. Checked on the
    unescaped text too, because the page spells one apostrophe `'` and another
    `&#x27;`, and a retired phrase must not pass by changing which.
    """
    assert retired not in page
    assert retired not in html.unescape(page)


def test_the_key_follows_the_rule_its_dispute_sentence_points_back_to(page: str) -> None:
    """The rule renders above the key, because the key says "by either path above".

    The key's dispute sentence names no path; it points back at the rule, which
    names both. Rendered above the rule -- or anywhere the rule is not above it
    -- "above" points at nothing, the same defect "the table below" was.
    """
    section = _section(page, "scope")
    rule = section.index('<p class="scope-rule">')
    key = section.index('<div class="grade-key">')
    assert rule < key, "the grade key renders before the rule it says is above it"
