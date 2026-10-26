# src/chd_atlas/build/methods.py
"""`methods.html`, "How to read this atlas": site-wide explanation, written once.

D50. Every gene page repeated the same explanation -- the rule that admits a
gene, the grade key and, in later sections, the reading notes for the validity,
rare-variant and expression evidence -- so a reader met the same paragraphs on
every page and the part particular to the gene was buried under them. This page
holds that explanation once; each gene page keeps a short form and a link to the
section here that states it in full.

**Every sentence below the introduction either was moved off a gene page or the
gene index, or was written here to state in full what a gene page now says in
one line -- and each was read in this position.** A sentence's truth conditions
travel with its position (CLAUDE.md section 4.27): "the table below" was true on
a gene page, false on the gene index and false here, and "this page" here names
a page with no gene on it.

* `#scope` was moved: the scope rule and the grade key. Those that were false
  here were reworded -- the scope rule's last clause, and the grade key's
  dispute and completeness sentences. The rest are unchanged because they are
  true anywhere: the key's opening sentences, and the rung glosses
  (`render.GRADE_GLOSS`), which stay verbatim on purpose because a gene page's
  card glosses its headline grade from the same constant, and one rung must not
  be glossed two ways.
* `#validity` was written new (2026-10-26). Its gene-page sentences did not
  leave: the ClinGen-copy notice stays on the page (P5), and the table's intro
  keeps a short form of its scope. The section states in full what those
  abbreviate, for a page with no table on it -- so "a gene page lists", never
  "listed here" or "the panel record above".
* `#rare-variants` is both: three paragraphs moved off every gene page (the
  folded reading notes and the composite note's second sentence), and two are
  the full forms of notices a gene page keeps in one line. Two more state in
  general what each forest plot's folded key says of its own marks (D57), and
  what a block with no effect measure shows instead of a plot.

`tests/unit/test_build_methods.py` pins the scope rule and the key's opening
paragraph against literals, sliced to the section; the glosses by a phrase
unique to each, as one label-chip-gloss triple per rung; `#validity` and
`#rare-variants` whole, as literals; and the phrases the move retired as absent
from the whole page.

Written through `Emitter.write_text`, for the reason `pages.py` gives: `Json`'s
union includes `str`, so `write_json` would type-check and publish the page as a
quoted JSON string. Nothing here reads a clock or iterates a set, so two builds
of one commit write the same bytes.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Final

from chd_atlas.build.emit import Emitter
from chd_atlas.build.paths import METHODS_PAGE
from chd_atlas.build.render import document, grade_key
from chd_atlas.models.cohort import Cohort

# What the published set *is*, stated in full here and in one line on both page
# kinds that show it (`pages._SCOPE_SHORT`, which links this section).
#
# No page on this site said it. `docs/data-api.md` states the rule and a reader
# of the HTML never sees it, so a browse page of 23 rows all reading
# `definitive`, under a title reading "CHD Atlas" and a tagline about congenital
# heart disease, said "these 23 genes are definitive for CHD". Measured
# 2026-08-04 against the committed mirrors, that was false for most of them: 10
# of the 23 were qualified by a panel that is not the Congenital Heart Disease
# GCEP -- SCID-CID for KMT2D and KDM6A, Hearing Loss for CHD7, Kidney Cystic and
# Ciliopathy Disorders for NOTCH2 and EVC, Intellectual Disability and Autism
# for ZEB2, Skeletal Disorders for EVC2, Syndromic Disorders for TBX5 and
# TFAP2B, Prenatal for FLT4 -- and 10 were definitive for a disease whose label
# names no cardiac feature at all (Kabuki syndrome 1 and 2, Mowat-Wilson,
# CHARGE, Alagille, Ellis-van Creveld twice, Holt-Oram, TARP, and NR2F2-related
# multiple congenital anomalies).
#
# The 2026-08-06 widening makes the sentence carry more, not less. Measured on
# the widened corpus: 92 genes, of which 76 are admitted by a ClinGen record at
# or above `Limited` and 16 by two or more GenCC submitters with no ClinGen
# record at all. A reader meeting `limited` needs to know both that a chartered
# panel wrote it and that this atlas did not.
#
# Moved from `pages._SCOPE_RULE` (2026-10-26, D50), where it was one constant
# for both page kinds. Its last clause said the disease "is shown beside the
# gene", true on a gene page and on the gene index and unanchored here, so it
# now names both places -- and says "for a gene a ClinGen panel graded", since
# the 21 genes admitted on submitter agreement have no panel-named disease to
# show (their card names the submitters instead).
SCOPE_RULE: Final = (
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

# How a gene page's validity table is listed (2026-10-26, D50). The page's own
# intro states the short form -- only in-scope records are listed -- and links
# here. The copy sentence describes `pages._repeats_the_panel` and is copied
# from `docs/data-api.md`'s account of it, read against the function that day:
# a copy is hidden only on a match on disease, mapped grade, inheritance and
# calendar date, and a copy with no date matches nothing.
#
# Two clauses are narrower than the text this was drafted from, each measured
# 2026-10-26 on the committed mirrors and a real build:
#
# * "Each gene's JSON keeps every record" read, beside the next paragraph's
#   "that record is not shown", as though the JSON held the out-of-scope
#   records. It holds none -- 591 bundle records, 0 for a disease outside
#   scope -- and all 121 copies, so it says "every copy".
# * "a gene page is not every classification the gene has" was false for 43 of
#   the 119 published genes, which have no ClinGen or GenCC record for any
#   disease outside scope. "may not be", as the gene page says.
#
# Reworded after review the same day. It opened "Every classification in this
# atlas is an upstream panel's or submitter's, republished with its provenance
# intact", which the grade key's first sentence already says two paragraphs up;
# it is said once, there. "the ClinGen expert panel's own records" became
# "ClinGen expert panels'", and "A panel may have classified" became "A panel or
# submitter": the scope filter drops GenCC rows too, and most out-of-scope
# records are GenCC's (76 published genes have one, 27 an out-of-scope ClinGen
# record).
_VALIDITY: Final = (
    '<h2 id="validity">Mirrored validity records</h2>'
    "<p>A gene page lists its classifications in two groups: ClinGen expert panels' own "
    "records first, then other curators' submissions to the Gene Curation Coalition "
    "(GenCC). ClinGen also submits its own classifications to GenCC; a gene page lists "
    "such a copy only where its mapped grade, date, inheritance or disease differs from "
    "the panel's own record, or where it carries no date or has no panel record to "
    "match. Each gene's JSON keeps every copy, listed or not.</p>"
    "<p>Only records for diseases in this atlas's remit are listed. A panel or submitter "
    "may have classified the same gene for a disease outside it, and that record is not "
    "shown, so a gene page may not be every classification the gene has.</p>"
)

# How to read a gene page's rare-variant section (2026-10-26, D50). The page
# keeps a one-line short form of the first two paragraphs (`pages._MIRRORED_NOTICE`
# and `pages._POOLING_NOTICE`) and links here; the last three were the page's
# folded reading notes and its composite note's second sentence, and left it.
#
# Narrower than the text this was drafted from, each measured 2026-10-26 on a
# real build of 4233eca:
#
# * "Every count and statistic in a gene page's rare-variant section" -> "in a
#   gene page's study tables". The evidence matrix in that section shows tallies
#   the atlas derives -- "n = N" in a cell on 86 pages, "N rows" on 48 -- and
#   each PMID:42230622 block carries the atlas's own Bonferroni arithmetic
#   (`0.05 / 138,609`); neither is an upstream study's number.
# * "republished exactly as published" -> "as published". The tables round to
#   three significant figures: 376 of the 1,210 published raw p-values render
#   rounded (PMID:40127276's 0.00583259259259259 reads 0.00583). Each gene's
#   JSON carries them in full.
# * The units sentence said `cases` and `controls` "each show observed / total".
#   On a row tested against a mutation model the `controls` cell is the model's
#   expected count (`pages._controls`, which keys on the comparator): 108 cells
#   on 36 pages read "0.0358 expected" and the like. It says so.
#
# Reworded after review the same day: "children" -> "people", because
# `shared_cohorts` counts shared control collections too and a control cohort
# may be adults (UK Biobank's note: "~500,000 UK adults recruited at ages
# 40-69"); the overlap clause says what the page can know -- collections shared
# *as this atlas has recorded them* -- since `shared_cohorts` sees only shared
# registry ids, and two studies may share individuals no entry names (review
# raised DDD and DECIPHER; the registry records no link between them); the
# synonymous paragraph
# lost "the same uncorrected scale" (true of the one study that publishes a
# synonymous row, stated in general) and "a scan of that size" (no antecedent
# here).
_RARE_VARIANTS: Final = (
    '<h2 id="rare-variants">Rare-variant association</h2>'
    "<p>Every count and statistic in a gene page's study tables is an upstream study's, "
    "republished as published. The atlas computes none of them and adds no assessment "
    "of them; where a study is by an author of this atlas, its block says so.</p>"
    # was pages._POOLING_NOTICE in full
    "<p>The atlas computes <strong>no pooled statistic across studies</strong>. Where two "
    "studies draw on the same sample collection they describe partly the same people, "
    "so combining their results would count those people twice; a gene page names any "
    "sample collection its studies share, as this atlas has recorded them.</p>"
    # was pages._BURDEN_PREAMBLE; "In each study's table," anchors what it describes
    "<p>In each study's table, <code>cases</code> and <code>controls</code> each show "
    "<em>observed / total</em> (for a row tested against a mutation model, "
    "<code>controls</code> is instead the count that model expected), and <strong>every "
    "cell names what it counted</strong> &mdash; carriers, alleles, or de novo mutations. "
    "The three are not interchangeable: a study counting alleles counts a person with two "
    "qualifying variants twice, and its denominator is roughly twice the number of people "
    "sequenced. A consequence class or case group with no row is one the study "
    "<strong>did not report</strong>, which is not a null result and does not always "
    "mean no such variant was seen &mdash; studies differ in which classes their "
    "gene-level tables cover.</p>"
    # was pages._SYNONYMOUS_NOTICE; "above it" -> "beside it"
    "<p>The <strong>synonymous</strong> row is a study's own negative control &mdash; "
    "synonymous variants should show no enrichment. Read it on the same scale as the rows "
    "beside it: a <em>strongly</em> enriched synonymous row is a warning about that "
    "gene's comparison, while a nominally significant one is what a scan over thousands "
    "of genes produces by chance. <strong>Only a study that publishes one has "
    "it</strong>; a table with no synonymous row offers no such check.</p>"
    # was the second sentence of pages._composite_note
    "<p>A <strong>damaging (LOF + missense)</strong> row is the <strong>union</strong> of "
    "the study's loss-of-function and damaging-missense rows, not a third independent "
    "result &mdash; the same variants are counted in both. It is shown because a study "
    "that defines its findings by the composite publishes a p-value for it that is not "
    "derivable from the other two.</p>"
    # The forest key in general; each panel's own key, conditional clause by
    # clause, is folded under "How to read this chart" (`pages._forest_caption`,
    # D57). Each clause measured against a real build of 9f24190 before it was
    # stated in general: all 29 right arrows are odds ratios with 0 control
    # carriers and a published lower bound; all 257 left arrows have 0 case
    # carriers; every plotted study corrects every row or none. "contains other
    # rows in the same panel" was false for 10 union-tagged rows, on 10 panels,
    # that contain exactly one other drawn row, so it says "at least one", as
    # each panel's key does.
    "<p>A forest plot shows each row's effect and interval on a logarithmic axis. A "
    "<strong>filled</strong> marker is a row that survived its own study's correction for "
    "multiple testing; a <strong>hollow</strong> marker is a row that did not, or whose "
    "study published no correction at all &mdash; the atlas computes none of its own. A "
    "forest plot usually ends in a summary diamond; these have none, because no pooled "
    "estimate is computed anywhere on this site. A <strong>hatched</strong> bar means no "
    "interval was published for that row, not a wide one. An arrow to the "
    "<strong>right</strong> is an effect unbounded above &mdash; no control carried the "
    "variant &mdash; so the study's lower bound is the whole finding; an arrow to the "
    "<strong>left</strong> is an effect of zero &mdash; no case carried the variant "
    "&mdash; which a logarithmic axis has no position for. A row tagged "
    "<strong>union</strong> contains at least one other row in the same panel and is not "
    "an independent finding. A <strong>muted</strong> row is a study's synonymous row, "
    "drawn apart from the results it calibrates.</p>"
    # The general form of `pages._NO_EFFECT_MEASURE`, which a block shows in
    # place of a plot (77 pages, every one in PMID:34324492's block).
    "<p>A study that publishes <strong>no effect measure</strong> for a gene &mdash; no "
    "odds ratio, enrichment or rate ratio &mdash; has nothing to place on an effect axis, "
    "so its block shows its counts and p-values in a table instead of a plot.</p>"
)

_INTRO: Final = (
    "<h1>How to read this atlas</h1>"
    "<p>What every gene page assumes, written once: which genes are here and what each "
    "classification means, how the validity records are listed, and how to read the "
    "rare-variant results and the expression figures.</p>"
)


def build_methods(emitter: Emitter, *, cohorts: Mapping[str, Cohort]) -> None:
    """Emit `methods.html`, the site's reading notes (D50).

    `cohorts` is the registry `runner.build_site` already hands the gene pages
    (`burden.cohort_registry`). No section reads it yet; it is accepted now so
    the runner's call does not change when the cohort section arrives.
    """
    sections = (
        '<h2 id="scope">Which genes are here, and what a classification means</h2>'
        + SCOPE_RULE
        + grade_key()
        + _VALIDITY
        + _RARE_VARIANTS
    )
    emitter.write_text(
        METHODS_PAGE,
        document(title="How to read this atlas — CHD Atlas", root="", body=_INTRO + sections),
    )
