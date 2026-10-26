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
  "listed here" or "the panel record above". Its third paragraph states in
  full what the card's two warning chips abbreviate to two words each.
* `#rare-variants` is both: three paragraphs moved off every gene page (the
  folded reading notes and the composite note's second sentence), and two are
  the full forms of notices a gene page keeps in one line. Two more state in
  general what each forest plot's folded key says of its own marks (D57), and
  what a block with no effect measure shows instead of a plot; one more, what
  the `all cases` case group is, which a forest key says in one clause
  (2026-10-26).
* `#cohorts` is curated data under one written paragraph: every record of
  `curation/cohorts.yaml`, the same ones `cohorts.json` publishes. A gene page
  keeps the notes of its own cohorts, folded, and links here (P4).

`tests/unit/test_build_methods.py` pins the scope rule and the key's opening
paragraph against literals, sliced to the section; the glosses by a phrase
unique to each, as one label-chip-gloss triple per rung; `#validity`,
`#rare-variants` and, on a fixture registry, `#cohorts` whole, as literals; the
cohort intro's figures against the mirror and the committed registry; and the
phrases the move retired as absent from the whole page, built with no cohorts
and with the committed registry.
`tests/test_built_site_is_consumable.py` checks `#cohorts` against the
`cohorts.json` of a real build.

Written through `Emitter.write_text`, for the reason `pages.py` gives: `Json`'s
union includes `str`, so `write_json` would type-check and publish the page as a
quoted JSON string. Nothing here reads a clock or iterates a set, so two builds
of one commit write the same bytes.
"""

from __future__ import annotations

import html
from collections.abc import Mapping
from typing import Final

from chd_atlas.build.emit import Emitter
from chd_atlas.build.paths import COHORTS, METHODS_PAGE
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
#
# Its first part was rewritten on 2026-10-26, after review measured it false
# read as a rule (each measured that day on the committed mirrors and the
# pinned HPO release):
#
# * It named a *property of a disease* where the code reads a *list*. A disease
#   is in scope only once it is in `curation/chd_scope.yaml`, on either path:
#   SCP007 reports 239 unnominated terms that would qualify on Orphanet's
#   frequencies (FLNA's periventricular nodular heterotopia among them), and
#   Orphanet's own GenCC submissions name Ebstein anomaly (MONDO:0009144),
#   which is not a scope term. So the sentence says what a term needs to be
#   listed, and that a qualifying disease is not in scope until it is.
# * "independently assert it" was weaker than the gate, which counts a GenCC
#   submitter only at `Limited` or better (`validity.agreeing_records`).
# * "Structural heart malformation" is now defined, from
#   `orphanet_scope.structural_phenotypes` -- and two clauses are narrower
#   than the text this was drafted from. The excluded branches are named as
#   *excluded*, not as all the branches that are not malformations: the kept
#   ones still hold aortic valve calcification, mitral valve prolapse and left
#   ventricular hypertrophy, and the sentence says so. And a subtype qualifies
#   on its syndrome's frequency (`Warrant.via == "subtype"`), which the
#   parenthesis states. `test_the_structural_set_the_build_applies_is_the_one_
#   methods_html_states` (`test_orphanet_scope.py`) ties every term named here
#   to the set the build computes, so changing the rule where this names it
#   fails there until this is reworded.
SCOPE_RULE: Final = (
    '<p class="scope-rule">A gene is published here when an external authority already '
    "treats it as a congenital heart disease gene: either a <strong>ClinGen expert "
    "panel</strong> classifies it <strong>Limited or better</strong> for a disease in scope "
    "here, or <strong>two or more Gene Curation Coalition submitters</strong> classify it "
    "Limited or better for one and no ClinGen panel disputes it. A disease is in scope when "
    "it is one of the terms this atlas lists, and a term is listed only if an external "
    "authority treats it as congenital heart disease, or if a curator nominates a syndrome "
    "and <strong>Orphanet records a structural heart malformation in at least 30% of its "
    "cases</strong> (for a subtype, of the syndrome it belongs to); a disease that would "
    "qualify is not in scope until it is listed. Here a structural heart malformation is a "
    "finding under the Human Phenotype Ontology's <em>Abnormal heart morphology</em> "
    "(HP:0001627), less branches of it excluded as not malformations (cardiomyopathy, "
    "tumours, the pericardium and the like), plus coarctation of the aorta and a patent "
    "ductus arteriosus. The branches kept still hold some findings that are not "
    "malformations, such as a calcified or prolapsing valve; lesions the ontology files "
    "elsewhere, among them pulmonic stenosis, truncus arteriosus and an interrupted aortic "
    "arch, are not counted. "
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
#
# The third paragraph glosses the card's two flags (2026-10-26), which no page
# explained: `methods.html` mentioned "sources disagree" zero times, and on all
# six pages carrying it ClinGen and GenCC agree on every gene-disease pair they
# both classify. Each clause is read off the code that sets the flag:
# `vocab.has_conflicting_evidence` over every in-scope record with a mapped
# classification (`derive.gene_facts`), and `validity._has_source_discordance`,
# whose sources are `ValiditySource.CLINGEN` and `ValiditySource.GENCC` -- so a
# ClinGen copy submitted to GenCC is GenCC's there. `NO_KNOWN_ASSOCIATION` is
# subtracted from the supporting side and is not in `CONTESTED`; `Supportive`
# maps to no classification and is read by neither. The flag is per gene by the
# owner's decision (2026-10-25), so the paragraph says what that permits rather
# than narrowing it. "Narrower" is exact: a gene with `has_source_discordance`
# has a contesting and a supporting record, which is `has_conflicting_evidence`.
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
    "<p>A gene page's card may carry two flags beside its grade. <em>conflicting "
    "evidence</em> means that among the gene's records for diseases in scope, at least one "
    "is <em>Disputed</em> or <em>Refuted</em> and at least one classifies the gene Limited "
    "or better, whoever made each and whatever disease each is for; the gene index marks "
    "the same genes beside their confidence. <em>sources disagree</em> is narrower: a "
    "record that disputes or refutes the gene and one that classifies it Limited or better "
    "come from different sources, one a ClinGen expert panel's own record and the other a "
    "submission to GenCC "
    "(ClinGen's own copies count as GenCC's here). It too is a flag on the gene, not on a "
    "disease: the two records may be for different diseases, and the two sources may agree "
    "on every disease they both classify. Neither flag counts a <em>No Known Disease "
    "Relationship</em> record on either side, nor GenCC's <em>Supportive</em>, which "
    "grades no evidence.</p>"
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
    # The case-group union, which nothing here stated until 2026-10-26: the
    # paragraph above explained `damaging` while a table's `all cases` row sat
    # beside the rows it contains. Measured that day on the committed mirror:
    # all three studies publish an `all` stratum and only PMID:42230622 also
    # splits its cases, into syndromic (1,471) and non-syndromic (2,405); on all
    # 470 cells where it publishes all three, the `all` row's cases and carriers
    # are the sum of the other two's and its controls are theirs. Stated in
    # general, so it names no study; `test_an_all_cases_row_is_the_union_of_
    # the_case_groups_its_study_splits` holds it to the mirror.
    #
    # The last clause said `all cases` was a study's "whole case set" until
    # review on 2026-10-26: it is not a denominator. PMID:34324492's rows count
    # 7,107 cases of the 7,958 its own note gives the study, and
    # PMID:40127276's count 3,887 trios or an allele total against 11,555
    # probands. What holds of every study is the narrower claim: not a subgroup.
    "<p>Where a study splits its cases into <strong>syndromic</strong> and "
    "<strong>non-syndromic</strong> groups, its <strong>all cases</strong> group is "
    "likewise the union of the two, and an <strong>all cases</strong> row is not a result "
    "independent of theirs. Where a study does not split its cases, <strong>all "
    "cases</strong> is its only case group: every case the row counts, not a subgroup.</p>"
    # The forest key in general; each panel's own key, conditional clause by
    # clause, is folded under "How to read this chart" (`pages._forest_caption`,
    # D57). Each clause measured against a real build of 9f24190 before it was
    # stated in general: all 29 right arrows are odds ratios with 0 control
    # carriers and a published lower bound; all 257 left arrows have 0 case
    # carriers; every plotted study corrects every row or none. "contains other
    # rows in the same panel" was false for 10 union-tagged rows, on 10 panels,
    # that contain exactly one other drawn row, so it says "at least one", as
    # each panel's key does.
    #
    # Two clauses reworded after review (2026-10-26): "these have none" had no
    # antecedent nearer than "a forest plot", so it names this atlas's; and the
    # muted row was singular, while every panel that draws one draws three
    # (measured on a real build: 113 panels, one synonymous row per case group
    # of PMID:42230622), so it is "one of a study's synonymous rows", true for
    # one or many.
    #
    # One sentence added 2026-10-26, when a marker's fill began reading
    # `concordance.survives_own_correction` for every row but a synonymous one
    # (`pages._survived_correction`): a depleted row that cleared its
    # correction is drawn hollow, and "did not" survive would be false of it.
    # 0 published rows are in that state (measured that day). "Non-synonymous"
    # since review the same day: a synonymous row that survived its correction is
    # drawn filled and muted (`pages._mark_class`), never hollow -- latent, as no
    # mirrored synonymous row carries a corrected p.
    "<p>A forest plot shows each row's effect and interval on a logarithmic axis. A "
    "<strong>filled</strong> marker is a row that survived its own study's correction for "
    "multiple testing; a <strong>hollow</strong> marker is a row that did not, or whose "
    "study published no correction at all &mdash; the atlas computes none of its own. A "
    "non-synonymous row that survived its correction without pointing to enrichment in "
    "cases is drawn hollow too: it is not evidence for the gene. A forest plot usually "
    "ends in a summary diamond; this atlas's have none, because no pooled estimate is "
    "computed anywhere on this site. A <strong>hatched</strong> bar "
    "means no interval was published for that row, not a wide one. An arrow to the "
    "<strong>right</strong> is an effect unbounded above &mdash; no control carried the "
    "variant &mdash; so the study's lower bound is the whole finding; an arrow to the "
    "<strong>left</strong> is an effect of zero &mdash; no case carried the variant "
    "&mdash; which a logarithmic axis has no position for. A row tagged "
    "<strong>union</strong> contains at least one other row in the same panel and is not "
    "an independent finding. A <strong>muted</strong> row is one of a study's synonymous "
    "rows, drawn apart from the results it calibrates.</p>"
    # The general form of `pages._NO_EFFECT_MEASURE`, which a block shows in
    # place of a plot (77 pages, every one in PMID:34324492's block).
    "<p>A study that publishes <strong>no effect measure</strong> for a gene &mdash; no "
    "odds ratio, enrichment or rate ratio &mdash; has nothing to place on an effect axis, "
    "so its block shows its counts and p-values in a table instead of a plot.</p>"
)

# Every collection the rare-variant rows cite, with its curated description
# (P4, 2026-10-26). A gene page folds the notes of only the collections its own
# rows cite (`pages._cohort_notes`, whose fold ends on a link here) and keeps
# them: CLAUDE.md section 9 forbids narrowing where a cohort caveat reaches, so
# this is a second place, not a move. It is the whole registry -- the records
# `cohorts.json` publishes -- and every registered collection is cited by a
# mirror row, so "every sample collection a rare-variant row names" describes
# the list both ways; `test_the_cohort_list_says_what_its_rows_name_and_what_
# they_do_not` holds that against the mirror.
#
# It read "every sample collection a rare-variant study in this atlas draws on,
# with what qualifies every number drawn from it" until review (2026-10-26)
# measured both halves false. PMID:34324492's rows name 4 collections
# (`decipher`, `isca`, `signature`, `taa_cases`), and the study's S2 Table lists
# 189 case sources holding 7,958 cases, of which those 4 hold 4,242 -- counted
# from the supplement's `s007.xlsx` (sheet "S2 Table": 189 source rows whose
# individuals sum to 7,958, matching the sheet's own "N" row; Signature 1,615,
# Decipher_database 1,221, PMID 21092924 777, ISCA 629). A reader told the list
# is complete concludes that two studies naming no collection in common share
# no people, which is the inference the next sentence must not license; the
# gene page's pooling sentence was already hedged ("as this atlas has recorded
# them"). And some descriptions qualify nothing (`isca`'s is a case count), so
# the list carries "the caveats its curator recorded", whatever they are. The
# sources the rows do not name are not named here either: an unregistered
# source has no record to resolve against. The 4 is pinned against the mirror,
# 4,242 against the four descriptions' own counts and 7,958 against
# `taa_cases`'s; 189 is the supplement's alone and nothing here can check it.
_COHORTS_INTRO: Final = (
    '<h2 id="cohorts">Cohorts</h2>'
    "<p>Every sample collection a rare-variant row in this atlas names, with the caveats "
    "its curator recorded. A row may not name every collection its study drew on: the "
    "rows for PMID:34324492 name 4 of the 189 case sources its S2 Table lists, and those "
    "4 hold 4,242 of its 7,958 cases. Two studies that name the same collection describe "
    "partly the same people, and two that name none in common may still share people "
    "this registry does not record. The same records are published as "
    f'<a href="{COHORTS}"><code>{COHORTS}</code></a>.</p>'
)


def _cohorts(cohorts: Mapping[str, Cohort]) -> str:
    """The registry as a definition list: name and id, then the description.

    Sorted by id here, because `burden.cohort_registry` returns the YAML's order
    and an array derived from a mapping is the caller's to sort (CLAUDE.md
    section 9) -- the same order `cohorts.json` publishes. Every curated field
    is escaped: names, descriptions and urls are free text. A url is the
    collection's public page, linked inside its own description; a cohort with
    none (4 of 13) gets no link rather than an empty one. The link reads
    "website", not "source" (2026-10-26): it leads to the collection's home
    page, and "source" beside a description reads as the citation for it.

    `cohort-list` is layout only. The base `dl` is a two-column grid whose label
    column grows to its widest label, which suits the grade key's chips and
    would give a cohort's full name most of the line.
    """
    items = "".join(
        f"<dt>{html.escape(cohort.name)} (<code>{html.escape(cohort.id)}</code>)</dt>"
        f"<dd>{html.escape(cohort.description)}"
        + (f' <a href="{html.escape(cohort.url)}">website</a>' if cohort.url else "")
        + "</dd>"
        for cohort in sorted(cohorts.values(), key=lambda cohort: cohort.id)
    )
    return f'{_COHORTS_INTRO}<dl class="cohort-list">{items}</dl>'


# The page's own account of what it holds (CLAUDE.md section 4.29); the cohort
# list joined it on 2026-10-26. It names only sections the page renders: the
# expression section is not on this page yet, and naming it here would promise
# a reader something the page does not hold. "the rare-variant studies draw on"
# became "the rare-variant results name" the same day, for the reason
# `_COHORTS_INTRO` records: the list is what the rows name, not all a study
# drew on.
_INTRO: Final = (
    "<h1>How to read this atlas</h1>"
    "<p>What every gene page assumes, written once: which genes are here and what each "
    "classification means, how the validity records are listed, how to read the "
    "rare-variant results, and every sample collection the "
    "rare-variant results name.</p>"
)


def build_methods(emitter: Emitter, *, cohorts: Mapping[str, Cohort]) -> None:
    """Emit `methods.html`, the site's reading notes (D50).

    `cohorts` is the registry `runner.build_site` already hands the gene pages
    (`burden.cohort_registry`); `#cohorts` lists all of it.
    """
    sections = (
        '<h2 id="scope">Which genes are here, and what a classification means</h2>'
        + SCOPE_RULE
        + grade_key()
        + _VALIDITY
        + _RARE_VARIANTS
        + _cohorts(cohorts)
    )
    emitter.write_text(
        METHODS_PAGE,
        document(title="How to read this atlas — CHD Atlas", root="", body=_INTRO + sections),
    )
