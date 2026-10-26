# src/chd_atlas/build/methods.py
"""`methods.html`, "How to read this atlas": site-wide explanation, written once.

D50. Every gene page repeated the same explanation -- the rule that admits a
gene, the grade key and, in later sections, the reading notes for the validity,
rare-variant and expression evidence -- so a reader met the same paragraphs on
every page and the part particular to the gene was buried under them. This page
holds that explanation once; each gene page keeps a short form and a link to the
section here that states it in full.

**Every sentence below the introduction was moved off a gene page or the gene
index, and each was re-read in this position.** A sentence's truth conditions
travel with its position (CLAUDE.md section 4.27): "the table below" was true on
a gene page, false on the gene index and false here, and "this page" here names
a page with no gene on it. Those that were false here were reworded: the scope
rule's last clause, and the grade key's dispute and completeness sentences.
The rest are unchanged because they are true anywhere -- the key's opening
sentences, and the rung glosses (`render.GRADE_GLOSS`), which stay verbatim on
purpose because a gene page's card glosses its headline grade from the same
constant, and one rung must not be glossed two ways.

`#validity` (2026-10-26) is the one section whose gene-page sentences did not
leave: the ClinGen-copy notice stays on the page as it is (P5), and the table's
intro keeps a short form of its scope. The section states in full what those
abbreviate, written for a page with no table on it -- so "a gene page lists",
never "listed here" or "the panel record above".

`tests/unit/test_build_methods.py` pins the scope rule and the key's opening
paragraph against literals, sliced to the section; the glosses by a phrase
unique to each, as one label-chip-gloss triple per rung; and the phrases the
move retired as absent from the whole page.

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
_VALIDITY: Final = (
    '<h2 id="validity">Mirrored validity records</h2>'
    "<p>Every classification in this atlas is an upstream panel's or submitter's, "
    "republished with its provenance intact. A gene page lists them in two groups: the "
    "ClinGen expert panel's own records first, then other curators' submissions to the "
    "Gene Curation Coalition (GenCC). ClinGen also submits its own classifications to "
    "GenCC; a gene page lists such a copy only where its mapped grade, date, inheritance "
    "or disease differs from the panel's own record, or where it carries no date or has "
    "no panel record to match. Each gene's JSON keeps every copy, listed or not.</p>"
    "<p>Only records for diseases in this atlas's remit are listed. A panel may have "
    "classified the same gene for a disease outside it, and that record is not shown, so "
    "a gene page may not be every classification the gene has.</p>"
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
    )
    emitter.write_text(
        METHODS_PAGE,
        document(title="How to read this atlas — CHD Atlas", root="", body=_INTRO + sections),
    )
