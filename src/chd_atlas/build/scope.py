# src/chd_atlas/build/scope.py
"""Publish `scope.json`, and the Orphanet warrants the gene pages explain.

A gene bundle names each validity record's disease by bare MONDO id, and the
reason that disease counts as congenital heart disease lived only in
`curation/chd_scope.yaml`. That was tolerable while every reason was an
authority's own curation. Since 2026-10-25 (issue #8 part 1) a syndrome can be in
scope because Orphanet *records* structural heart malformations in at least 30%
of its cases -- a frequency, not a curation -- and a reader who meets PTPN11 in a
CHD atlas needs that sentence to read the page correctly. So the registry is
published whole, the way `cohorts.json` publishes the cohort registry.

The warrant's phenotypes are re-derived from the Orphanet mirrors on every build
rather than stored in YAML, so `scope.json` and the gene-page note cannot drift
from what Orphanet records. They are derived by `orphanet_scope.warrant`, the
same function SCP006 refused on, so the build and the gate are one reading of
one rule.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping

from chd_atlas.build.emit import Emitter, Json
from chd_atlas.build.paths import SCOPE
from chd_atlas.models.scope import ScopeAuthority, ScopeEntry
from chd_atlas.orphanet_scope import OrphanetTables, Warrant, warrant


def scope_warrants(
    entries: Iterable[ScopeEntry],
    tables: OrphanetTables | None,
    structural: frozenset[str] | None,
) -> dict[str, Warrant]:
    """Each Orphanet-scoped term's warrant, keyed by MONDO id.

    **Raises** rather than skipping an unwarranted term. SCP006 is an error, so
    `validate_repository` has already refused a corpus with one, and reaching it
    here means the gate was bypassed. Publishing the term anyway would print, on
    every gene page it admits, a scope sentence the data does not support. A
    bypassed gate must fail rather than publish (CLAUDE.md section 9), and it is
    `raise`, never `assert`: `-O` strips `assert`.
    """
    nominated = [
        entry
        for entry in entries
        if entry.admitted_by is ScopeAuthority.ORPHANET_PHENOTYPE_FREQUENCY
    ]
    if not nominated:
        return {}
    if tables is None or structural is None:
        raise ValueError(
            f"{len(nominated)} Orphanet-scoped term(s) but no Orphanet mirrors or HPO release "
            f"to derive their warrants from; SCP006 should have refused this corpus"
        )
    warrants: dict[str, Warrant] = {}
    for entry in nominated:
        result = warrant(entry.id, entry.orphanet or "", tables, structural)
        if isinstance(result, str):
            raise ValueError(f"{entry.id}: {result}; SCP006 should have refused this corpus")
        warrants[str(entry.id)] = result
    return warrants


def _orphanet_payload(found: Warrant | None) -> Json:
    if found is None:
        return None
    return {
        "code": found.orpha,
        "name": found.name,
        "via": found.via,
        "omim": found.omim,
        "phenotypes": [
            {"hpo": annotation.hpo, "label": annotation.label, "frequency": annotation.frequency}
            for annotation in found.phenotypes
        ],
    }


def build_scope(
    entries: Iterable[ScopeEntry], warrants: Mapping[str, Warrant], emitter: Emitter
) -> None:
    """Write `scope.json`: every scope term, sorted by id, with its authority and warrant.

    Sorted by id because `encode_json`'s `sort_keys` orders dict keys and never
    list elements. `orphanet` is `null`, not absent, on a term an authority
    curated, so every object has one shape; `syndrome` is always a boolean, and
    `review_status` and `review_note` are `null` together on a term with no open
    question recorded.
    """
    emitter.write_json(
        SCOPE,
        {
            "terms": [
                {
                    "id": str(entry.id),
                    "label": entry.label,
                    "admitted_by": entry.admitted_by.value,
                    "attributed_to": entry.attributed_to,
                    "orphanet": _orphanet_payload(warrants.get(str(entry.id))),
                    # A curated mark, `false` unless set (R5): the term names a
                    # syndrome, so a grade for it is a grade for the syndrome.
                    "syndrome": entry.syndrome,
                    # An open question about the term (owner's standing rule,
                    # 2026-10-27), `null` and `null` where none is recorded.
                    # Published here only: no page renders it.
                    "review_status": (
                        None if entry.review_status is None else entry.review_status.value
                    ),
                    "review_note": entry.review_note,
                }
                for entry in sorted(entries, key=lambda entry: entry.id)
            ]
        },
    )
