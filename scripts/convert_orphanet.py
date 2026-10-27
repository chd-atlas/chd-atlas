# scripts/convert_orphanet.py
"""Convert Orphanet products 4 and 1 and MONDO's SSSOM into the three scope mirrors.

Run locally by a curator, never by CI or the build -- design decision D6.

    .venv/bin/python scripts/convert_orphanet.py \\
        local/data/orphanet/en_product4.xml local/data/orphanet/en_product1.xml \\
        local/data/mondo/mondo.sssom.tsv mirrors/

Writes `mirrors/orphanet_phenotypes.tsv`, `mirrors/orphanet_omim.tsv` and
`mirrors/mondo_xrefs.tsv`, which `chd_atlas.orphanet_scope` reads to decide
whether a nominated syndrome is in scope (issue #8 part 1).

LICENCES
--------
Both Orphanet files declare `CC-BY-4.0` in their own XML header (`<Licence>`),
and Orphadata's legal notice says the same of its Orphadata Science files,
asking to be cited as "Orphadata Science: Free access data from Orphanet.
(c) INSERM 1999. Available on http://sciences.orphadata.com/. Data version
[date]". Its AI-use clause covers LLM/AI tools offered to third parties; this
atlas is a static data site and publishes no such tool. MONDO is CC-BY-4.0.
Each file's sha256 is recorded in `mirrors/sources.yaml`.

THE DATA VERSION IS WRITTEN, NOT PRINTED
----------------------------------------
Both Orphanet files open on `<JDBOR date="..." version="...">`, the data
version Orphadata's citation asks for. Since owner decision 2026-10-26 (R11)
this converter reads that header and writes it into the file's own entry in
`mirrors/sources.yaml` -- `version` (the sha256), `file_version` and
`file_date` -- every time it runs, so the checksum, the version and the date
cannot come from three different downloads. The landing page's Orphadata
citation prints them from there. MONDO's sha256 is still printed only.

WHAT IS KEPT, AND WHY IT IS WIDER THAN THE RULE
-----------------------------------------------
* **Product 4: every annotation under HP:0001626 (cardiovascular).** The rule
  counts only the structural subset of HP:0001627, but keeping the whole
  cardiovascular branch means a change to that subset needs no re-mirror.
  Measured 2026-10-25: 7,902 of 116,664 annotations, over 2,182 disorders. (An
  estimate of 7,895 that ignored alt ids was 7 rows short: some annotations carry
  a retired HPO id, which the closure must include.)
* **Product 1: every Orphanet-to-OMIM reference**, whatever its relation or
  validation status. The rule keeps only validated `E`/`BTNT`; the mirror keeps
  the rest so the refusals are visible. 8,745 references.
* **MONDO: exact matches to Orphanet or OMIM, for terms either validity mirror
  names.** 10,011 rows. Orphanet codes are normalised from MONDO's `Orphanet:648`
  to Orphanet's own `ORPHA:648`, the one spelling every mirror uses.

Product 4 annotates 38 (orpha, hpo) pairs twice, over all its annotations.
Where the band also repeats, the row is an exact duplicate and is collapsed
(2, within the cardiovascular branch); where it differs, both are kept, which is
why the band is in the mirror's key. The report names the collapsed count.

Runnable without the package installed, matching the other converters, which is
why the OBO reader below is a few lines of its own rather than an import of
`chd_atlas.orphanet_scope`.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path
from typing import Final

SOURCE_PRODUCT4: Final = "orphanet_product4"
SOURCE_PRODUCT1: Final = "orphanet_product1"
SOURCE_MONDO: Final = "mondo_sssom"
CARDIOVASCULAR_ROOT: Final = "HP:0001626"
HPO_RELEASE: Final = "ontologies/hp-2026-06-23.obo"

# `TABLE_SCHEMAS`' column orders, written out for the reason `burden_mirror.ORDER`
# is; `test_convert_orphanet.py` asserts each equals its schema.
PHENOTYPE_ORDER: Final = ("orpha", "disorder_name", "hpo", "hpo_label", "frequency", "source")
OMIM_ORDER: Final = ("orpha", "omim", "relation", "validation", "source")
XREF_ORDER: Final = ("mondo", "target", "source")


def _cardiovascular(obo: Path) -> tuple[frozenset[str], frozenset[str]]:
    """HP:0001626's closure with alt ids, and every id the release knows."""
    children: dict[str, set[str]] = defaultdict(set)
    alternates: dict[str, set[str]] = defaultdict(set)
    known: set[str] = set()
    current: str | None = None
    in_term = False
    with obo.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("["):
                in_term = line == "[Term]"
                current = None
            elif not in_term:
                continue
            elif line.startswith("id: "):
                current = line[4:].strip()
                known.add(current)
            elif current is not None and line.startswith("is_a: "):
                children[line[6:].split(" ", 1)[0]].add(current)
            elif current is not None and line.startswith("alt_id: "):
                alternates[current].add(line[8:].strip())
                known.add(line[8:].strip())
    seen, stack = {CARDIOVASCULAR_ROOT}, [CARDIOVASCULAR_ROOT]
    while stack:
        for child in children[stack.pop()]:
            if child not in seen:
                seen.add(child)
                stack.append(child)
    return frozenset(seen | {alt for term in seen for alt in alternates[term]}), frozenset(known)


def _text(element: ET.Element, path: str) -> str:
    return (element.findtext(path) or "").strip()


def _phenotype_rows(element: ET.Element, cardiovascular: frozenset[str]) -> list[dict[str, str]]:
    """One product 4 `HPODisorderSetStatus`, as rows for its cardiovascular annotations."""
    disorder = element.find("Disorder")
    if disorder is None:
        return []
    orpha = f"ORPHA:{_text(disorder, 'OrphaCode')}"
    name = _text(disorder, "Name")
    rows = []
    for association in disorder.iter("HPODisorderAssociation"):
        hpo = _text(association, "HPO/HPOId")
        if hpo not in cardiovascular:
            continue
        rows.append(
            {
                "orpha": orpha,
                "disorder_name": name,
                "hpo": hpo,
                "hpo_label": _text(association, "HPO/HPOTerm"),
                # Verbatim, never mapped to a number: the band is Orphanet's claim.
                "frequency": _text(association, "HPOFrequency/Name"),
                "source": SOURCE_PRODUCT4,
            }
        )
    return rows


def _omim_rows(disorder: ET.Element) -> list[dict[str, str]]:
    """One product 1 `Disorder`, as rows for its OMIM references."""
    orpha = f"ORPHA:{_text(disorder, 'OrphaCode')}"
    rows = []
    for reference in disorder.iter("ExternalReference"):
        if _text(reference, "Source") != "OMIM":
            continue
        rows.append(
            {
                "orpha": orpha,
                "omim": f"OMIM:{_text(reference, 'Reference')}",
                # "BTNT (ORPHAcode is broader than ...)" -> "BTNT"
                "relation": _text(reference, "DisorderMappingRelation/Name").split(" ", 1)[0],
                "validation": _text(reference, "DisorderMappingValidationStatus/Name"),
                "source": SOURCE_PRODUCT1,
            }
        )
    return rows


def _xref_rows(triples: Iterable[tuple[str, str, str]], terms: set[str]) -> list[dict[str, str]]:
    """MONDO exact matches to Orphanet or OMIM for `terms`, normalised and sorted."""
    kept: set[tuple[str, str]] = set()
    for subject, predicate, target in triples:
        if predicate != "skos:exactMatch" or subject not in terms:
            continue
        if target.startswith("Orphanet:"):
            kept.add((subject, "ORPHA:" + target.split(":", 1)[1]))
        elif target.startswith("OMIM:"):
            kept.add((subject, target))
    return [
        {"mondo": mondo, "target": target, "source": SOURCE_MONDO} for mondo, target in sorted(kept)
    ]


def _mirrored_terms(root: Path) -> set[str]:
    terms: set[str] = set()
    for name in ("clingen_gene_validity.tsv", "gencc_submissions.tsv"):
        with (root / "mirrors" / name).open(encoding="utf-8") as handle:
            terms.update(row["disease"] for row in csv.DictReader(handle, delimiter="\t"))
    return terms


def jdbor_release(path: Path) -> tuple[str, str]:
    """The `version` and `date` an Orphadata file states on its root `<JDBOR>` element.

    Read from the first start event only, so a 50 MB file costs one chunk.
    Refuses a file whose root is not `JDBOR` or that states either one blank:
    a release recorded as nothing would cite as nothing.
    """
    for _, element in ET.iterparse(path, events=("start",)):
        if element.tag != "JDBOR":
            raise ValueError(f"{path}: root element is <{element.tag}>, not <JDBOR>")
        version, date = element.get("version", "").strip(), element.get("date", "").strip()
        if not version or not date:
            raise ValueError(f"{path}: <JDBOR> states no version or no date")
        return version, date
    raise ValueError(f"{path}: no XML element")


def record_release(sources: Path, source_id: str, digest: str, version: str, date: str) -> None:
    """Write one file's sha256 and header version and date into its `sources.yaml` entry.

    A line edit of the three keys inside `- id: <source_id>`'s block, so every
    comment in the file survives, as no YAML round-trip here would guarantee.
    Each key must already be there exactly once: the converter keeps the
    curated entry current, it does not create one, and a missing key is a
    malformed entry a human must look at. Values are written double-quoted
    (JSON string syntax is valid YAML), as the file writes `version`.
    """
    lines = sources.read_text(encoding="utf-8").splitlines(keepends=True)
    start = lines.index(f"  - id: {source_id}\n")
    end = next(
        (i for i in range(start + 1, len(lines)) if lines[i].startswith("  - id: ")), len(lines)
    )
    for key, value in (("version", digest), ("file_version", version), ("file_date", date)):
        found = [i for i in range(start, end) if lines[i].startswith(f"    {key}: ")]
        if len(found) != 1:
            raise ValueError(f"{source_id}: expected one `{key}:` line, found {len(found)}")
        lines[found[0]] = f"    {key}: {json.dumps(value)}\n"
    sources.write_text("".join(lines), encoding="utf-8")


def _write(path: Path, order: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    lines = ["\t".join(order), *("\t".join(row[column] for column in order) for row in rows)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str]) -> int:
    if not 4 <= len(argv) <= 5:
        print(__doc__, file=sys.stderr)
        return 2
    product4, product1, sssom, out = (Path(arg) for arg in argv[:4])
    root = Path(argv[4]) if len(argv) == 5 else Path.cwd()
    cardiovascular, known = _cardiovascular(root / HPO_RELEASE)

    phenotypes: list[dict[str, str]] = []
    unknown_hpo: set[str] = set()
    for _, element in ET.iterparse(product4):
        if element.tag == "HPODisorderSetStatus":
            for association in element.iter("HPODisorderAssociation"):
                hpo = _text(association, "HPO/HPOId")
                if hpo not in known:
                    unknown_hpo.add(hpo)
            phenotypes.extend(_phenotype_rows(element, cardiovascular))
            element.clear()
    unique = {tuple(row[c] for c in PHENOTYPE_ORDER): row for row in phenotypes}
    collapsed = len(phenotypes) - len(unique)
    phenotypes = [unique[key] for key in sorted(unique, key=lambda k: (k[0], k[2], k[4]))]

    omim: list[dict[str, str]] = []
    for _, element in ET.iterparse(product1):
        if element.tag == "Disorder":
            omim.extend(_omim_rows(element))
            element.clear()
    omim.sort(key=lambda row: (row["orpha"], row["omim"]))

    with sssom.open(encoding="utf-8") as handle:
        body = (line for line in handle if not line.startswith("#"))
        reader = csv.DictReader(body, delimiter="\t")
        triples = [(r["subject_id"], r["predicate_id"], r["object_id"]) for r in reader]
    xrefs = _xref_rows(triples, _mirrored_terms(root))

    out.mkdir(parents=True, exist_ok=True)
    _write(out / "orphanet_phenotypes.tsv", PHENOTYPE_ORDER, phenotypes)
    _write(out / "orphanet_omim.tsv", OMIM_ORDER, omim)
    _write(out / "mondo_xrefs.tsv", XREF_ORDER, xrefs)

    disorders = len({row["orpha"] for row in phenotypes})
    print(f"orphanet_phenotypes: {len(phenotypes)} row(s) over {disorders} disorder(s)")
    print(f"  {collapsed} exact-duplicate annotation(s) collapsed")
    if unknown_hpo:
        print(
            f"  WARNING: {len(unknown_hpo)} HPO id(s) unknown to {HPO_RELEASE}, so they cannot "
            f"be placed in the cardiovascular branch: {', '.join(sorted(unknown_hpo))}"
        )
    print(f"orphanet_omim: {len(omim)} row(s)")
    print(f"mondo_xrefs: {len(xrefs)} row(s)")
    sources = root / "mirrors" / "sources.yaml"
    for source_id, path in ((SOURCE_PRODUCT4, product4), (SOURCE_PRODUCT1, product1)):
        version, date = jdbor_release(path)
        digest = f"sha256:{hashlib.sha256(path.read_bytes()).hexdigest()}"
        record_release(sources, source_id, digest, version, date)
        print(f"{source_id}: {digest}, version {version!r}, date {date!r} -> {sources}")
    print(f"sha256 of {sssom.name}: {hashlib.sha256(sssom.read_bytes()).hexdigest()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
