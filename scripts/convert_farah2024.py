#!/usr/bin/env python3
"""Per-heart pseudobulk of Farah et al. 2024 -> mirrors/cell_type_profiles.

Input is the TSV `scripts/hpc/farah_extract.py` writes on the cluster: one row
per (gene, heart, cell type, region) with the group's summed counts, the number
of its cells with a non-zero count, its cell count, and its library size (the
`sum` layer summed over every gene the source measured, not only the registry's).

Per (gene, cell_type, region) this publishes the median across hearts of
counts-per-million, with quartiles when at least three hearts support them, the
number of hearts, the cells pooled, and the share of those cells expressing the
gene. A heart is `{age}pcw-{replicate}`: two hearts of one gestational age are
two donors.

Nothing is guessed. A symbol that more than one gene claims is dropped and
reported; two symbols that resolve to one gene, a heart counted twice in one
group, and any malformed input row stop the run before the output is touched.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import math
import sys
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import NamedTuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

from convert_cardoso_moreira import quantile_linear  # noqa: E402

from chd_atlas.fs import write_bytes_atomically  # noqa: E402
from chd_atlas.genes import GeneRegistry  # noqa: E402

UNIT = "cpm"
REGIONS = frozenset({"LA", "RA", "LV", "RV", "IVS", "all"})
MIN_HEARTS_FOR_QUARTILES = 3
COLUMNS = (
    "dataset",
    "gene",
    "cell_type",
    "region",
    "median_abundance",
    "unit",
    "q25",
    "q75",
    "n_donors",
    "n_cells",
    "fraction_expressing",
    "source",
)
INPUT_COLUMNS = (
    "symbol",
    "heart",
    "cell_type",
    "region",
    "counts",
    "count_nonzero",
    "n_cells",
    "library_size",
)


class PseudobulkRow(NamedTuple):
    symbol: str
    heart: str
    cell_type: str
    region: str
    counts: float
    count_nonzero: int
    n_cells: int
    library_size: float


def _fail(message: str) -> SystemExit:
    return SystemExit(f"ERROR: {message}")


def _finite(text: str, what: str, where: str) -> float:
    try:
        value = float(text)
    except ValueError:
        raise _fail(f"{where}: {what} {text!r} is not a number") from None
    if not math.isfinite(value):
        raise _fail(f"{where}: {what} {text!r} is not finite")
    return value


def read_pseudobulk(path: Path) -> list[PseudobulkRow]:
    """Every row of the extract, or SystemExit naming the first bad one."""
    rows: list[PseudobulkRow] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None or tuple(header) != INPUT_COLUMNS:
            raise _fail(f"{path}: header is {header!r}, expected {list(INPUT_COLUMNS)}")
        for number, record in enumerate(reader, start=2):
            where = f"{path}:{number}"
            if len(record) != len(INPUT_COLUMNS):
                raise _fail(f"{where}: {len(record)} fields, expected {len(INPUT_COLUMNS)}")
            symbol, heart, cell_type, region, counts, nonzero, cells, library = record
            if not (symbol and heart and cell_type):
                raise _fail(f"{where}: blank symbol, heart or cell_type")
            if region not in REGIONS:
                raise _fail(f"{where}: region {region!r} is not one of {sorted(REGIONS)}")
            count_value = _finite(counts, "counts", where)
            library_value = _finite(library, "library_size", where)
            nonzero_value = _finite(nonzero, "count_nonzero", where)
            cells_value = _finite(cells, "n_cells", where)
            if count_value < 0:
                raise _fail(f"{where}: counts {counts!r} is negative")
            if library_value <= 0:
                raise _fail(f"{where}: library_size {library!r} must be positive")
            if nonzero_value != int(nonzero_value) or cells_value != int(cells_value):
                raise _fail(f"{where}: count_nonzero and n_cells must be whole numbers")
            if not 0 <= nonzero_value <= cells_value or cells_value < 1:
                raise _fail(
                    f"{where}: count_nonzero {nonzero!r} must lie within 0..n_cells {cells!r}"
                )
            rows.append(
                PseudobulkRow(
                    symbol,
                    heart,
                    cell_type,
                    region,
                    count_value,
                    int(nonzero_value),
                    int(cells_value),
                    library_value,
                )
            )
    return rows


def map_symbols(symbols: Iterable[str], registry: GeneRegistry) -> tuple[dict[str, str], list[str]]:
    """symbol -> HGNC id, plus report lines.

    A symbol that no gene claims, or that several do, is left out and reported;
    it is never resolved to one candidate. Two symbols resolving to the same id
    would pool one gene's counts under two names, so that stops the run. Every
    registry gene that no symbol reached is named, so a gene absent from the
    source is visible rather than missing without comment.
    """
    mapping: dict[str, str] = {}
    report: list[str] = []
    claimed_by: dict[str, str] = {}
    for symbol in sorted(set(symbols)):
        resolved = registry.resolve(symbol)
        if not resolved.ok:
            kind = "ambiguous" if resolved.candidates else "unknown"
            report.append(f"{symbol}: {kind} {list(resolved.candidates)}, dropped")
            continue
        assert resolved.gene is not None
        if resolved.gene in claimed_by:
            raise _fail(
                f"{claimed_by[resolved.gene]!r} and {symbol!r} both resolve to {resolved.gene}"
            )
        claimed_by[resolved.gene] = symbol
        mapping[symbol] = resolved.gene
    for approved in registry.approved_symbols():
        gene = registry.resolve(approved).gene
        if gene is not None and gene not in claimed_by:
            report.append(f"{gene} ({approved}): no symbol in the pseudobulk")
    return mapping, report


def summarise(
    rows: Sequence[PseudobulkRow],
    symbol_to_gene: Mapping[str, str],
    *,
    dataset: str,
    source: str,
) -> list[dict[str, object]]:
    """One row per (gene, cell_type, region), sorted to the table's sort key.

    Rows whose symbol is not in `symbol_to_gene` are skipped.
    """
    cpm: dict[tuple[str, str, str], dict[str, float]] = defaultdict(dict)
    nonzero: dict[tuple[str, str, str], int] = defaultdict(int)
    cells: dict[tuple[str, str, str], int] = defaultdict(int)
    for row in rows:
        gene = symbol_to_gene.get(row.symbol)
        if gene is None:
            continue
        key = (gene, row.cell_type, row.region)
        if row.heart in cpm[key]:
            raise _fail(f"heart {row.heart} appears twice for {key}")
        cpm[key][row.heart] = row.counts / row.library_size * 1e6
        nonzero[key] += row.count_nonzero
        cells[key] += row.n_cells

    out: list[dict[str, object]] = []
    for key in sorted(cpm, key=lambda k: (dataset, k[0], k[1], k[2])):
        gene, cell_type, region = key
        values = sorted(cpm[key].values())
        enough = len(values) >= MIN_HEARTS_FOR_QUARTILES
        out.append(
            {
                "dataset": dataset,
                "gene": gene,
                "cell_type": cell_type,
                "region": region,
                "median_abundance": quantile_linear(values, 50),
                "unit": UNIT,
                "q25": quantile_linear(values, 25) if enough else None,
                "q75": quantile_linear(values, 75) if enough else None,
                "n_donors": len(values),
                "n_cells": cells[key],
                "fraction_expressing": nonzero[key] / cells[key],
                "source": source,
            }
        )
    return out


def _render(rows: Sequence[Mapping[str, object]]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n", extrasaction="raise"
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({k: ("" if v is None else v) for k, v in row.items()})
    return buffer.getvalue().encode("utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pseudobulk", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    digest = hashlib.sha256(args.pseudobulk.read_bytes()).hexdigest()
    print(f"pseudobulk sha256: {digest}", file=sys.stderr)
    rows = read_pseudobulk(args.pseudobulk)
    with args.registry.open(newline="", encoding="utf-8") as handle:
        registry = GeneRegistry(csv.DictReader(handle, delimiter="\t"))
    mapping, report = map_symbols({r.symbol for r in rows}, registry)
    for line in report:
        print(f"  {line}", file=sys.stderr)
    summary = summarise(rows, mapping, dataset=args.dataset, source=args.source)
    if not summary:
        raise _fail("no pseudobulk row maps to a registry gene; nothing to write")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_bytes_atomically(args.out, _render(summary))
    print(f"input rows: {len(rows)}  genes mapped: {len(mapping)}", file=sys.stderr)
    print(f"output rows: {len(summary)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
