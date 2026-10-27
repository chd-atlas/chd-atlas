# tests/unit/test_convert_farah2024.py
"""The Farah 2024 converter's contracts, checked on synthetic per-heart records.

Not written against the committed mirror: these call the builders directly, for
the reason `test_convert_cardoso_moreira.py` gives. Neither h5py nor anndata is
needed -- `extract_rows` is pure, and the file reading stays in its `main()`.

The risks pinned here:

- **The denominator of a CPM is the whole library**, summed over every gene the
  source measured. Dividing by the registry subset's total would inflate every
  value by a factor that depends on how many genes the registry holds.
- **A donor is a heart, not an age.** Two hearts at one gestational age are two
  donors; merging them would halve `n_donors` and hide the spread.
- **A quartile of fewer than three hearts is a coincidence of the sample**, so
  it is withheld rather than published.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))
sys.path.insert(0, str(REPO_ROOT / "scripts" / "hpc"))

from convert_farah2024 import (  # noqa: E402
    PseudobulkRow,
    main,
    map_symbols,
    summarise,
)
from farah_extract import extract_rows  # noqa: E402

from chd_atlas.genes import GeneRegistry  # noqa: E402

HEADER = "symbol\theart\tcell_type\tregion\tcounts\tcount_nonzero\tn_cells\tlibrary_size"


def _row(
    heart: str = "9pcw-1",
    *,
    counts: float = 50.0,
    library: float = 1_000_000.0,
    nonzero: int = 5,
    cells: int = 10,
    cell_type: str = "cardiomyocyte",
    region: str = "LV",
    symbol: str = "AAA",
) -> PseudobulkRow:
    return PseudobulkRow(symbol, heart, cell_type, region, counts, nonzero, cells, library)


def _summarise(rows: list[PseudobulkRow]) -> list[dict[str, object]]:
    return summarise(rows, {"AAA": "HGNC:1"}, dataset="farah_2024", source="farah.h5ad")


def test_cpm_divides_by_the_whole_library_not_the_registry_genes() -> None:
    (out,) = _summarise([_row(counts=50.0, library=1_000_000.0)])
    assert out["median_abundance"] == 50.0
    assert out["unit"] == "cpm"
    assert out["gene"] == "HGNC:1"
    assert out["dataset"] == "farah_2024"
    assert out["source"] == "farah.h5ad"
    assert list(out) == [
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
    ]
    (halved,) = _summarise([_row(counts=50.0, library=2_000_000.0)])
    assert halved["median_abundance"] == 25.0


def test_a_heart_is_age_and_replicate() -> None:
    rows = [_row("9pcw-1"), _row("9pcw-2")]
    (out,) = _summarise(rows)
    assert out["n_donors"] == 2
    assert out["n_cells"] == 20


def test_quartiles_need_three_hearts() -> None:
    (two,) = _summarise([_row("9pcw-1", counts=10), _row("9pcw-2", counts=20)])
    assert two["q25"] is None and two["q75"] is None
    assert two["median_abundance"] == 15.0
    hearts = [_row("9pcw-1", counts=10), _row("10pcw-1", counts=20), _row("11pcw-1", counts=40)]
    (three,) = _summarise(hearts)
    assert (three["q25"], three["median_abundance"], three["q75"]) == (15.0, 20.0, 30.0)
    assert three["n_donors"] == 3


def test_share_of_cells_pools_across_hearts() -> None:
    rows = [_row("9pcw-1", nonzero=3, cells=10), _row("9pcw-2", nonzero=1, cells=30)]
    (out,) = _summarise(rows)
    assert out["fraction_expressing"] == pytest.approx(0.1)


def test_groups_are_split_by_cell_type_and_region_and_sorted() -> None:
    rows = [
        _row(region="RV"),
        _row(region="LV"),
        _row(cell_type="aaa", region="LV"),
    ]
    out = _summarise(rows)
    assert [(r["cell_type"], r["region"]) for r in out] == [
        ("aaa", "LV"),
        ("cardiomyocyte", "LV"),
        ("cardiomyocyte", "RV"),
    ]


def test_a_heart_counted_twice_in_one_group_is_refused() -> None:
    with pytest.raises(SystemExit):
        _summarise([_row("9pcw-1"), _row("9pcw-1")])


def test_an_ambiguous_symbol_is_refused_and_a_missing_gene_named() -> None:
    registry = GeneRegistry(
        [
            {"hgnc_id": "HGNC:1", "symbol": "AAA", "aliases": "DUP", "prev_symbols": ""},
            {"hgnc_id": "HGNC:2", "symbol": "BBB", "aliases": "DUP", "prev_symbols": ""},
            {"hgnc_id": "HGNC:3", "symbol": "CCC", "aliases": "", "prev_symbols": "OLDC"},
        ]
    )
    mapping, report = map_symbols(["AAA", "DUP", "OLDC"], registry)
    assert mapping == {"AAA": "HGNC:1", "OLDC": "HGNC:3"}
    text = "\n".join(report)
    assert "DUP" in text and "ambiguous" in text
    assert "BBB" in text and "HGNC:2" in text
    assert "CCC" not in text

    with pytest.raises(SystemExit):
        map_symbols(["CCC", "OLDC"], registry)


def test_extract_rows_keeps_only_wanted_symbols_and_tags_the_region() -> None:
    groups = [
        {"Age": 9, "Rep": 1, "major_cell_class": "cardiomyocyte", "Region": "LV", "n_cells": 12},
        {"Age": 10, "Rep": 2, "major_cell_class": "fibroblast", "Region": "RA", "n_cells": 30},
    ]
    genes = ["AAA", "BBB", "CCC"]
    sums = [[1.0, 2.0, 7.0], [4.0, 4.0, 2.0]]
    nonzero = [[1, 2, 5], [3, 3, 2]]
    by_region = extract_rows(groups, genes, sums, nonzero, {"AAA", "CCC"}, None)
    assert by_region == [
        {
            "symbol": "AAA",
            "heart": "9pcw-1",
            "cell_type": "cardiomyocyte",
            "region": "LV",
            "counts": 1.0,
            "count_nonzero": 1,
            "n_cells": 12,
            "library_size": 10.0,
        },
        {
            "symbol": "CCC",
            "heart": "9pcw-1",
            "cell_type": "cardiomyocyte",
            "region": "LV",
            "counts": 7.0,
            "count_nonzero": 5,
            "n_cells": 12,
            "library_size": 10.0,
        },
        {
            "symbol": "AAA",
            "heart": "10pcw-2",
            "cell_type": "fibroblast",
            "region": "RA",
            "counts": 4.0,
            "count_nonzero": 3,
            "n_cells": 30,
            "library_size": 10.0,
        },
        {
            "symbol": "CCC",
            "heart": "10pcw-2",
            "cell_type": "fibroblast",
            "region": "RA",
            "counts": 2.0,
            "count_nonzero": 2,
            "n_cells": 30,
            "library_size": 10.0,
        },
    ]
    pooled = extract_rows(groups, genes, sums, nonzero, {"BBB"}, "all")
    assert {r["region"] for r in pooled} == {"all"}
    assert [r["library_size"] for r in pooled] == [10.0, 10.0]


def _registry_file(tmp_path: Path) -> Path:
    path = tmp_path / "genes.tsv"
    path.write_text(
        "hgnc_id\tsymbol\tname\taliases\tprev_symbols\tensembl_gene\tncbi_gene\tlocus\tuniprot\t"
        "mane_select\nHGNC:1\tAAA\tn\t\t\t\t\t\t\t\n"
    )
    return path


def test_main_writes_the_mirror_columns(tmp_path: Path) -> None:
    pseudobulk = tmp_path / "pb.tsv"
    pseudobulk.write_text(
        HEADER + "\nAAA\t9pcw-1\tcardiomyocyte\tLV\t50\t5\t10\t1000000\n", encoding="utf-8"
    )
    out = tmp_path / "out" / "cell_type_profiles.tsv"
    argv = [
        "--pseudobulk",
        str(pseudobulk),
        "--registry",
        str(_registry_file(tmp_path)),
        "--dataset",
        "farah_2024",
        "--source",
        "farah.h5ad",
        "--out",
        str(out),
    ]
    assert main(argv) == 0
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == [
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
    ]
    assert lines[1].split("\t") == [
        "farah_2024",
        "HGNC:1",
        "cardiomyocyte",
        "LV",
        "50.0",
        "cpm",
        "",
        "",
        "1",
        "10",
        "0.5",
        "farah.h5ad",
    ]


@pytest.mark.parametrize(
    "line",
    [
        "AAA\t9pcw-1\tcardiomyocyte\tLV\tNaN\t5\t10\t1000000",
        "AAA\t9pcw-1\tcardiomyocyte\tLV\t50\t5\t10\t0",
        "AAA\t9pcw-1\tcardiomyocyte\tLV\t50\t11\t10\t1000000",
        "AAA\t9pcw-1\tcardiomyocyte\tLV\t50\t5\t10",
        "AAA\t9pcw-1\tcardiomyocyte\tleft atrium\t50\t5\t10\t1000000",
    ],
)
def test_nothing_is_written_when_a_row_is_malformed(tmp_path: Path, line: str) -> None:
    pseudobulk = tmp_path / "pb.tsv"
    good = "AAA\t9pcw-2\tcardiomyocyte\tLV\t50\t5\t10\t1000000"
    pseudobulk.write_text(f"{HEADER}\n{good}\n{line}\n", encoding="utf-8")
    out = tmp_path / "out.tsv"
    argv = [
        "--pseudobulk",
        str(pseudobulk),
        "--registry",
        str(_registry_file(tmp_path)),
        "--dataset",
        "d",
        "--source",
        "s",
        "--out",
        str(out),
    ]
    with pytest.raises(SystemExit):
        main(argv)
    assert not out.exists()
