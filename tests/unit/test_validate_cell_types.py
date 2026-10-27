# tests/unit/test_validate_cell_types.py
"""CTP001-CTP007: a cell-type row checked against its dataset record.

Each case breaks exactly one rule on top of a clean base, so a failure names
the rule that stopped isolating.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from chd_atlas.models.dataset import Dataset
from chd_atlas.validate.cell_types import validate_cell_type_profiles

_HEADER = (
    "dataset\tgene\tcell_type\tregion\tmedian_abundance\tunit\tq25\tq75\t"
    "n_donors\tn_cells\tfraction_expressing\tsource\n"
)


def _dataset() -> Dataset:
    return Dataset(
        id="hoc",
        archive="ucsc_cell_browser",
        technology="scrnaseq",
        design="cell_type_profile",
        tissue="Heart",
        developmental_stage="adult",
        organism="NCBITaxon:9606",
        n_samples=10,
        licence="CC BY 4.0",
        cell_types=[
            {"token": "Cardiomyocyte", "label": "Cardiomyocyte", "compartment": "myocardial"},
        ],
        regions=["LV", "RV"],
        detection_floor=0.1,
        floor_source="source methods, section 2",
        quantile_estimator="linear",
    )


def _row(
    *,
    dataset: str = "hoc",
    gene: str = "HGNC:1100",
    cell_type: str = "Cardiomyocyte",
    region: str = "all",
    median: str = "2.0",
    q25: str = "1.0",
    q75: str = "3.0",
    n_donors: str = "4",
) -> str:
    return (
        f"{dataset}\t{gene}\t{cell_type}\t{region}\t{median}\tcpm\t{q25}\t{q75}\t"
        f"{n_donors}\t900\t0.4\ts\n"
    )


def _run(tmp_path: Path, rows: list[str], genes: set[str] | None) -> list[str]:
    shard = tmp_path / "mirrors" / "cell_type_profiles"
    shard.mkdir(parents=True)
    (shard / "hoc.tsv").write_text(_HEADER + "".join(rows), encoding="utf-8")
    issues = validate_cell_type_profiles(tmp_path, [_dataset()], genes)
    return [issue.code for issue in issues]


_GENES = {"HGNC:1100"}


@pytest.mark.parametrize(
    ("rows", "genes", "expected"),
    [
        pytest.param([_row()], _GENES, [], id="clean"),
        pytest.param([_row(dataset="other")], _GENES, ["CTP001"], id="CTP001"),
        pytest.param([_row(cell_type="Mast")], _GENES, ["CTP002"], id="CTP002"),
        pytest.param([_row(), _row(region="RA")], _GENES, ["CTP003"], id="CTP003"),
        pytest.param([_row(region="LV")], _GENES, ["CTP004"], id="CTP004"),
        pytest.param([_row(median="5.0")], _GENES, ["CTP005"], id="CTP005"),
        pytest.param([_row(n_donors="11")], _GENES, ["CTP006"], id="CTP006"),
        pytest.param([_row()], {"HGNC:9"}, ["CTP007"], id="CTP007"),
        pytest.param([_row(q25="")], _GENES, ["CTP005"], id="one-quartile-null"),
        pytest.param([_row()], None, [], id="registry-unloaded-skips-CTP007"),
    ],
)
def test_each_cell_type_check_isolates_its_own_error(
    tmp_path: Path, rows: list[str], genes: set[str] | None, expected: list[str]
) -> None:
    assert _run(tmp_path, rows, genes) == expected
