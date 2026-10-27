"""Cross-checks between ``mirrors/cell_type_profiles/`` and the dataset records.

`validate_table` checks each shard against its own column declarations and
cannot see the dataset a row belongs to. A cell-type row is only meaningful
against that dataset's declarations: the cell types it names, the regions it
sampled, and how many donors it had. Seven ways a row can disagree with them:

- **CTP001** the row's dataset has no `cell_type_profile` record.
- **CTP002** the cell type is not one the dataset declares.
- **CTP003** the region is neither one the dataset declares nor `all`.
- **CTP004** a region row has no `all` row for the same
  (dataset, gene, cell_type); `all` is the pooled value every page falls back to.
- **CTP005** `q25 <= median_abundance <= q75` does not hold, or exactly one of
  q25 and q75 is null (an interval with one end is not an interval).
- **CTP006** `n_donors` exceeds the dataset's `n_samples`.
- **CTP007** the gene is not in the gene registry. Skipped when `known_genes`
  is None, which means the registry did not load and the runner has already
  said why; an empty-but-loaded registry is not None and flags every gene.

CTP002 and CTP003 need the dataset's declarations, so they are skipped for a
row whose dataset is unrecorded (CTP001 already names that). A shard that
cannot be read, or lacks a column a rule needs, is skipped: TBL001/TBL002
already report it. Each rule reports once per offending value per shard, not
once per row.
"""

from __future__ import annotations

from collections.abc import Collection, Sequence
from pathlib import Path

from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.models.dataset import Dataset
from chd_atlas.tables import CELL_TYPE_PROFILES, mirror_paths, read_table
from chd_atlas.vocab import DatasetDesign

_COLUMNS = [
    "dataset",
    "gene",
    "cell_type",
    "region",
    "median_abundance",
    "q25",
    "q75",
    "n_donors",
]
_POOLED = "all"


def validate_cell_type_profiles(
    root: Path,
    datasets: Sequence[Dataset],
    known_genes: Collection[str] | None,
) -> list[ValidationIssue]:
    """CTP001-CTP007 over every `cell_type_profiles` shard under `root`."""
    records: dict[str, Dataset] = {
        str(d.id): d for d in datasets if d.design is DatasetDesign.CELL_TYPE_PROFILE
    }
    issues: list[ValidationIssue] = []
    for path, name in mirror_paths(root):
        if name != CELL_TYPE_PROFILES.name:
            continue
        frame, _ = read_table(path, CELL_TYPE_PROFILES)
        if frame is None or not set(_COLUMNS) <= set(frame.columns):
            continue
        rows = list(frame.select(_COLUMNS).iter_rows(named=True))
        issues.extend(_shard_issues(str(path), rows, records, known_genes))
    return issues


def _shard_issues(
    location: str,
    rows: list[dict[str, object]],
    records: dict[str, Dataset],
    known_genes: Collection[str] | None,
) -> list[ValidationIssue]:
    pooled = {(r["dataset"], r["gene"], r["cell_type"]) for r in rows if r["region"] == _POOLED}
    found: dict[tuple[str, str], str] = {}

    def report(code: str, key: str, message: str) -> None:
        found.setdefault((code, key), message)

    for row in rows:
        dataset, gene = str(row["dataset"]), str(row["gene"])
        cell_type, region = str(row["cell_type"]), str(row["region"])
        record = records.get(dataset)
        if record is None:
            report(
                "CTP001",
                dataset,
                f"dataset '{dataset}' has no cell_type_profile record in datasets/",
            )
        else:
            if cell_type not in {c.token for c in record.cell_types}:
                report(
                    "CTP002",
                    f"{dataset}/{cell_type}",
                    f"cell type '{cell_type}' is not declared by dataset '{dataset}'",
                )
            if region != _POOLED and region not in record.regions:
                report(
                    "CTP003",
                    f"{dataset}/{region}",
                    f"region '{region}' is not declared by dataset '{dataset}' "
                    f"(declared: {sorted(record.regions)}, plus '{_POOLED}')",
                )
            n_donors = row["n_donors"]
            if isinstance(n_donors, int) and n_donors > record.n_samples:
                report(
                    "CTP006",
                    f"{dataset}/{gene}/{cell_type}/{region}",
                    f"{gene}/{cell_type}/{region}: n_donors {n_donors} exceeds "
                    f"dataset '{dataset}' n_samples {record.n_samples}",
                )
        if region != _POOLED and (dataset, gene, cell_type) not in pooled:
            report(
                "CTP004",
                f"{dataset}/{gene}/{cell_type}",
                f"{gene}/{cell_type} in dataset '{dataset}' has region rows but no '{_POOLED}' row",
            )
        broken = _quartile_problem(row)
        if broken:
            report(
                "CTP005",
                f"{dataset}/{gene}/{cell_type}/{region}",
                f"{gene}/{cell_type}/{region} in dataset '{dataset}': {broken}",
            )
        if known_genes is not None and gene not in known_genes:
            report("CTP007", gene, f"gene '{gene}' is not in the gene registry")
    return [
        ValidationIssue(code, Severity.ERROR, location, message)
        for (code, _), message in sorted(found.items())
    ]


def _quartile_problem(row: dict[str, object]) -> str | None:
    """What is wrong with the row's q25/median/q75, or None if it is fine."""
    q25, median, q75 = row["q25"], row["median_abundance"], row["q75"]
    if (q25 is None) != (q75 is None):
        return "exactly one of q25 and q75 is null"
    if q25 is None or q75 is None or median is None:
        return None
    if not (q25 <= median <= q75):  # type: ignore[operator]
        return f"q25 {q25} <= median {median} <= q75 {q75} does not hold"
    return None
