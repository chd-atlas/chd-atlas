# src/chd_atlas/build/cell_type_chart.py
"""The gene page's cell-type block: one dot-and-whisker row per declared cell type.

Static SVG, no script. One shared log axis from `_AXIS_MIN` to `_AXIS_MAX` CPM
with a gridline at each decade and a shaded band below the dataset's detection
floor. A row's dot sits at the pooled median, its area is proportional to the
share of the cell type's cells with at least one read, and its whisker spans the
quartiles across hearts. A "By region" fold repeats the cell types against the
dataset's declared regions.

**Every declared cell type gets a row, in the dataset's own order.**
`cell_types.py` omits a cell type no heart reached 10 cells of, so this module
walks `Dataset.cell_types` and writes "no data" where `CellTypeDataset` has no
entry: an absent row would read as a measured zero. A region the payload has no
row for draws "–" for the same reason. A median at or below `_AXIS_MIN`
(including 0) is drawn at the axis minimum; its number column still says what it
is.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Final

from chd_atlas.build.cell_types import CellTypeDataset, CellTypeEntry, CellTypeValues
from chd_atlas.models.dataset import Dataset

_HEADING: Final = "Cell types (Farah 2024)"
_CAPTION: Final = "Log axis, 0.1 to 10,000 CPM; a dot's size is its share of cells."
_READING_SUMMARY: Final = "How to read this chart"
_READING_TEXT: Final = (
    "The dot is the median across hearts; the line spans the middle half of hearts (shown "
    "from three hearts). A heart counts only with 10 or more cells of that type. The shaded "
    'band is below this atlas\'s 1 CPM detection floor. "No data" means no heart reached 10 '
    "cells; it is not a low value. Values pool 9 to 15 weeks. "
)
_READING_LINK: Final = '<a href="../methods.html#expression">How to read cell-type figures</a>.'
_NO_DATA: Final = "no data"
_NO_REGION: Final = "–"

_AXIS_MIN: Final = 0.1
_AXIS_MAX: Final = 10_000.0
_LOG_MIN: Final = math.log10(_AXIS_MIN)
_LOG_MAX: Final = math.log10(_AXIS_MAX)

_ROW_H: Final = 26
_TOP: Final = 40
_LEFT: Final = 190
_AXIS_W: Final = 430
_NUMBERS_W: Final = 230
_REGION_COL_W: Final = 62


def _round(value: float) -> int:
    """Half up, as the figures' reader expects; `round` would go to even."""
    return int(value + 0.5)


def format_median(median: float) -> str:
    """`0` reads `0`; from 100 up an integer with thousands separators; from 10 an
    integer; from 1 one decimal; below 1 two decimals."""
    if median == 0:
        return "0"
    if median >= 100:
        return f"{_round(median):,}"
    if median >= 10:
        return str(_round(median))
    if median >= 1:
        return f"{median:.1f}"
    return f"{median:.2f}"


def _numbers(values: CellTypeValues) -> str:
    hearts = values["n_donors"]
    noun = "heart" if hearts == 1 else "hearts"
    share = _round(values["fraction_expressing"] * 100)
    return f"{format_median(values['median'])} CPM · {share}% · {hearts} {noun}"


def _title(label: str, values: CellTypeValues, where: str = "") -> str:
    spread = (
        ""
        if values["q25"] is None or values["q75"] is None
        else f", quartiles {format_median(values['q25'])} to {format_median(values['q75'])}"
    )
    return (
        f"{label}{where}: median {format_median(values['median'])} CPM{spread}, "
        f"{_round(values['fraction_expressing'] * 100)}% of {values['n_cells']:,} cells, "
        f"{values['n_donors']} hearts"
    )


def _x(value: float, left: float, width: float) -> float:
    """Log position of `value`; anything at or below the axis minimum clamps to it."""
    clamped = min(max(value, _AXIS_MIN), _AXIS_MAX)
    return left + (math.log10(clamped) - _LOG_MIN) / (_LOG_MAX - _LOG_MIN) * width


def _radius(fraction: float, scale: float) -> float:
    return 2.5 + scale * math.sqrt(min(max(fraction, 0.0), 1.0))


def _text(x: float, y: float, css: str, body: str, extra: str = "") -> str:
    return f'<text x="{x:.1f}" y="{y:.1f}" class="{css}"{extra}>{body}</text>'


def _rows(dataset: Dataset) -> list[tuple[str, str | None, str | None]]:
    """(compartment, token, label) per drawn row; a header row has no token."""
    rows: list[tuple[str, str | None, str | None]] = []
    previous: str | None = None
    for declared in dataset.cell_types:
        if declared.compartment != previous:
            rows.append((declared.compartment, None, None))
            previous = declared.compartment
        rows.append((declared.compartment, declared.token, declared.label))
    return rows


def _svg(width: float, height: float, label: str, body: str) -> str:
    return (
        f'<div class="ct-scroll"><svg class="ct-svg" viewBox="0 0 {width:.0f} {height:.0f}" '
        f'width="{width:.0f}" height="{height:.0f}" role="img" aria-label="{label}">'
        f"{body}</svg></div>"
    )


def _axis(top: float, bottom: float) -> str:
    parts = []
    for power in range(int(_LOG_MIN), int(_LOG_MAX) + 1):
        value = 10.0**power
        x = _x(value, _LEFT, _AXIS_W)
        parts.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{top}" y2="{bottom}" class="ct-grid"/>')
        label = "0.1" if value < 1 else f"{int(value):,}"
        parts.append(_text(x, top - 6, "ct-small ct-mid", label))
    parts.append(
        _text(_x(_AXIS_MAX, _LEFT, _AXIS_W), top - 20, "ct-small ct-end", "CPM, log scale")
    )
    return "".join(parts)


def _chart(entry: CellTypeDataset, dataset: Dataset, by_token: Mapping[str, CellTypeEntry]) -> str:
    rows = _rows(dataset)
    bottom = _TOP + len(rows) * _ROW_H
    column = _LEFT + _AXIS_W + 40
    floor_x = _x(entry["detection_floor"], _LEFT, _AXIS_W)
    left_x = _x(_AXIS_MIN, _LEFT, _AXIS_W)
    parts = [
        f'<rect x="{left_x:.1f}" y="{_TOP}" width="{floor_x - left_x:.1f}" '
        f'height="{bottom - _TOP}" class="ct-below"/>',
        _axis(_TOP, bottom),
        f'<line x1="{floor_x:.1f}" x2="{floor_x:.1f}" y1="{_TOP}" y2="{bottom}" class="ct-floor"/>',
    ]
    for index, (compartment, token, label) in enumerate(rows):
        y = _TOP + index * _ROW_H + _ROW_H / 2
        if token is None or label is None:
            parts.append(_text(0, y + 4, "ct-group", compartment))
            continue
        parts.append(_text(12, y + 4, "ct-label", label))
        found = by_token.get(token)
        if found is None:
            parts.append(_text(_LEFT + 4, y + 4, "ct-muted", _NO_DATA))
            continue
        values = found["pooled"]
        if values["q25"] is not None and values["q75"] is not None:
            parts.append(
                f'<line x1="{_x(values["q25"], _LEFT, _AXIS_W):.1f}" '
                f'x2="{_x(values["q75"], _LEFT, _AXIS_W):.1f}" y1="{y:.1f}" y2="{y:.1f}" '
                f'class="ct-whisker" data-cell-type="{token}"/>'
            )
        parts.append(
            f'<circle cx="{_x(values["median"], _LEFT, _AXIS_W):.1f}" cy="{y:.1f}" '
            f'r="{_radius(values["fraction_expressing"], 7.5):.2f}" class="ct-dot" '
            f'data-cell-type="{token}"><title>{_title(label, values)}</title></circle>'
        )
        parts.append(_text(column, y + 4, "ct-num", _numbers(values), f' data-cell-type="{token}"'))
    parts.append(_text(column, _TOP - 6, "ct-small", "median · share · hearts"))
    height = bottom + 28
    return _svg(
        _LEFT + _AXIS_W + _NUMBERS_W,
        height,
        f"{dataset.id} expression by cell type",
        "".join(parts),
    )


def _opacity(values: CellTypeValues) -> float:
    """Fill opacity stepped by CPM decade; a value below the floor is an open ring."""
    median = values["median"]
    for bound, opacity in ((10, 0.3), (100, 0.55), (1000, 0.8)):
        if median < bound:
            return opacity
    return 1.0


_LEGEND: Final = (
    ("below floor", None),
    ("under 10", 0.3),
    ("10–100", 0.55),
    ("100–1,000", 0.8),
    ("1,000+ CPM", 1.0),
)


def _region_grid(dataset: Dataset, by_token: Mapping[str, CellTypeEntry]) -> str:
    rows = _rows(dataset)
    top = 30
    parts = []
    for column, region in enumerate(dataset.regions):
        cx = _LEFT + column * _REGION_COL_W + _REGION_COL_W / 2
        parts.append(_text(cx, top - 10, "ct-small ct-mid", region))
    for index, (compartment, token, label) in enumerate(rows):
        y = top + index * _ROW_H + _ROW_H / 2
        if token is None or label is None:
            parts.append(_text(0, y + 4, "ct-group", compartment))
            continue
        parts.append(_text(12, y + 4, "ct-label", label))
        found = by_token.get(token)
        by_region = {} if found is None else {r["region"]: r for r in found["regions"]}
        for column, region in enumerate(dataset.regions):
            cx = _LEFT + column * _REGION_COL_W + _REGION_COL_W / 2
            cell = by_region.get(region)
            tag = f' data-cell-type="{token}" data-region="{region}"'
            if cell is None:
                parts.append(_text(cx, y + 4, "ct-muted ct-mid", _NO_REGION, tag))
                continue
            css = "ct-ring" if cell["below_floor"] else "ct-dot"
            opacity = "" if cell["below_floor"] else f' fill-opacity="{_opacity(cell)}"'
            radius = _radius(cell["fraction_expressing"], 8)
            title = _title(label, cell, f" ({region})")
            parts.append(
                f'<circle cx="{cx:.1f}" cy="{y:.1f}" r="{radius:.2f}" class="{css}"{opacity}{tag}>'
                f"<title>{title}</title></circle>"
            )
    legend_y = top + len(rows) * _ROW_H + 18
    for index, (name, level) in enumerate(_LEGEND):
        cx = 40 + index * 100
        mark = 'class="ct-ring"' if level is None else f'class="ct-dot" fill-opacity="{level}"'
        parts.append(f'<circle cx="{cx}" cy="{legend_y}" r="6" {mark}/>')
        parts.append(_text(cx + 10, legend_y + 4, "ct-small", name))
    return _svg(
        _LEFT + len(dataset.regions) * _REGION_COL_W + 10,
        legend_y + 24,
        f"{dataset.id} expression by cell type and region",
        "".join(parts),
    )


def _tau_line(entry: CellTypeDataset, dataset: Dataset) -> str:
    spec = entry["specificity"]
    if spec is None:
        return ""
    labels = {c.token: c.label for c in dataset.cell_types}
    head = f"τ = {spec['tau']:.2f} ({spec['n_cell_types']} cell types, log2(x+1)), "
    highest = spec["highest_in"]
    tail = "with no single highest." if highest is None else f"highest in {labels[highest]}."
    return f'<p class="ct-tau">{head}{tail}</p>'


def cell_type_block(entry: CellTypeDataset, dataset: Dataset) -> str:
    """One dataset's cell-type block: heading, caption, chart, and two folds."""
    by_token = {item["cell_type"]: item for item in entry["cell_types"]}
    return (
        '<section class="cell-types">'
        f"<h3>{_HEADING}</h3>"
        f'<p class="ct-caption">{_CAPTION}</p>'
        + _tau_line(entry, dataset)
        + _chart(entry, dataset, by_token)
        + "<details><summary>By region</summary>"
        + _region_grid(dataset, by_token)
        + "</details>"
        + f"<details><summary>{_READING_SUMMARY}</summary>"
        + f"<p>{_READING_TEXT}{_READING_LINK}</p></details>"
        + "</section>"
    )
