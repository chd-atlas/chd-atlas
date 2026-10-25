# tests/unit/test_convert_hgnc.py
"""What the HGNC converter must not lose on the way into the mirror.

Every column here is one a downstream join reads. `uniprot` in particular is
what any future PTM or proteomics join keys on, so a gene silently arriving
with a null accession is evidence loss rather than a cosmetic gap.
"""

from pathlib import Path

import polars as pl
import pytest
from scripts.convert_hgnc import convert

_HEADER = (
    "hgnc_id\tsymbol\tname\talias_symbol\tensembl_gene_id\tentrez_id\t"
    "location\tuniprot_ids\tmane_select\tstatus\tprev_symbol"
)
_WIDTH = _HEADER.count("\t") + 1


def _source(tmp_path: Path, *rows: str) -> Path:
    """The fixture source file, padding a row that declares no previous symbol.

    `prev_symbol` is last in `_HEADER` and every row here predates it, so a row
    one field short means "this gene has no retired name" rather than a typo.
    Padded rather than spelled out on each row because `convert` selects
    columns by name; the position carries no meaning and a trailing tab on four
    unrelated fixtures would just be noise. **Only one field of slack** -- a row
    two short is a mistake and raises here rather than silently mirroring a
    shifted value.
    """
    padded: list[str] = []
    for row in rows:
        missing = _WIDTH - (row.count("\t") + 1)
        if missing not in (0, 1):
            raise ValueError(f"row has {_WIDTH - missing} of {_WIDTH} fields: {row!r}")
        padded.append(row + "\t" * missing)
    path = tmp_path / "hgnc.tsv"
    path.write_text("\n".join((_HEADER, *padded)) + "\n", encoding="utf-8")
    return path


def test_a_retired_symbol_is_mirrored_verbatim_and_never_folded_into_aliases(tmp_path):
    """`prev_symbols` is the column issue #33 asked for and the reason the
    standing `GEN003` on HGNC:1152 could not be cleared by regenerating the
    mirror: `BVES` is that gene's *previous* symbol, and the converter read
    `alias_symbol` alone.

    Two genes rather than one, because a fixture whose rows all carry the same
    shape cannot tell "copied the column" from "copied every column"
    (CLAUDE.md 4.36). One gene retires two names, so the pipe survives; the
    other retires none, so the null does.

    The `aliases` assertions are the load-bearing half. `genes.py` ranks a live
    alias *above* a retired name, and that precedence is what makes `ODD`
    resolve to GJA1 rather than to nothing -- a converter that merged the two
    columns would erase the distinction upstream of the module that depends on
    it, and every resolution test would still pass.
    """
    source = _source(
        tmp_path,
        "HGNC:1152\tPOPDC1\tpopeye domain cAMP effector 1\tPOP1|HBVES\tENSG00000112276\t"
        '11149\t6q21\tQ8NE79\t"ENST00000265242.9|NM_001199563.2"\tApproved\tBVES|POPDC',
        "HGNC:11604\tTBX5\tT-box transcription factor 5\t\tENSG00000089225\t6910\t"
        '12q24.21\tQ99593\t"ENST00000405440.7|NM_181486.4"\tApproved',
    )
    out = tmp_path / "genes.tsv"
    convert(source, out, keep={"HGNC:1152", "HGNC:11604"})

    rows = {
        row["hgnc_id"]: row
        for row in pl.read_csv(out, separator="\t", infer_schema_length=0).iter_rows(named=True)
    }
    assert rows["HGNC:1152"]["prev_symbols"] == "BVES|POPDC"
    assert rows["HGNC:1152"]["aliases"] == "POP1|HBVES"
    assert rows["HGNC:11604"]["prev_symbols"] is None
    assert rows["HGNC:11604"]["aliases"] is None


def test_the_refseq_half_of_mane_select_is_kept_and_the_ensembl_half_dropped(tmp_path):
    source = _source(
        tmp_path,
        "HGNC:11604\tTBX5\tT-box transcription factor 5\t\tENSG00000089225\t6910\t"
        '12q24.21\tQ99593\t"ENST00000405440.7|NM_181486.4"\tApproved',
    )
    out = tmp_path / "genes.tsv"
    convert(source, out, keep={"HGNC:11604"})

    row = pl.read_csv(out, separator="\t", infer_schema_length=0).row(0, named=True)
    assert row["mane_select"] == "NM_181486.4"
    assert row["uniprot"] == "Q99593"
    assert row["ncbi_gene"] == "6910"
    assert row["locus"] == "12q24.21"
    assert row["aliases"] is None


def _rbm10(uniprot_ids: str) -> str:
    """RBM10's complete-set row, with its `uniprot_ids` cell supplied by the test."""
    return (
        "HGNC:9896\tRBM10\tRNA binding motif protein 10\t\tENSG00000182872\t8241\t"
        f'Xp11.3\t"{uniprot_ids}"\t"ENST00000377604.9|NM_005676.5"\tApproved'
    )


def test_rbm10_mirrors_its_own_protein_rather_than_the_microprotein_hgnc_lists_first(
    tmp_path,
):
    """RBM10 is the one live multi-valued cell: HGNC's `uniprot_ids` reads
    `P0DW28|P98175` -- alphabetical, not ranked -- and no other of the 185
    registry genes lists more than one (measured 2026-10-25 over the 2026-10-16
    complete set). UniProt REST, the same day: P98175 is RBM10_HUMAN, "RNA-binding
    protein 10", 930 aa; P0DW28 is MNS60_HUMAN, 130 aa, the MINAS-60 microprotein
    UniProt describes as "the product of an alternative open reading frame
    (alt-ORF) of transcripts coding for the RBM10 (AC P98175) protein". Both are
    reviewed Swiss-Prot entries and both carry gene name RBM10, so neither field
    could have told a script which one is the gene's own protein.

    Taking the first accession would have published the microprotein as RBM10's
    identity in schema 2.17's bundle, and nothing would have failed: a
    well-formed wrong accession matches `UNIPROT_PATTERN` and breaks no join that
    has data to join. TBX5 rides along so a converter that pinned *every* gene
    to P98175 fails too.
    """
    source = _source(
        tmp_path,
        _rbm10("P0DW28|P98175"),
        "HGNC:11604\tTBX5\tT-box transcription factor 5\t\tENSG00000089225\t6910\t"
        '12q24.21\tQ99593\t"ENST00000405440.7|NM_181486.4"\tApproved',
    )
    out = tmp_path / "genes.tsv"
    convert(source, out, keep={"HGNC:9896", "HGNC:11604"})

    frame = pl.read_csv(out, separator="\t", infer_schema_length=0)
    assert dict(frame.select("hgnc_id", "uniprot").rows()) == {
        "HGNC:11604": "Q99593",
        "HGNC:9896": "P98175",
    }


@pytest.mark.parametrize(
    ("row", "match"),
    [
        pytest.param(
            "HGNC:2488\tNKX2-5\tNK2 homeobox 5\t\tENSG00000183072\t1482\t"
            '5q35.1\t"P52952|Q00000"\t"ENST00000329198.5|NM_004387.4"\tApproved',
            r"HGNC:2488 \(NKX2-5\) lists 2 UniProt accessions \(P52952\|Q00000\)",
            id="a-new-multi-valued-cell",
        ),
        pytest.param(
            _rbm10("P0DW28"),
            r"HGNC:9896 \(RBM10\): the pinned UniProt accession P98175 is not among "
            r"those HGNC lists \(P0DW28\)",
            id="a-stale-override",
        ),
    ],
)
def test_a_uniprot_cell_no_human_has_resolved_is_refused_rather_than_guessed(tmp_path, row, match):
    """The two ways the next RBM10 could publish a wrong protein unseen.

    A second gene gaining a multi-valued cell must reach a curator, never be
    resolved by position -- position is how RBM10 went wrong. And a pinned
    accession HGNC no longer lists for that gene is stale evidence: publishing
    it would assert a mapping the authority has withdrawn. Both refuse before
    anything is written, so a failed refresh leaves no half-mirror behind.
    """
    source = _source(tmp_path, row)
    out = tmp_path / "genes.tsv"
    with pytest.raises(ValueError, match=match):
        convert(source, out, keep={"HGNC:2488", "HGNC:9896"})
    assert not out.exists()


def test_genes_outside_the_keep_set_are_not_mirrored(tmp_path):
    source = _source(
        tmp_path,
        "HGNC:11604\tTBX5\tT-box transcription factor 5\t\tENSG00000089225\t6910\t"
        '12q24.21\tQ99593\t"ENST00000405440.7|NM_181486.4"\tApproved',
        "HGNC:5\tA1BG\talpha-1-B glycoprotein\t\tENSG00000121410\t1\t"
        '19q13.43\tP04217\t"ENST00000263100.8|NM_130786.4"\tApproved',
    )
    out = tmp_path / "genes.tsv"
    assert convert(source, out, keep={"HGNC:11604"}) == 1

    frame = pl.read_csv(out, separator="\t", infer_schema_length=0)
    assert frame["hgnc_id"].to_list() == ["HGNC:11604"]


def test_rows_are_sorted_by_hgnc_id_as_a_string(tmp_path):
    """`TABLE_SCHEMAS["genes"].sort_key` is `("hgnc_id",)` and `sort_order.py`
    compares the published strings, so "HGNC:11604" precedes "HGNC:2488".
    Asserted against a literal rather than by sorting the output again, which
    would only prove the test agrees with itself.
    """
    source = _source(
        tmp_path,
        "HGNC:2488\tNKX2-5\tNK2 homeobox 5\t\tENSG00000183072\t1482\t"
        '5q35.1\tP52952\t"ENST00000329198.5|NM_004387.4"\tApproved',
        "HGNC:11604\tTBX5\tT-box transcription factor 5\t\tENSG00000089225\t6910\t"
        '12q24.21\tQ99593\t"ENST00000405440.7|NM_181486.4"\tApproved',
    )
    out = tmp_path / "genes.tsv"
    convert(source, out, keep={"HGNC:2488", "HGNC:11604"})

    frame = pl.read_csv(out, separator="\t", infer_schema_length=0)
    assert frame["hgnc_id"].to_list() == ["HGNC:11604", "HGNC:2488"]
