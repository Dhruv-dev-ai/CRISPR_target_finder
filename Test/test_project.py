import os
import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from project import find_crispr_targets, check_Fasta, plot_targets

def test_find_crispr_targets_valid_sequence():
    seq = "ATGCGTACGTACGTACGTACGG"
    df = find_crispr_targets(seq)
    assert not df.empty
    assert "gRNA" in df.columns
    assert df.iloc[0]["PAM"] == "GG"
    assert len(df.iloc[0]["gRNA"]) == 20

def test_find_crispr_targets_no_sites():
    seq = "ATGCGTATATATATATATAAA"
    df = find_crispr_targets(seq)
    assert df.empty


def test_check_fasta_reads(tmp_path):
    # Create temporary FASTA file
    fasta_file = tmp_path / "test.fasta"
    record = SeqRecord(Seq("ATGCGTACGTACGTACGTACGG"), id="test_seq")
    SeqIO.write(record, fasta_file, "fasta")

    seq = check_Fasta(str(fasta_file))
    assert str(seq) == "ATGCGTACGTACGTACGTACGG"

def test_check_fasta_file_not_found():
    with pytest.raises(SystemExit):
        check_Fasta("nonexistent.fasta")


def test_plot_targets_creates_file(tmp_path):
    df = pd.DataFrame({
        "gRNA": ["ATGCGTACGTACGTACGTAC"],
        "PAM": ["GG"],
        "Start": [0],
        "End": [22],
        "GC_Content(%)": [50.0]
    })

    os.chdir(tmp_path)

    plot_targets(df, seq_length=25)

    assert os.path.exists("crispr_targets_plot.png")
