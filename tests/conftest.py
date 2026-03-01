"""
Shared pytest fixtures for Clans-3D unit tests.
"""
import os
import pytest
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")


@pytest.fixture
def fixtures_dir():
    return FIXTURES_DIR


@pytest.fixture
def small_fasta_path():
    return os.path.join(FIXTURES_DIR, "small.fasta")


@pytest.fixture
def small_clans_path():
    return os.path.join(FIXTURES_DIR, "small.clans")


@pytest.fixture
def small_scores_df():
    """A minimal pairwise scores DataFrame with 3 sequences."""
    return pd.DataFrame(
        {
            "Sequence_ID_1": ["P11111", "P11111", "P22222"],
            "Sequence_ID_2": ["P22222", "P33333", "P33333"],
            "score": [1.5e-30, 2.0e-28, 3.5e-25],
        }
    )


@pytest.fixture
def three_seq_records():
    """Three minimal SeqRecord objects matching the small.clans fixture."""
    return [
        SeqRecord(
            Seq("ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLM"),
            id="sp|P11111|PROT1_HUMAN/1-50",
            description="sp|P11111|PROT1_HUMAN/1-50",
        ),
        SeqRecord(
            Seq("WYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIK"),
            id="sp|P22222|PROT2_MOUSE/10-60",
            description="sp|P22222|PROT2_MOUSE/10-60",
        ),
        SeqRecord(
            Seq("MNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYAC"),
            id="sp|P33333|PROT3_RAT",
            description="sp|P33333|PROT3_RAT",
        ),
    ]


@pytest.fixture
def three_coordinates():
    return [
        ("P11111", 0.1, 0.2, 0.3),
        ("P22222", 0.4, 0.5, 0.6),
        ("P33333", 0.7, 0.8, 0.9),
    ]
