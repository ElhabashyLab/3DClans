"""Integration tests verifying A2M/A3M cleaning produces true UniProt sequences.

These tests fetch actual sequences from UniProt to verify that our cleaning
logic produces the REAL protein sequence, not a truncated version.

Key insight: In A2M/A3M format, lowercase letters are REAL amino acids
(insertions relative to query), not alignment artifacts to be removed.
"""

import pytest
import requests
from pathlib import Path

from clans3d.utils.fasta_utils import clean_aligned_sequence, generate_fasta_from_alignment_file
from Bio import SeqIO


# Real A2M entries from examples/big_fasta_files/C5BN68_TERTT_528-601_b0.3.a2m
# Format: (uid, region_start, region_end, a2m_sequence)
REAL_A2M_ENTRIES = [
    (
        "A0A2M8UWB6",
        524,
        595,
        "..sPQEFASQLAKPLGAKTAQREVEQALRDLHLPFDERRPYALRRLRDRIEANLSGLMGPSVAQDMVESFlpyk"
    ),
    (
        "A0A836ZK00",
        73,
        144,
        "..sPQEFATQLAKPLGAKAAQKEVEQALRDLYLPFDERRPYALRRLRDRIEANLSGLMGPSVAQDMVETFlpyk"
    ),
]


def fetch_uniprot_sequence(uid: str) -> str | None:
    """Fetch full protein sequence from UniProt REST API."""
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200 and response.text.startswith(">"):
            # Parse FASTA and return sequence
            lines = response.text.strip().split("\n")
            return "".join(lines[1:])  # Skip header, join sequence lines
        return None
    except requests.RequestException:
        return None


def extract_region(full_seq: str, start: int, end: int) -> str:
    """Extract region from sequence (1-indexed, inclusive)."""
    return full_seq[start - 1 : end]


@pytest.mark.network
class TestCleanedSequenceMatchesUniProt:
    """Verify cleaned A2M sequences match actual UniProt sequences."""

    @pytest.mark.parametrize("uid,start,end,a2m_seq", REAL_A2M_ENTRIES)
    def test_cleaned_matches_uniprot_region(self, uid, start, end, a2m_seq):
        """Cleaned A2M sequence should exactly match UniProt region."""
        # Fetch real sequence from UniProt
        full_uniprot_seq = fetch_uniprot_sequence(uid)
        
        if full_uniprot_seq is None:
            pytest.skip(f"Could not fetch {uid} from UniProt (network issue or obsolete entry)")
        
        # Extract the region we expect
        expected_region = extract_region(full_uniprot_seq, start, end)
        
        # Clean the A2M sequence
        cleaned = clean_aligned_sequence(a2m_seq)
        
        # They should match exactly
        assert cleaned == expected_region, (
            f"\nUID: {uid}, Region: {start}-{end}\n"
            f"Expected (UniProt): {expected_region}\n"
            f"Got (cleaned A2M):  {cleaned}\n"
            f"Length expected: {len(expected_region)}, got: {len(cleaned)}"
        )

    @pytest.mark.parametrize("uid,start,end,a2m_seq", REAL_A2M_ENTRIES)
    def test_cleaned_length_matches_region_span(self, uid, start, end, a2m_seq):
        """Cleaned sequence length should equal region span."""
        cleaned = clean_aligned_sequence(a2m_seq)
        expected_len = end - start + 1
        
        assert len(cleaned) == expected_len, (
            f"Region {start}-{end} should have {expected_len} residues, got {len(cleaned)}"
        )


@pytest.mark.network
class TestRealA2MFileAgainstUniProt:
    """Test with actual A2M file from examples."""

    def test_example_a2m_file_first_entries(self, tmp_path):
        """Verify first few entries of real A2M file match UniProt."""
        a2m_file = Path(__file__).parent.parent.parent / "examples" / "big_fasta_files" / "C5BN68_TERTT_528-601_b0.3.a2m"
        
        if not a2m_file.exists():
            pytest.skip("Example A2M file not found")
        
        # Parse first 3 entries (after query)
        records = list(SeqIO.parse(str(a2m_file), "fasta"))[:5]
        
        mismatches = []
        for record in records[1:]:  # Skip query sequence
            # Extract UID and region from header like "tr|A0A2M8UWB6|A0A2M8UWB6_PSESP/524-595"
            header_parts = record.id.split("/")
            if len(header_parts) != 2:
                continue
                
            uid_part = header_parts[0]
            region_part = header_parts[1]
            
            # Extract UID (between pipes or the whole thing)
            if "|" in uid_part:
                uid = uid_part.split("|")[1]
            else:
                uid = uid_part
            
            # Parse region
            try:
                start, end = map(int, region_part.split("-"))
            except ValueError:
                continue
            
            # Fetch from UniProt
            full_seq = fetch_uniprot_sequence(uid)
            if full_seq is None:
                continue  # Skip if can't fetch
            
            # Clean and compare
            cleaned = clean_aligned_sequence(str(record.seq))
            expected = extract_region(full_seq, start, end)
            
            if cleaned != expected:
                mismatches.append({
                    "uid": uid,
                    "region": f"{start}-{end}",
                    "expected": expected,
                    "got": cleaned,
                })
        
        assert not mismatches, f"Mismatches found:\n{mismatches}"


class TestLowercaseAreRealResidues:
    """Demonstrate that lowercase letters in A2M are real amino acids."""

    def test_lowercase_s_is_serine(self):
        """The 's' in '..sPQE...' is a real Serine, not an artifact."""
        # In A2M: lowercase = insertion relative to query
        # But it's still a REAL amino acid that exists in UniProt!
        a2m_seq = "..sPQEFAS"
        cleaned = clean_aligned_sequence(a2m_seq)
        
        # The 's' should become 'S' (Serine)
        assert cleaned == "SPQEFAS"
        assert cleaned[0] == "S"  # First residue is Serine

    def test_lowercase_at_end_is_real(self):
        """Lowercase 'lpyk' at end are real residues (Leu-Pro-Tyr-Lys)."""
        a2m_seq = "MVESFlpyk"
        cleaned = clean_aligned_sequence(a2m_seq)
        
        # All 9 residues should be preserved
        assert cleaned == "MVESFLPYK"
        assert len(cleaned) == 9

    def test_removing_lowercase_would_lose_real_residues(self):
        """If we removed lowercase, we'd lose actual protein residues."""
        a2m_seq = "..sPQEFASQLAKPLGAKTAQREVEQALRDLHLPFDERRPYALRRLRDRIEANLSGLMGPSVAQDMVESFlpyk"
        
        # Correct cleaning: keep all letters
        correct_cleaned = clean_aligned_sequence(a2m_seq)
        
        # WRONG: if we removed lowercase (like old buggy implementation)
        wrong_cleaned = ''.join(c for c in a2m_seq if c.isupper())
        
        # The wrong approach loses 5 real residues: s, l, p, y, k
        assert len(correct_cleaned) == 72  # Full region
        assert len(wrong_cleaned) == 67    # Missing 5 residues!
        
        # Region 524-595 must have exactly 72 residues
        expected_len = 595 - 524 + 1
        assert len(correct_cleaned) == expected_len
        assert len(wrong_cleaned) != expected_len  # Bug would cause mismatch
