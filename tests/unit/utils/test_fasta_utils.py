"""Unit tests for clans3d.utils.fasta_utils."""
import os
import pytest
import requests
from io import StringIO
from unittest.mock import patch, MagicMock
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from clans3d.utils.fasta_utils import (
    extract_uid_from_recordID,
    extract_region_from_record,
    add_region_to_record,
    create_mock_up_record,
    extract_uids_from_fasta,
    extract_records_from_fasta,
    clean_fasta_file,
    copy_records_from_fasta,
    download_fasta_record,
    remove_non_existing_uniprot_accessions,
    generate_fasta_from_uids_with_regions,
    clean_aligned_sequence,
    generate_fasta_from_alignment_file,
)

MODULE = "clans3d.utils.fasta_utils"


def _fasta_response(uid="P11111", seq="ACDEFGHIKLM"):
    """Mock requests.Response for a successful FASTA reply."""
    resp = MagicMock()
    resp.status_code = 200
    resp.text = f">sp|{uid}|PROT_HUMAN\n{seq}\n"
    return resp


def _not_found_response():
    resp = MagicMock()
    resp.status_code = 404
    resp.text = ""
    return resp


# ---------------------------------------------------------------------------
# extract_uid_from_recordID
# ---------------------------------------------------------------------------

class TestExtractUidFromRecordID:
    @pytest.mark.parametrize("record_id, expected", [
        ("sp|P49811|MYOD1_PIG", "P49811"),
        ("tr|A0A2M8UWB6|A0A2M8UWB6_PSESP", "A0A2M8UWB6"),
        ("sp|P49811|MYOD1_PIG/2-300", "P49811"),
        ("tr|A0A2M8UWB6.1|A0A2M8UWB6_PSESP/524-595", "A0A2M8UWB6"),
        # version suffix stripped
        ("P49811.1", "P49811"),
        ("P49811", "P49811"),
    ])
    def test_common_formats(self, record_id, expected):
        assert extract_uid_from_recordID(record_id) == expected

    def test_plain_uid_returned_unchanged(self):
        assert extract_uid_from_recordID("Q9Y5Y9") == "Q9Y5Y9"

    def test_region_only_header(self):
        # header of the form "P12345/10-100" (no pipe, no underscore)
        result = extract_uid_from_recordID("P12345/10-100")
        assert result == "P12345"

    def test_raises_for_gene_name_style(self):
        # gene-name style IDs (e.g. "MYOD1_PIG") contain no UniProt accession
        with pytest.raises(ValueError):
            extract_uid_from_recordID("MYOD1_PIG")


# ---------------------------------------------------------------------------
# extract_region_from_record
# ---------------------------------------------------------------------------

class TestExtractRegionFromRecord:
    def _make_record(self, header):
        return SeqRecord(Seq("ACDE"), id=header, description="")

    def test_region_parsed_correctly(self):
        record = self._make_record("sp|P49811|MYOD1_PIG/2-300")
        assert extract_region_from_record(record) == (2, 300)

    def test_large_region_numbers(self):
        record = self._make_record("A0A819D2C1_9BILA/1286-1334")
        assert extract_region_from_record(record) == (1286, 1334)

    def test_no_region_returns_none(self):
        record = self._make_record("sp|P49811|MYOD1_PIG")
        assert extract_region_from_record(record) is None

    def test_region_starting_at_one(self):
        record = self._make_record("Y1502_ARCFU/1-68")
        assert extract_region_from_record(record) == (1, 68)


# ---------------------------------------------------------------------------
# add_region_to_record
# ---------------------------------------------------------------------------

class TestAddRegionToRecord:
    def _make_record(self, seq_str="ABCDEFGHIJ", uid="P11111"):
        return SeqRecord(Seq(seq_str), id=uid, description="")

    def test_no_region_leaves_id_unchanged(self):
        record = self._make_record()
        result = add_region_to_record(record, None)
        assert result.id == "P11111"
        assert str(result.seq) == "ABCDEFGHIJ"

    def test_region_appended_to_id(self):
        record = self._make_record("ABCDEFGHIJ")
        result = add_region_to_record(record, (2, 5))
        assert result.id == "P11111/2-5"

    def test_sequence_trimmed_to_region(self):
        record = self._make_record("ABCDEFGHIJ")  # 10 residues
        result = add_region_to_record(record, (3, 7))
        # 1-based inclusive: positions 3-7 = characters at indices 2-6
        assert str(result.seq) == "CDEFG"

    def test_region_out_of_bounds_raises(self):
        record = self._make_record("ABCDE")  # length 5
        with pytest.raises(ValueError):
            add_region_to_record(record, (2, 10))

    def test_start_greater_than_end_raises(self):
        record = self._make_record("ABCDE")
        with pytest.raises(ValueError):
            add_region_to_record(record, (4, 2))

    def test_description_cleared(self):
        record = self._make_record()
        record.description = "old description"
        result = add_region_to_record(record, None)
        assert result.description == ""


# ---------------------------------------------------------------------------
# create_mock_up_record
# ---------------------------------------------------------------------------

class TestCreateMockUpRecord:
    def test_sequence_is_not_found(self):
        record = create_mock_up_record("P11111", None)
        assert str(record.seq) == "not_found"

    def test_id_without_region(self):
        record = create_mock_up_record("P11111", None)
        assert "P11111" in record.id

    def test_id_with_region(self):
        record = create_mock_up_record("P11111", (10, 50))
        assert "P11111" in record.id
        assert "10" in record.id
        assert "50" in record.id


# ---------------------------------------------------------------------------
# extract_uids_from_fasta  (uses the small.fasta fixture on disk)
# ---------------------------------------------------------------------------

class TestExtractUidsFromFasta:
    def test_returns_correct_uids(self, small_fasta_path):
        uids = extract_uids_from_fasta(small_fasta_path)
        assert uids == ["P11111", "P22222", "P33333"]

    def test_returns_list(self, small_fasta_path):
        uids = extract_uids_from_fasta(small_fasta_path)
        assert isinstance(uids, list)


# ---------------------------------------------------------------------------
# extract_records_from_fasta
# ---------------------------------------------------------------------------

class TestExtractRecordsFromFasta:
    def test_returns_three_records(self, small_fasta_path):
        records = extract_records_from_fasta(small_fasta_path)
        assert len(records) == 3

    def test_record_ids_are_uids(self, small_fasta_path):
        records = extract_records_from_fasta(small_fasta_path)
        assert records[0].id == "P11111"
        assert records[1].id == "P22222"
        assert records[2].id == "P33333"


# ---------------------------------------------------------------------------
# clean_fasta_file
# ---------------------------------------------------------------------------

class TestCleanFastaFile:
    def test_keeps_only_listed_ids(self, tmp_path):
        fasta = tmp_path / "plain.fasta"
        fasta.write_text(">P11111\nACDE\n>P22222\nFGHI\n>P33333\nKLMN\n")
        out = str(tmp_path / "cleaned.fasta")
        result_path = clean_fasta_file(str(fasta), ["P11111", "P33333"], out)
        uids = extract_uids_from_fasta(result_path)
        assert "P11111" in uids
        assert "P33333" in uids
        assert "P22222" not in uids

    def test_output_file_created(self, tmp_path):
        fasta = tmp_path / "plain.fasta"
        fasta.write_text(">P11111\nACDE\n")
        out = str(tmp_path / "cleaned.fasta")
        clean_fasta_file(str(fasta), ["P11111"], out)
        assert os.path.exists(out)


# ---------------------------------------------------------------------------
# copy_records_from_fasta
# ---------------------------------------------------------------------------

class TestCopyRecordsFromFasta:
    def test_copies_only_specified_uids(self, small_fasta_path, tmp_path):
        out_path = str(tmp_path / "out.fasta")
        copy_records_from_fasta(small_fasta_path, ["P22222"], out_path)
        uids = extract_uids_from_fasta(out_path)
        assert "P22222" in uids
        assert "P11111" not in uids
        assert "P33333" not in uids

    def test_output_file_is_created(self, small_fasta_path, tmp_path):
        out_path = str(tmp_path / "out.fasta")
        copy_records_from_fasta(small_fasta_path, ["P11111", "P33333"], out_path)
        assert os.path.exists(out_path)


# ---------------------------------------------------------------------------
# download_fasta_record
# ---------------------------------------------------------------------------

class TestDownloadFastaRecord:
    def test_returns_record_on_uniprot_success(self):
        with patch(f"{MODULE}.requests.get", return_value=_fasta_response("P11111")):
            record = download_fasta_record("P11111")
        assert isinstance(record, SeqRecord)

    def test_queries_uniprot_first(self):
        with patch(f"{MODULE}.requests.get", return_value=_fasta_response()) as mock_get:
            download_fasta_record("P11111", upi="UPI000")
        first_url = mock_get.call_args_list[0][0][0]
        assert "uniprotkb" in first_url

    def test_only_uniprot_queried_when_it_succeeds(self):
        with patch(f"{MODULE}.requests.get", return_value=_fasta_response()) as mock_get:
            download_fasta_record("P11111", upi="UPI000")
        assert mock_get.call_count == 1

    def test_falls_back_to_uniparc_on_uniprot_404(self):
        responses = [_not_found_response(), _fasta_response("P11111")]
        with patch(f"{MODULE}.requests.get", side_effect=responses):
            record = download_fasta_record("P11111", upi="UPI000")
        assert isinstance(record, SeqRecord)

    def test_uniparc_url_contains_upi(self):
        responses = [_not_found_response(), _fasta_response()]
        with patch(f"{MODULE}.requests.get", side_effect=responses) as mock_get:
            download_fasta_record("P11111", upi="UPI000ABC")
        second_url = mock_get.call_args_list[1][0][0]
        assert "UPI000ABC" in second_url

    def test_raises_value_error_when_uniprot_404_and_no_upi(self):
        with patch(f"{MODULE}.requests.get", return_value=_not_found_response()):
            with pytest.raises(ValueError, match="not found"):
                download_fasta_record("P11111")

    def test_raises_value_error_when_both_endpoints_404(self):
        with patch(f"{MODULE}.requests.get", return_value=_not_found_response()):
            with pytest.raises(ValueError, match="not found"):
                download_fasta_record("P11111", upi="UPI000")

    def test_network_error_propagates_without_upi(self):
        with patch(f"{MODULE}.requests.get", side_effect=requests.ConnectionError()):
            with pytest.raises(requests.exceptions.RequestException):
                download_fasta_record("P11111")

    def test_region_applied_to_downloaded_record(self):
        seq = "A" * 50
        with patch(f"{MODULE}.requests.get", return_value=_fasta_response("P11111", seq)):
            record = download_fasta_record("P11111", region=(1, 10))
        assert len(record.seq) == 10

    def test_no_region_returns_full_sequence(self):
        seq = "ACDEFGHIKLM"
        with patch(f"{MODULE}.requests.get", return_value=_fasta_response("P11111", seq)):
            record = download_fasta_record("P11111")
        assert len(record.seq) == len(seq)

    def test_response_body_not_starting_with_gt_treated_as_not_found(self):
        # UniProt sometimes returns error messages without FASTA format
        resp = MagicMock()
        resp.status_code = 200
        resp.text = "Error: entry not found"
        with patch(f"{MODULE}.requests.get", return_value=resp):
            with pytest.raises(ValueError):
                download_fasta_record("P11111")


# ---------------------------------------------------------------------------
# remove_non_existing_uniprot_accessions
# ---------------------------------------------------------------------------

class TestRemoveNonExistingUniprotAccessions:
    API_MODULE = "clans3d.utils.fasta_utils.uniprot_accessions_to_uniparc_accessions"

    def _uid_to_upi(self, uids):
        return {uid: f"UPI_{uid}" for uid in uids}

    def test_retains_all_when_all_found(self):
        uids = ["P11111", "P22222"]
        with patch(self.API_MODULE, return_value=self._uid_to_upi(uids)), \
             patch(f"{MODULE}.download_fasta_record", return_value=MagicMock()):
            result = remove_non_existing_uniprot_accessions(uids)
        assert result == ["P11111", "P22222"]

    def test_removes_uid_when_value_error_raised(self):
        uids = ["P11111", "P22222"]
        def download_side_effect(uid, **kwargs):
            if uid == "P22222":
                raise ValueError("not found")
            return MagicMock()
        with patch(self.API_MODULE, return_value=self._uid_to_upi(uids)), \
             patch(f"{MODULE}.download_fasta_record", side_effect=download_side_effect):
            result = remove_non_existing_uniprot_accessions(uids)
        assert "P11111" in result
        assert "P22222" not in result

    def test_removes_uid_when_network_error_raised(self):
        uids = ["P11111", "P22222"]
        def download_side_effect(uid, **kwargs):
            if uid == "P11111":
                raise requests.ConnectionError()
            return MagicMock()
        with patch(self.API_MODULE, return_value=self._uid_to_upi(uids)), \
             patch(f"{MODULE}.download_fasta_record", side_effect=download_side_effect):
            result = remove_non_existing_uniprot_accessions(uids)
        assert "P11111" not in result
        assert "P22222" in result

    def test_returns_empty_list_when_all_fail(self):
        uids = ["P11111", "P22222"]
        with patch(self.API_MODULE, return_value=self._uid_to_upi(uids)), \
             patch(f"{MODULE}.download_fasta_record", side_effect=ValueError("not found")):
            result = remove_non_existing_uniprot_accessions(uids)
        assert result == []

    def test_empty_input_returns_empty_list(self):
        with patch(self.API_MODULE, return_value={}):
            result = remove_non_existing_uniprot_accessions([])
        assert result == []

    def test_upi_is_passed_to_download(self):
        uids = ["P11111"]
        calls = []
        def capture(uid, **kwargs):
            calls.append(kwargs.get("upi"))
            return MagicMock()
        with patch(self.API_MODULE, return_value={"P11111": "UPI_ABC"}), \
             patch(f"{MODULE}.download_fasta_record", side_effect=capture):
            remove_non_existing_uniprot_accessions(uids)
        assert calls == ["UPI_ABC"]


# ---------------------------------------------------------------------------
# generate_fasta_from_uids_with_regions
# ---------------------------------------------------------------------------

class TestGenerateFastaFromUidsWithRegions:
    """Tests for downloading sequences from UniProt (used for TSV input)."""

    API_MODULE = "clans3d.utils.fasta_utils.uniprot_accessions_to_uniparc_accessions"

    def _no_upi(self, uids):
        return {uid: None for uid in uids}

    def test_returns_out_path(self, tmp_path):
        out = str(tmp_path / "out.fasta")
        uids = {"P11111": None}
        with patch(self.API_MODULE, return_value=self._no_upi(uids)), \
             patch(f"{MODULE}.download_fasta_record", return_value=SeqRecord(Seq("ACDE"), id="P11111", description="")):
            result = generate_fasta_from_uids_with_regions(uids, out)
        assert result == out

    def test_writes_file(self, tmp_path):
        out = str(tmp_path / "out.fasta")
        uids = {"P11111": None}
        with patch(self.API_MODULE, return_value=self._no_upi(uids)), \
             patch(f"{MODULE}.download_fasta_record", return_value=SeqRecord(Seq("ACDE"), id="P11111", description="")):
            generate_fasta_from_uids_with_regions(uids, out)
        assert os.path.exists(out)

    def test_record_id_set_to_uid(self, tmp_path):
        out = str(tmp_path / "out.fasta")
        uids = {"P11111": None}
        with patch(self.API_MODULE, return_value=self._no_upi(uids)), \
             patch(f"{MODULE}.download_fasta_record", return_value=SeqRecord(Seq("ACDE"), id="sp|P11111|PROT_HUMAN", description="")):
            generate_fasta_from_uids_with_regions(uids, out)
        records = list(SeqIO.parse(out, "fasta"))
        assert records[0].id == "P11111"

    def test_uses_mock_up_when_download_fails(self, tmp_path):
        out = str(tmp_path / "out.fasta")
        uids = {"P11111": None}
        with patch(self.API_MODULE, return_value=self._no_upi(uids)), \
             patch(f"{MODULE}.download_fasta_record", side_effect=ValueError("not found")):
            generate_fasta_from_uids_with_regions(uids, out)
        records = list(SeqIO.parse(out, "fasta"))
        assert len(records) == 1
        assert str(records[0].seq) == "not_found"

    def test_mock_up_used_on_network_error(self, tmp_path):
        out = str(tmp_path / "out.fasta")
        uids = {"P11111": None}
        with patch(self.API_MODULE, return_value=self._no_upi(uids)), \
             patch(f"{MODULE}.download_fasta_record", side_effect=requests.ConnectionError()):
            generate_fasta_from_uids_with_regions(uids, out)
        records = list(SeqIO.parse(out, "fasta"))
        assert str(records[0].seq) == "not_found"

    def test_partial_failure_writes_both_real_and_mockup(self, tmp_path):
        out = str(tmp_path / "out.fasta")
        uids = {"P11111": None, "P22222": None}
        def download_side_effect(uid, *args, **kwargs):
            if uid == "P22222":
                raise ValueError("not found")
            return SeqRecord(Seq("ACDE"), id=uid, description="")
        with patch(self.API_MODULE, return_value=self._no_upi(uids)), \
             patch(f"{MODULE}.download_fasta_record", side_effect=download_side_effect):
            generate_fasta_from_uids_with_regions(uids, out)
        records = list(SeqIO.parse(out, "fasta"))
        assert len(records) == 2
        seqs = {r.id: str(r.seq) for r in records}
        assert seqs["P11111"] == "ACDE"
        assert seqs[next(k for k in seqs if "P22222" in k)] == "not_found"

    def test_empty_input_writes_empty_file(self, tmp_path):
        out = str(tmp_path / "out.fasta")
        with patch(self.API_MODULE, return_value={}):
            generate_fasta_from_uids_with_regions({}, out)
        records = list(SeqIO.parse(out, "fasta"))
        assert records == []


# ---------------------------------------------------------------------------
# clean_aligned_sequence
# ---------------------------------------------------------------------------

class TestCleanAlignedSequence:
    """Tests for gap removal from A2M/A3M sequences (keeps all real residues)."""

    def test_keeps_lowercase_insertions(self):
        # Lowercase letters are real residues (insertions relative to query)
        assert clean_aligned_sequence("MKLVFFaavvAED") == "MKLVFFAAVVAED"

    def test_removes_dots(self):
        assert clean_aligned_sequence("..MKLVFF") == "MKLVFF"

    def test_removes_dashes(self):
        assert clean_aligned_sequence("-MKLVFF--AED") == "MKLVFFAED"

    def test_uppercases_all_letters(self):
        # All letters (both cases) are kept and uppercased
        assert clean_aligned_sequence("MKlvFF") == "MKLVFF"

    def test_mixed_alignment_characters(self):
        # Real example from A2M file: dots are removed, all letters kept
        assert clean_aligned_sequence("..sPQEFASqlakPLGAK") == "SPQEFASQLAKPLGAK"

    def test_empty_string(self):
        assert clean_aligned_sequence("") == ""

    def test_only_gaps(self):
        # Only dots and dashes - no residues
        assert clean_aligned_sequence("..--..") == ""

    def test_only_lowercase(self):
        # Lowercase letters are real residues, should be kept
        assert clean_aligned_sequence("aabb") == "AABB"

    def test_real_a2m_sequence(self):
        # From examples/big_fasta_files/C5BN68_TERTT_528-601_b0.3.a2m
        # The sequence has leading dots (gaps) and lowercase letters (real insertions)
        seq = "..sPQEFASQLAKPLGAKTAQREVEQALRDLHLPFDERRPYALRRLRDRIEANLSGLMGPSVAQDMVESFlpyk"
        # True sequence: all letters kept (s, lpyk are real residues), dots removed
        expected = "SPQEFASQLAKPLGAKTAQREVEQALRDLHLPFDERRPYALRRLRDRIEANLSGLMGPSVAQDMVESFLPYK"
        assert clean_aligned_sequence(seq) == expected


# ---------------------------------------------------------------------------
# generate_fasta_from_alignment_file
# ---------------------------------------------------------------------------

class TestGenerateFastaFromAlignmentFile:
    """Tests for generating clean FASTA from A2M/A3M alignment files."""

    def test_produces_true_sequences(self, tmp_path):
        """Verify gaps removed but all real residues (including lowercase) kept."""
        # Create a mock A2M file
        a2m_content = """>tr|P11111|P11111_HUMAN/10-20
..sPQEFASQLAK
>tr|P22222|P22222_HUMAN/5-15
MKLVFFaavvAED
"""
        a2m_file = tmp_path / "test.a2m"
        a2m_file.write_text(a2m_content)
        out_path = tmp_path / "out.fasta"

        uids_with_regions = {
            "P11111": (10, 20),
            "P22222": (5, 15),
        }

        generate_fasta_from_alignment_file(uids_with_regions, str(out_path), str(a2m_file))

        records = {r.id: str(r.seq) for r in SeqIO.parse(str(out_path), "fasta")}
        # Lowercase 's' is kept (real residue), dots removed
        assert records["P11111/10-20"] == "SPQEFASQLAK"
        # Lowercase 'aavv' are kept (real insertion residues)
        assert records["P22222/5-15"] == "MKLVFFAAVVAED"

    def test_skips_uids_not_in_mapping(self, tmp_path):
        """Only include UIDs that are in the uids_with_regions mapping."""
        a2m_content = """>tr|P11111|P11111_HUMAN/10-20
MKLVFF
>tr|P22222|P22222_HUMAN/5-15
ACDEFG
"""
        a2m_file = tmp_path / "test.a2m"
        a2m_file.write_text(a2m_content)
        out_path = tmp_path / "out.fasta"

        # Only P11111 is in the mapping
        uids_with_regions = {"P11111": (10, 20)}

        generate_fasta_from_alignment_file(uids_with_regions, str(out_path), str(a2m_file))

        records = list(SeqIO.parse(str(out_path), "fasta"))
        assert len(records) == 1
        assert records[0].id == "P11111/10-20"

    def test_handles_no_region(self, tmp_path):
        """When region is None, output just the UID."""
        a2m_content = """>tr|P11111|P11111_HUMAN
MKLVFF
"""
        a2m_file = tmp_path / "test.a2m"
        a2m_file.write_text(a2m_content)
        out_path = tmp_path / "out.fasta"

        uids_with_regions = {"P11111": None}

        generate_fasta_from_alignment_file(uids_with_regions, str(out_path), str(a2m_file))

        records = list(SeqIO.parse(str(out_path), "fasta"))
        assert len(records) == 1
        assert records[0].id == "P11111"

    def test_skips_invalid_headers(self, tmp_path):
        """Skip entries without valid UniProt accession."""
        a2m_content = """>invalid_header
MKLVFF
>tr|P11111|P11111_HUMAN
ACDEFG
"""
        a2m_file = tmp_path / "test.a2m"
        a2m_file.write_text(a2m_content)
        out_path = tmp_path / "out.fasta"

        uids_with_regions = {"P11111": None}

        generate_fasta_from_alignment_file(uids_with_regions, str(out_path), str(a2m_file))

        records = list(SeqIO.parse(str(out_path), "fasta"))
        assert len(records) == 1
        assert records[0].id == "P11111"

    def test_empty_alignment_file(self, tmp_path):
        """Handle empty alignment file."""
        a2m_file = tmp_path / "empty.a2m"
        a2m_file.write_text("")
        out_path = tmp_path / "out.fasta"

        generate_fasta_from_alignment_file({}, str(out_path), str(a2m_file))

        records = list(SeqIO.parse(str(out_path), "fasta"))
        assert records == []


class TestCleanedSequenceMatchesRegionLength:
    """Verify cleaned sequences have correct length for specified regions."""

    def test_cleaned_length_matches_region_span(self):
        """Cleaned sequence should have exactly (end - start + 1) residues."""
        # Real A2M sequence from C5BN68_TERTT_528-601_b0.3.a2m
        # Region: 524-595 (72 residues expected)
        a2m_seq = "..sPQEFASQLAKPLGAKTAQREVEQALRDLHLPFDERRPYALRRLRDRIEANLSGLMGPSVAQDMVESFlpyk"
        
        cleaned = clean_aligned_sequence(a2m_seq)
        
        # 595 - 524 + 1 = 72 residues
        assert len(cleaned) == 72
        # Should start with S (from lowercase 's') and end with K (from uppercase 'K' or lowercase 'k')
        assert cleaned.startswith("S")
        assert cleaned.endswith("K")

    def test_multiple_regions_correct_length(self):
        """Test multiple A2M entries have correct cleaned lengths."""
        test_cases = [
            # (a2m_sequence, region_start, region_end)
            ("iraPEEIVERLTHRIGETTARAEVNRALEQLGFNVGESRPYALRRLRDEVEANLSGLMGISMATEIMDSElpyk", 528, 601),
            ("..sPQEFASQLAKPLGAKTAQREVEQALRDLHLPFDERRPYALRRLRDRIEANLSGLMGPSVAQDMVESFlpyk", 524, 595),
            ("..sPQEFATQLAKPLGAKAAQKEVEQALRDLYLPFDERRPYALRRLRDRIEANLSGLMGPSVAQDMVETFlpyk", 73, 144),
        ]
        
        for a2m_seq, start, end in test_cases:
            cleaned = clean_aligned_sequence(a2m_seq)
            expected_len = end - start + 1
            assert len(cleaned) == expected_len, (
                f"Region {start}-{end} should have {expected_len} residues, "
                f"got {len(cleaned)}: {cleaned}"
            )

    def test_only_gaps_are_removed(self):
        """Dots and dashes removed, all letters (upper+lower) kept."""
        # Construct a sequence where we know exactly what should remain
        # 5 uppercase + 3 lowercase + 2 dots + 2 dashes = 12 chars input
        # Expected: 8 uppercase letters
        test_seq = "AB..cd--EFgh"
        cleaned = clean_aligned_sequence(test_seq)
        
        assert cleaned == "ABCDEFGH"
        assert len(cleaned) == 8

    def test_real_a2m_file_entry_structure(self, tmp_path):
        """Full pipeline test: A2M file → cleaned FASTA with correct lengths."""
        # Multi-line A2M entry (as in real files)
        a2m_content = """>tr|A0A2M8UWB6|A0A2M8UWB6_PSESP/524-595
..sPQEFASQLAKPLGAKTAQREVEQALRDLHLPFDERRPYALRRLRDRIEANLSGLMGP
SVAQDMVESFlpyk
"""
        a2m_file = tmp_path / "test.a2m"
        a2m_file.write_text(a2m_content)
        out_path = tmp_path / "out.fasta"

        uids_with_regions = {"A0A2M8UWB6": (524, 595)}

        generate_fasta_from_alignment_file(uids_with_regions, str(out_path), str(a2m_file))

        records = list(SeqIO.parse(str(out_path), "fasta"))
        assert len(records) == 1
        
        seq = str(records[0].seq)
        expected_len = 595 - 524 + 1  # 72
        assert len(seq) == expected_len, f"Expected {expected_len} residues, got {len(seq)}"
        
        # Verify all uppercase (no lowercase or gaps leaked through)
        assert seq == seq.upper()
        assert "." not in seq
        assert "-" not in seq
