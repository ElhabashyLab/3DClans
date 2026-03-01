"""Unit tests for clans3d.utils.fasta_utils."""
import os
import pytest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from clans3d.utils.fasta_utils import (
    extract_uid_from_recordID,
    extract_region_from_record,
    add_region_to_record,
    create_mock_up_record,
    extract_uids_from_fasta,
    extract_records_from_fasta,
    clean_fasta_file,
    copy_records_from_fasta,
)


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
        # underscore-style: everything before first underscore
        ("MYOD1_PIG", "MYOD1"),
    ])
    def test_common_formats(self, record_id, expected):
        assert extract_uid_from_recordID(record_id) == expected

    def test_plain_uid_returned_unchanged(self):
        assert extract_uid_from_recordID("Q9Y5Y9") == "Q9Y5Y9"

    def test_region_only_header(self):
        # header of the form "P12345/10-100" (no pipe, no underscore)
        result = extract_uid_from_recordID("P12345/10-100")
        assert result == "P12345"


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
