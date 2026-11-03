"""Samplesheet pooling ([S132], ``D03`` amendment).

Several single-end rows sharing a ``sample`` value are an explicit
request to pool their files; paired or mixed duplicates still abort.
"""

# COVERAGE: [S132]

import pytest

from parse_samplesheet import SamplesheetError, parse_rows, to_tsv


FIELDS = ["sample", "fastq_1", "fastq_2", "run"]


def _rows(*dicts: dict) -> list[tuple[int, dict]]:
    return [(i, d) for i, d in enumerate(dicts, start=2)]


def test_single_end_rows_sharing_a_sample_are_kept() -> None:
    profile, records = parse_rows(
        FIELDS,
        _rows(
            {"sample": "A", "fastq_1": "/d/A_0.fq", "fastq_2": "",
             "run": "r1"},
            {"sample": "A", "fastq_1": "/d/A_1.fq", "fastq_2": "",
             "run": "r2"},
            {"sample": "B", "fastq_1": "/d/B.fq", "fastq_2": "",
             "run": ""},
        ),
    )
    assert profile == "fastq"
    lines = to_tsv(profile, records).splitlines()
    assert lines[1:] == [
        "A\t/d/A_0.fq\t\tr1",
        "A\t/d/A_1.fq\t\tr2",
        "B\t/d/B.fq\t\t",
    ]


def test_paired_rows_sharing_a_sample_are_rejected() -> None:
    rows = _rows(
        {"sample": "A", "fastq_1": "A1_1.fq", "fastq_2": "A1_2.fq",
         "run": ""},
        {"sample": "A", "fastq_1": "A2_1.fq", "fastq_2": "A2_2.fq",
         "run": ""},
    )
    with pytest.raises(SamplesheetError) as excinfo:
        parse_rows(FIELDS, rows)
    message = str(excinfo.value)
    assert "duplicate sample" in message
    assert "2" in message and "3" in message


def test_paired_and_single_end_rows_sharing_a_sample_are_rejected() -> None:
    rows = _rows(
        {"sample": "A", "fastq_1": "A_1.fq", "fastq_2": "A_2.fq",
         "run": ""},
        {"sample": "A", "fastq_1": "A.fq", "fastq_2": "", "run": ""},
    )
    with pytest.raises(SamplesheetError) as excinfo:
        parse_rows(FIELDS, rows)
    message = str(excinfo.value)
    assert "duplicate sample" in message
    assert "2" in message and "3" in message
