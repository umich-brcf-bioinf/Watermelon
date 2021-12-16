import pytest

import scripts.split_bam_percentage as split_bam_percentage


def test_parse_lines_basecase():
    lines = ["Total records:                                         41325059",
        "analysis_20211207/alignment_results/05-rseqc/Sample_716-AS-10_ribo_content.in.bam (Alignments consumed by input gene list):1165241",
        "analysis_20211207/alignment_results/05-rseqc/Sample_716-AS-10_ribo_content.ex.bam (Alignments not consumed by input gene list):38329643",
        "analysis_20211207/alignment_results/05-rseqc/Sample_716-AS-10_ribo_content.junk.bam (qcfailed, unmapped reads):1830175"
    ]
    expected_dict = {'Total records': 41325059,
        'in': 1165241,
        'ex': 38329643,
        'junk': 1830175
    }
    actual_dict = split_bam_percentage._parse_lines(lines)
    assert(actual_dict == expected_dict)

def test_parse_lines_emptyfile():
    lines = []
    expected_dict = {}
    actual_dict = split_bam_percentage._parse_lines(lines)
    assert(actual_dict == expected_dict)

def test_calculate_percentages_basecase():
    info = {'Total records': 41325059,
        'in': 1165241,
        'ex': 38329643,
        'junk': 1830175
    }
    expected_tuple = (2.82, 92.75, 4.43)
    actual_tuple = split_bam_percentage._calculate_percentages(info)
    assert(actual_tuple == expected_tuple)

def test_calculate_percentages_zerovals():
    info = {'Total records': 0,
        'in': 0,
        'ex': 0,
        'junk': 0
    }
    # Zero vals - warn and return zeroes
    with pytest.warns(UserWarning, match="Zero value"):
        actual_tuple = split_bam_percentage._calculate_percentages(info)
        assert((0,0,0) == actual_tuple)
