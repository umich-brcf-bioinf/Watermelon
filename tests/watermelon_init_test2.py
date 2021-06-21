# cd /nfs/mm-isilon/bioinfcore/ActiveProjects/cgates/Watermelon
# pytest tests/watermelon_init_test2.py 


import pytest

import scripts.watermelon_init as watermelon_init


def test_set_up_refs_validKey():
    ref_dict = watermelon_init._set_up_refs(watermelon_init._DEFAULT_GENOME_REFERENCES, 'GRCh38', 'align_qc')
    assert(sorted(ref_dict.keys()) == ['fastq_screen', 'genome', 'references'])


def test_set_up_refs_invalidKey():
    with pytest.raises(Exception) as e_info:
        watermelon_init._set_up_refs(watermelon_init._DEFAULT_GENOME_REFERENCES, 'foobar', 'align_qc')

