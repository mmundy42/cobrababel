import cobrababel
import pytest
from requests import HTTPError


class TestKEGG:
    def test_list_ids(self):
        id_list = cobrababel.list_kegg_ids('rn')
        assert len(id_list) >= 10550
        assert id_list[0].startswith('rn:R00001')

    def test_list_id_bad_name(self):
        with pytest.raises(HTTPError):
            cobrababel.list_kegg_ids('xx')

    def test_get_reactions(self):
        rxn_list = cobrababel.get_kegg_reactions(['R00001', 'R00013', 'R01175'])
        assert len(rxn_list) == 3
        assert rxn_list[0].id == 'R00001'
        assert rxn_list[1].id == 'R00013'
        assert rxn_list[2].id == 'R01175'

    def test_get_enzymes(self):
        enzyme_list = cobrababel.get_kegg_enzymes(['1.1.1.1', '1.1.1.101', '1.5.3.7'])
        assert len(enzyme_list) == 3
        assert enzyme_list[0].id == '1.1.1.1'
        assert enzyme_list[1].id == '1.1.1.101'
        assert enzyme_list[2].id == '1.5.3.7'

    def test_get_amino_seq(self):
        aa_seq = cobrababel.get_kegg_amino_acid_seq('T00122', ['BT_0001', 'BT_0009', 'BT_0015'])
        assert len(aa_seq) == 22
        assert aa_seq[0].startswith('>bth:BT_0001')
        assert aa_seq[1].startswith('MVSTS')
        assert aa_seq[5].startswith('>bth:BT_0009')
        assert aa_seq[6].startswith('MATIR')
        assert aa_seq[13].startswith('>bth:BT_0015')
        assert aa_seq[14].startswith('MKSIL')

    def test_get_dna_seq(self):
        dna_seq = cobrababel.get_kegg_dna_seq('T00122', ['BT_0002', 'BT_0033', 'BT_0042'])
        assert len(dna_seq) == 45
        assert dna_seq[0].startswith('>bth:BT_0002')
        assert dna_seq[1].startswith('atgatgaa')
        assert dna_seq[18].startswith('>bth:BT_0033')
        assert dna_seq[19].startswith('atgatta')
        assert dna_seq[24].startswith('>bth:BT_0042')
        assert dna_seq[25].startswith('atgaaaa')
