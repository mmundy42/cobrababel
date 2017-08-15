import cobrababel
import pytest
from os.path import join
from cobra.io import read_sbml_model


class TestTranslate:
    def test_vmh_translate(self, data_folder, vmh_reaction_xref, vmh_metabolite_xref):
        model_file_name = join(data_folder, 'Btheta.xml')
        model = read_sbml_model(model_file_name)
        original_reactions = set([rxn.id for rxn in model.reactions])
        original_metabolites = set([met.id for met in model.metabolites])
        ms_model = cobrababel.translate(model, vmh_reaction_xref, vmh_metabolite_xref, 'vmh', 'seed')
        vmh_model = cobrababel.translate(ms_model, vmh_reaction_xref, vmh_metabolite_xref, 'seed', 'vmh')
        translate_reactions = set([rxn.id for rxn in vmh_model.reactions])
        translate_metabolites = set([met.id for met in vmh_model.metabolites])
        assert original_reactions == translate_reactions
        assert original_metabolites == translate_metabolites
        return

    def test_bad_reaction_xref(self, data_folder, vmh_metabolite_xref):
        model_file_name = join(data_folder, 'Btheta.xml')
        model = read_sbml_model(model_file_name)
        with pytest.raises(IOError):
            ms_model = cobrababel.translate(model, 'unknown.tsv', vmh_metabolite_xref, 'vmh', 'seed')

    def test_bad_metabolite_xref(self, data_folder, vmh_reaction_xref):
        model_file_name = join(data_folder, 'Btheta.xml')
        model = read_sbml_model(model_file_name)
        with pytest.raises(IOError):
            ms_model = cobrababel.translate(model, vmh_reaction_xref, 'unknown.tsv', 'vmh', 'seed')

    def test_bad_from_namespace(self, data_folder, vmh_reaction_xref, vmh_metabolite_xref):
        model_file_name = join(data_folder, 'Btheta.xml')
        model = read_sbml_model(model_file_name)
        with pytest.raises(ValueError):
            ms_model = cobrababel.translate(model, vmh_reaction_xref, vmh_metabolite_xref, 'bad', 'seed')

    def test_bad_to_namespace(self, data_folder, vmh_reaction_xref, vmh_metabolite_xref):
        model_file_name = join(data_folder, 'Btheta.xml')
        model = read_sbml_model(model_file_name)
        with pytest.raises(ValueError):
            ms_model = cobrababel.translate(model, vmh_reaction_xref, vmh_metabolite_xref, 'vmh', 'bad')