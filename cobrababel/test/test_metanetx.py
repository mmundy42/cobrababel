import cobrababel
from cobra.io import write_sbml_model
from cobra.io.sbml3 import validate_sbml_model
from os.path import join
from os import unlink
import pytest


class TestMetaNetX:
    def test_create_universal(self, test_folder):
        universal = cobrababel.create_metanetx_universal_model()
        assert universal.id == 'metanetx_universal'
        assert len(universal.reactions) >= 42952
        assert len(universal.metabolites) >= 31130
        file_name = join(test_folder, 'metanetx.xml')
        write_sbml_model(universal, file_name)
        model, errors = validate_sbml_model(file_name)
        assert len(errors['other']) == 0
        assert len(errors['SBML errors']) == 0
        assert len(errors['warnings']) == 0
        assert len(errors['validator']) >= 325
        unlink(file_name)

    def test_create_metabolite_xref(self, test_folder):
        file_name = join(test_folder, 'metanetx_metabolite_xref.tsv')
        cobrababel.create_metanetx_metabolite_xref('bigg', file_name)
        assert open(file_name).read().count('\n') >= 5175
        unlink(file_name)

    def test_create_reaction_xref(self, test_folder):
        file_name = join(test_folder, 'metanetx_reaction_xref.tsv')
        cobrababel.create_metanetx_reaction_xref('bigg', file_name)
        assert open(file_name).read().count('\n') >= 15288
        unlink(file_name)

    def test_bad_xref_namespace(self, test_folder):
        file_name = join(test_folder, 'metanetx_reaction_xref.tsv')
        with pytest.raises(ValueError):
            cobrababel.create_metanetx_reaction_xref('foobar', file_name)
