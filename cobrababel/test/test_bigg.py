import pytest
import cobrababel
from cobra import Model


class TestBigg:
    def test_model_list(self):
        model_list = cobrababel.get_bigg_model_list()
        assert len(model_list) >= 80
        assert 'bigg_id' in model_list[0]
        assert 'organism' in model_list[0]
        assert 'metabolite_count' in model_list[0]
        assert 'reaction_count' in model_list[0]
        assert 'gene_count' in model_list[0]

    def test_create_model(self):
        model = cobrababel.create_cobra_model_from_bigg_model('iAF1260')
        assert len(model.metabolites) == 1668
        assert len(model.reactions) == 2382
        assert len(model.genes) == 1261
        solution = model.optimize()
        assert solution.f == pytest.approx(0.736701)

    def test_add_metabolite(self):
        model = Model('bigg_test')
        cobrababel.add_bigg_metabolites([cobrababel.get_bigg_metabolite('h2o_c', 'iAF1260')], model)
        assert model.metabolites[0].id == 'h2o_c'
        assert model.metabolites[0].name == 'H2O'
        assert model.metabolites[0].compartment == 'c'
        cobrababel.add_bigg_metabolites([cobrababel.get_bigg_metabolite('10fthf')], model)
        assert model.metabolites[1].id == '10fthf_c'
        assert model.metabolites[1].name == '10-Formyltetrahydrofolate'
        assert model.metabolites[1].compartment == 'c'
        assert model.compartments['c'] == 'cytosol'

    def test_add_reaction(self):
        model = Model('bigg_test')
        cobrababel.add_bigg_metabolites([cobrababel.get_bigg_metabolite('cynt_p', 'iAF1260'),
                                         cobrababel.get_bigg_metabolite('h_p', 'iAF1260'),
                                         cobrababel.get_bigg_metabolite('cynt_c', 'iAF1260'),
                                         cobrababel.get_bigg_metabolite('h_c', 'iAF1260')], model)
        cobrababel.add_bigg_reactions([cobrababel.get_bigg_reaction('CYNTt2pp', 'iAF1260')], model)
        assert model.reactions[0].id == 'CYNTt2pp'
        assert model.reactions[0].name == 'Cyanate transport via proton symport (periplasm)'
        assert model.reactions[0].reaction == 'cynt_p + h_p --> cynt_c + h_c'
        assert len(model.reactions[0].metabolites) == 4
        assert len(model.compartments) == 2
