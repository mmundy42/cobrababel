from six import iteritems
import json
import pandas as pd
from tabulate import tabulate
from cobra.core import DictList
from warnings import warn

from .util import format_long_string

"""
Notes
-----
    For ModelSEED models, I don't understand the boundary reactions rxn13782_c, rxn13783_c,
    and rxn13784_c which I think are sink reactions for protein biosynthesis, DNA replication,
    and RNA transcrption.  The compounds show up as being consumed but are produced by the
    biomass reaction.
"""


def compare_models(model1, model2, details=None, boundary=False):
    """ Compare two models and report differences.

    For a useful comparison, the models must use the same ID types. For example,
    comparing models that use ModelSEED IDs is valid but comparing a model that
    uses ModelSEED IDs with a model that uses BiGG IDs does not work.

    @todo Maybe a parameter with checks, metabolites in reactions, metabolite formulas
        or just do the comparisons
        
    Parameters
    ----------
    model1 : cobra.core.Model
        First model to analyze
    model2 : cobra.core.Model
        Second model to analyze
    details : bool, optional
        When true, print details on differences
    boundary : bool, optional
        When true, print info about boundary reactions
    """

    if details is None:
        details = set()

    reaction_header = ['ID', 'NAME', 'REACTION']
    metabolite_header = ['ID', 'NAME']
    difference_header = ['ID', 'MODEL_1', 'MODEL_2']

    # Compare reactions.
    print('REACTIONS\n' + '---------')
    print('{0} reactions in {1}'.format(len(model1.reactions), model1.id))
    print('{0} reactions in {1}\n'.format(len(model2.reactions), model2.id))

    # See if reactions from first model are in the second model.
    num_matched = 0
    reaction_only_in_one = DictList()
    different_name = DictList()
    different_bounds = DictList()
    different_definition = DictList()
    different_genes = DictList()
    for rxn1 in model1.reactions:
        try:
            rxn2 = model2.reactions.get_by_id(rxn1.id)
            num_matched += 1
            if rxn1.name != rxn2.name:
                different_name.append(rxn1)
            if rxn1.bounds != rxn2.bounds:
                different_bounds.append(rxn1)
            if rxn1.reaction != rxn2.reaction:
                different_definition.append(rxn1)
            if rxn1.gene_reaction_rule != rxn2.gene_reaction_rule:
                different_genes.append(rxn1)
        except KeyError:
            reaction_only_in_one.append(rxn1)
    print('{0} reactions in {1} and {2}'.format(num_matched, model1.id, model2.id))
    print('{0} reactions only in {1}\n'.format(len(reaction_only_in_one), model1.id))

    # If requested, show the details on reactions only in the first model.
    if 'reaction_id' in details and len(reaction_only_in_one) > 0:
        reaction_only_in_one.sort(key=lambda x: x.id)
        output = [[rxn.id, format_long_string(rxn.name, 20), rxn.reaction]
                  for rxn in reaction_only_in_one]
        print(tabulate(output, tablefmt='simple', headers=reaction_header) + '\n')

    # See if reactions from second model are in the first model.
    num_matched = 0
    reaction_only_in_two = DictList()
    for rxn2 in model2.reactions:
        if model1.reactions.has_id(rxn2.id):
            num_matched += 1
        else:
            reaction_only_in_two.append(rxn2)
    print('{0} reactions in both {1} and {2}'.format(num_matched, model2.id, model1.id))
    print('{0} reactions only in {1}\n'.format(len(reaction_only_in_two), model2.id))

    # If requested, show the details on reactions only in the second model.
    if 'reaction_id' in details and len(reaction_only_in_two) > 0:
        reaction_only_in_two.sort(key=lambda x: x.id)
        output = [[rxn.id, format_long_string(rxn.name, 20), rxn.reaction]
                  for rxn in reaction_only_in_two]
        print(tabulate(output, tablefmt='simple', headers=reaction_header) + '\n')

    # Display details on reaction attribute differences.
    print('{0} reactions with different names'.format(len(different_name)))
    if 'reaction_name' in details and len(different_name) > 0:
        different_name.sort(key=lambda x: x.id)
        output = [[rxn.id, rxn.name, model2.reactions.get_by_id(rxn.id).name]
                  for rxn in different_name]
        print(tabulate(output, tablefmt='simple', headers=difference_header) + '\n')
    print('{0} reactions with different bounds'.format(len(different_bounds)))
    if 'reaction_bounds' in details and len(different_bounds) > 0:
        different_bounds.sort(key=lambda x: x.id)
        output = [[rxn.id, rxn.bounds, model2.reactions.get_by_id(rxn.id).bounds]
                  for rxn in different_bounds]
        print(tabulate(output, tablefmt='simple', headers=difference_header) + '\n')
    print('{0} reactions with different definitions'.format(len(different_definition)))
    if 'reaction_definition' in details and len(different_definition) > 0:
        different_definition.sort(key=lambda x: x.id)
        output = [[rxn.id, rxn.reaction, model2.reactions.get_by_id(rxn.id).reaction]
                  for rxn in different_definition]
        print(tabulate(output, tablefmt='simple', headers=difference_header) + '\n')
    print('{0} reactions with different genes'.format(len(different_genes)))
    if 'reaction_gene' in details and len(different_genes) > 0:
        different_genes.sort(key=lambda x: x.id)
        output = [[rxn.id, rxn.gene_reaction_rule, model2.reactions.get_by_id(rxn.id).gene_reaction_rule]
                  for rxn in different_genes]
        print(tabulate(output, tablefmt='simple', headers=difference_header) + '\n')

    # Compare metabolites.
    print('\nMETABOLITES\n' + '-----------')
    print('{0} metabolites in {1}'.format(len(model1.metabolites), model1.id))
    print('{0} metabolites in {1}\n'.format(len(model2.metabolites), model2.id))

    # See if metabolites from first model are in the second model.
    num_matched = 0
    metabolite_only_in_one = DictList()
    different_name = DictList()
    different_formula = DictList()
    different_charge = DictList()
    different_compartment = DictList()
    for m1 in model1.metabolites:
        try:
            m2 = model2.metabolites.get_by_id(m1.id)
            num_matched += 1
            if m1.name != m2.name:
                different_name.append(m1)
            if m1.formula != m2.formula:
                different_formula.append(m1)
            if m1.charge != m2.charge:
                different_charge.append(m1)
            if m1.compartment != m2.compartment:
                different_compartment.append(m1)
        except KeyError:
            metabolite_only_in_one.append(m1)
    print('{0} metabolites in both {1} and {2}'.format(num_matched, model1.id, model2.id))
    print('{0} metabolites only in {1}\n'.format(len(metabolite_only_in_one), model1.id))
    if 'metabolite_id' in details and len(metabolite_only_in_one) > 0:
        metabolite_only_in_one.sort(key=lambda x: x.id)
        output = [[met.id, format_long_string(met.name, 70)] for met in metabolite_only_in_one]
        print(tabulate(output, tablefmt='simple', headers=metabolite_header) + '\n')

    # See if metabolites from second model are in the first model.
    num_matched = 0
    metabolite_only_in_two = DictList()
    for m2 in model2.metabolites:
        if model1.metabolites.has_id(m2.id):
            num_matched += 1
        else:
            metabolite_only_in_two.append(m2)
    print('{0} metabolites in both {1} and {2}'.format(num_matched, model2.id, model1.id))
    print('{0} metabolites only in {1}\n'.format(len(metabolite_only_in_two), model2.id))
    if 'metabolite_id' in details and len(metabolite_only_in_two) > 0:
        metabolite_only_in_two.sort(key=lambda x: x.id)
        output = [[met.id, format_long_string(met.name, 70)] for met in metabolite_only_in_two]
        print(tabulate(output, tablefmt='simple', headers=metabolite_header) + '\n')

    # Display details on metabolite attribute differences.
    print('{0} metabolites with different names'.format(len(different_name)))
    if 'metabolite_name' in details and len(different_name) > 0:
        different_name.sort(key=lambda x: x.id)
        output = [[met.id, met.name, model2.metabolites.get_by_id(met.id).name]
                  for met in different_name]
        print(tabulate(output, tablefmt='simple', headers=difference_header) + '\n')
    print('{0} metabolites with different formulas'.format(len(different_formula)))
    if 'metabolite_formula' in details and len(different_formula) > 0:
        different_formula.sort(key=lambda x: x.id)
        output = [[met.id, met.formula, model2.metabolites.get_by_id(met.id).formula]
                  for met in different_formula]
        print(tabulate(output, tablefmt='simple', headers=difference_header) + '\n')
    print('{0} metabolites with different charges'.format(len(different_charge)))
    if 'metabolite_charge' in details and len(different_charge) > 0:
        different_charge.sort(key=lambda x: x.id)
        output = [[met.id, met.charge, model2.metabolites.get_by_id(met.id).charge]
                  for met in different_charge]
        print(tabulate(output, tablefmt='simple', headers=difference_header) + '\n')
    print('{0} metabolites with different compartments'.format(len(different_compartment)))
    if 'metabolite_compartment' in details and len(different_compartment) > 0:
        different_compartment.sort(key=lambda x: x.id)
        output = [[met.id, met.compartment, model2.metabolites.get_by_id(met.id).compartment]
                  for met in different_compartment]
        print(tabulate(output, tablefmt='simple', headers=difference_header) + '\n')

    # See about system boundary reactions.
    if boundary:
        # Get the list of system boundary reactions from first model.
        model1_boundary = model1.reactions.query(lambda x: x, 'boundary')
        print('{0} reactions are system boundary reactions in {1}'.format(len(model1_boundary), model1.id))

        # Get the list of system boundary reactions from second model.
        model2_boundary = model2.reactions.query(lambda x: x, 'boundary')
        print('{0} reactions are system boundary reactions in {1}'.format(len(model2_boundary), model2.id))

    return
