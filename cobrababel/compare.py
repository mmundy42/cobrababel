from six import iteritems
import json
import pandas as pd
from tabulate import tabulate
from cobra.core import DictList
from warnings import warn
from numpy import isclose

from .util import format_long_string

"""
Notes
-----
    For ModelSEED models, I don't understand the boundary reactions rxn13782_c, rxn13783_c,
    and rxn13784_c which I think are sink reactions for protein biosynthesis, DNA replication,
    and RNA transcrption.  The compounds show up as being consumed but are produced by the
    biomass reaction.
"""

# Header lines for tabulated output
reaction_header = ['ID', 'NAME', 'REACTION']
metabolite_header = ['ID', 'NAME']
gene_header = ['ID', 'NAME']
difference_header = ['ID', 'FIRST', 'SECOND']


def compare_models(model1, model2, details=None, boundary=False):
    """ Compare two models and report differences.

    For a useful comparison, the models must use the same ID types. For example,
    comparing models that use ModelSEED IDs is valid but comparing a model that
    uses ModelSEED IDs with a model that uses BiGG IDs does not work.

    Parameters
    ----------
    model1 : cobra.core.Model
        First model to analyze
    model2 : cobra.core.Model
        Second model to analyze
    details : set, optional
        When specified, print details on given types of differences (see other compare functions)
    boundary : bool, optional
        When true, print info about boundary reactions
    """

    # Compare reactions, metabolites, and genes.
    compare_reactions(model1.reactions, model2.reactions, details=details, id1=model1.id, id2=model2.id)
    compare_metabolites(model1.metabolites, model2.metabolites, details=details, id1=model1.id, id2=model2.id)
    compare_genes(model1.genes, model2.genes, details=details, id1=model1.id, id2=model2.id)

    # See about system boundary reactions.
    if boundary:
        # Get the list of system boundary reactions from first model.
        model1_boundary = model1.reactions.query(lambda x: x, 'boundary')
        print('{0} reactions are system boundary reactions in {1}'.format(len(model1_boundary), model1.id))

        # Get the list of system boundary reactions from second model.
        model2_boundary = model2.reactions.query(lambda x: x, 'boundary')
        print('{0} reactions are system boundary reactions in {1}'.format(len(model2_boundary), model2.id))

    return


def compare_reactions(reaction1, reaction2, details=None, id1='first', id2='second'):
    """ Compare two lists of cobra.core.Reaction objects and report differences.

    To determine if two reactions are the same, the function compares the following 
    attributes: (1) ID {'reaction_id'}, (2) name {'reaction_name'}, (3) bounds
    {'reaction_bounds'}, (4) definition {'reaction_definition'}, (5) gene reaction 
    rule {'reaction_gpr'}. Include the value in {} in the details parameter to
    display the details of reactions where the values are different.
    
    Parameters
    ----------
    reaction1 : cobra.core.DictList
        First list of cobra.core.Reaction objects to analyze
    reaction2 : cobra.core.DictList
        Second list of cobra.core.Reaction objects to analyze
    details : set, optional
        When specified, print details on given types of differences
    id1 : str, optional
        ID for labeling first list of reactions
    id2 : str, optional
        ID for labeling second list of reactions
    """

    if details is None:
        details = set()

    print('REACTIONS\n' + '---------')
    print('{0} reactions in {1}'.format(len(reaction1), id1))
    print('{0} reactions in {1}\n'.format(len(reaction2), id2))

    # See if reactions from first model are in the second model.
    num_matched = 0
    reaction_only_in_one = DictList()
    different_name = DictList()
    different_bounds = DictList()
    different_definition = DictList()
    different_genes = DictList()
    for r1 in reaction1:
        try:
            r2 = reaction2.get_by_id(r1.id)
            num_matched += 1
            if r1.name != r2.name:
                different_name.append(r1)
            if r1.bounds != r2.bounds:
                different_bounds.append(r1)
            if r1.reaction != r2.reaction:
                different = False
                for met, coefficient in iteritems(r1.metabolites):
                    if not isclose(r2.get_coefficient(met.id), coefficient):
                        different = True
                if different:
                    different_definition.append(r1)
            if r1.gene_reaction_rule != r2.gene_reaction_rule:
                different_genes.append(r1)
        except KeyError:
            reaction_only_in_one.append(r1)
    print('{0} reactions in {1} and {2}'.format(num_matched, id1, id2))
    print('{0} reactions only in {1}\n'.format(len(reaction_only_in_one), id1))

    # If requested, show the details on reactions only in the first model.
    if 'reaction_id' in details and len(reaction_only_in_one) > 0:
        reaction_only_in_one.sort(key=lambda x: x.id)
        output = [[rxn.id, format_long_string(rxn.name, 20), rxn.reaction]
                  for rxn in reaction_only_in_one]
        print(tabulate(output, tablefmt='simple', headers=reaction_header) + '\n')

    # See if reactions from second model are in the first model.
    num_matched = 0
    reaction_only_in_two = DictList()
    for r2 in reaction2:
        if reaction1.has_id(r2.id):
            num_matched += 1
        else:
            reaction_only_in_two.append(r2)
    print('{0} reactions in both {1} and {2}'.format(num_matched, id1, id2))
    print('{0} reactions only in {1}\n'.format(len(reaction_only_in_two), id2))

    # If requested, show the details on reactions only in the second model.
    if 'reaction_id' in details and len(reaction_only_in_two) > 0:
        reaction_only_in_two.sort(key=lambda x: x.id)
        output = [[rxn.id, format_long_string(rxn.name, 20), rxn.reaction]
                  for rxn in reaction_only_in_two]
        print('\n' + tabulate(output, tablefmt='simple', headers=reaction_header) + '\n')

    # Display details on reaction attribute differences.
    print('{0} reactions with different names'.format(len(different_name)))
    if 'reaction_name' in details and len(different_name) > 0:
        different_name.sort(key=lambda x: x.id)
        output = [[rxn.id, rxn.name, reaction2.get_by_id(rxn.id).name]
                  for rxn in different_name]
        print('\n' + tabulate(output, tablefmt='simple', headers=difference_header) + '\n')
    print('{0} reactions with different bounds'.format(len(different_bounds)))
    if 'reaction_bounds' in details and len(different_bounds) > 0:
        different_bounds.sort(key=lambda x: x.id)
        output = [[rxn.id, rxn.bounds, reaction2.get_by_id(rxn.id).bounds]
                  for rxn in different_bounds]
        print('\n' + tabulate(output, tablefmt='simple', headers=difference_header) + '\n')
    print('{0} reactions with different definitions'.format(len(different_definition)))
    if 'reaction_definition' in details and len(different_definition) > 0:
        different_definition.sort(key=lambda x: x.id)
        output = [[rxn.id, rxn.reaction, reaction2.get_by_id(rxn.id).reaction]
                  for rxn in different_definition]
        print('\n' + tabulate(output, tablefmt='simple', headers=difference_header) + '\n')
    print('{0} reactions with different genes'.format(len(different_genes)))
    if 'reaction_gpr' in details and len(different_genes) > 0:
        different_genes.sort(key=lambda x: x.id)
        output = [[rxn.id, rxn.gene_reaction_rule, reaction2.get_by_id(rxn.id).gene_reaction_rule]
                  for rxn in different_genes]
        print('\n' + tabulate(output, tablefmt='simple', headers=difference_header) + '\n')

    return


def compare_metabolites(metabolite1, metabolite2, details=None, id1='first', id2='second'):
    """ Compare two lists of cobra.core.Metabolite objects and report differences.

    To determine if two metabolites are the same, the function compares the following 
    attributes: (1) ID {'metabolite_id'}, (2) name {'metabolite_name'}, (3) formula
    {'metabolite_formula'}, (4) charge {'metabolite_charge'}, (5) compartment 
    {'metabolite_compartment'}. Include the value in {} in the details parameter to 
    display the details of metabolites where the values are different.
    
    Parameters
    ----------
    metabolite1 : cobra.core.DictList
        First list of cobra.core.Metabolite objects to analyze
    metabolite2 : cobra.core.DictList
        Second list of cobra.core.Metabolite objects to analyze
    details : set, optional
        When specified, print details on given types of differences
    id1 : str, optional
        ID for labeling first list of metabolites
    id2 : str, optional
        ID for labeling second list of metabolites
    """

    if details is None:
        details = set()

    print('\nMETABOLITES\n' + '-----------')
    print('{0} metabolites in {1}'.format(len(metabolite1), id1))
    print('{0} metabolites in {1}\n'.format(len(metabolite2), id2))

    # See if metabolites from first model are in the second model.
    num_matched = 0
    metabolite_only_in_one = DictList()
    different_name = DictList()
    different_formula = DictList()
    different_charge = DictList()
    different_compartment = DictList()
    for m1 in metabolite1:
        try:
            m2 = metabolite2.get_by_id(m1.id)
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
    print('{0} metabolites in both {1} and {2}'.format(num_matched, id1, id2))
    print('{0} metabolites only in {1}\n'.format(len(metabolite_only_in_one), id1))
    if 'metabolite_id' in details and len(metabolite_only_in_one) > 0:
        metabolite_only_in_one.sort(key=lambda x: x.id)
        output = [[met.id, format_long_string(met.name, 70)] for met in metabolite_only_in_one]
        print('\n' + tabulate(output, tablefmt='simple', headers=metabolite_header) + '\n')

    # See if metabolites from second model are in the first model.
    num_matched = 0
    metabolite_only_in_two = DictList()
    for m2 in metabolite2:
        if metabolite1.has_id(m2.id):
            num_matched += 1
        else:
            metabolite_only_in_two.append(m2)
    print('{0} metabolites in both {1} and {2}'.format(num_matched, id1, id2))
    print('{0} metabolites only in {1}\n'.format(len(metabolite_only_in_two), id2))
    if 'metabolite_id' in details and len(metabolite_only_in_two) > 0:
        metabolite_only_in_two.sort(key=lambda x: x.id)
        output = [[met.id, format_long_string(met.name, 70)] for met in metabolite_only_in_two]
        print('\n' + tabulate(output, tablefmt='simple', headers=metabolite_header) + '\n')

    # Display details on metabolite attribute differences.
    print('{0} metabolites with different names'.format(len(different_name)))
    if 'metabolite_name' in details and len(different_name) > 0:
        different_name.sort(key=lambda x: x.id)
        output = [[met.id, met.name, metabolite2.get_by_id(met.id).name]
                  for met in different_name]
        print('\n' + tabulate(output, tablefmt='simple', headers=difference_header) + '\n')
    print('{0} metabolites with different formulas'.format(len(different_formula)))
    if 'metabolite_formula' in details and len(different_formula) > 0:
        different_formula.sort(key=lambda x: x.id)
        output = [[met.id, met.formula, metabolite2.get_by_id(met.id).formula]
                  for met in different_formula]
        print('\n' + tabulate(output, tablefmt='simple', headers=difference_header) + '\n')
    print('{0} metabolites with different charges'.format(len(different_charge)))
    if 'metabolite_charge' in details and len(different_charge) > 0:
        different_charge.sort(key=lambda x: x.id)
        output = [[met.id, met.charge, metabolite2.get_by_id(met.id).charge]
                  for met in different_charge]
        print('\n' + tabulate(output, tablefmt='simple', headers=difference_header) + '\n')
    print('{0} metabolites with different compartments'.format(len(different_compartment)))
    if 'metabolite_compartment' in details and len(different_compartment) > 0:
        different_compartment.sort(key=lambda x: x.id)
        output = [[met.id, met.compartment, metabolite2.get_by_id(met.id).compartment]
                  for met in different_compartment]
        print('\n' + tabulate(output, tablefmt='simple', headers=difference_header) + '\n')

    return


def compare_genes(gene1, gene2, details=None, id1='first', id2='second'):
    """ Compare two lists of cobra.core.Gene objects and report differences.

    To determine if two genes are the same, the function compares the following 
    attributes: (1) ID {'gene_id'}, (2) name {'gene_name'}. Include the value 
    in {} in the details parameter to display the details of genes where the 
    values are different.
    
    Parameters
    ----------
    gene1 : cobra.core.DictList
        First list of cobra.core.Gene objects to analyze
    gene2 : cobra.core.DictList
        Second list of cobra.core.Gene objects to analyze
    details : set, optional
        When specified, print details on given types of differences
    id1 : str, optional
        ID for labeling first list of genes
    id2 : str, optional
        ID for labeling second list of genes
    """

    if details is None:
        details = set()

    print('\nGENES\n' + '------')
    print('{0} genes in {1}'.format(len(gene1), id1))
    print('{0} genes in {1}\n'.format(len(gene2), id2))

    # See if genes from first list are in the second list.
    num_matched = 0
    gene_only_in_one = DictList()
    different_name = DictList()
    for g1 in gene1:
        try:
            g2 = gene2.get_by_id(g1.id)
            num_matched += 1
            if g1.name.lower() != g2.name.lower():
                different_name.append(g1)
        except KeyError:
            gene_only_in_one.append(g1)
    print('{0} genes in both {1} and {2}'.format(num_matched, id1, id2))
    print('{0} genes only in {1}\n'.format(len(gene_only_in_one), id1))
    if 'gene_id' in details and len(gene_only_in_one) > 0:
        gene_only_in_one.sort(key=lambda x: x.id)
        output = [[gene.id, format_long_string(gene.name, 90)] for gene in gene_only_in_one]
        print('\n' + tabulate(output, tablefmt='simple', headers=gene_header) + '\n')

    # See if genes from second list are in the first list.
    num_matched = 0
    gene_only_in_two = DictList()
    for g2 in gene2:
        if gene1.has_id(g2.id):
            num_matched += 1
        else:
            gene_only_in_two.append(g2)
    print('{0} genes in both {1} and {2}'.format(num_matched, id1, id2))
    print('{0} genes only in {1}\n'.format(len(gene_only_in_two), id2))
    if 'gene_id' in details and len(gene_only_in_two) > 0:
        gene_only_in_two.sort(key=lambda x: x.id)
        output = [[gene.id, format_long_string(gene.name, 90)] for gene in gene_only_in_two]
        print('\n' + tabulate(output, tablefmt='simple', headers=gene_header) + '\n')

    # Display details on gene attribute differences.
    print('{0} genes with different names'.format(len(different_name)))
    if 'gene_name' in details and len(different_name) > 0:
        different_name.sort(key=lambda x: x.id)
        output = [[gene.id, gene.name, gene2.get_by_id(gene.id).name]
                  for gene in different_name]
        print('\n' + tabulate(output, tablefmt='simple', headers=difference_header) + '\n')

    return
