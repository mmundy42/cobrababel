from __future__ import absolute_import

import json
import jsonschema
import re
from warnings import warn
from math import log10


def save_visual_json_model(model, file_name, pathways=None, pretty=False, threshold=1E-8):
    """ Save a model to a JSON file for MetabolicNetworks visualization.

    Parameters
    ----------
        model
            Model to analyze (must be optimized)
        file_name : str
            Path to file for saving MetabolicNetworks JSON
        pathways : dict, optional
            Dictionary that identifies what pathways a reaction participates in
        pretty : bool, optional
            When True, pretty print the JSON
        threshold : float, optional
            Tolerance for determining if a flux is zero
    """

    # First check for ModelSEED format compartment suffix. A ModelSEED model uses
    # compartment suffixes in the format "_c0" or "_c".
    compartment_suffix_re = re.compile(r'_([a-z])\d*$')
    match = re.search(compartment_suffix_re, model.reactions[0].id)
    if match is None:
        # Second check for BiGG format compartment suffix. A BiGG model uses compartment
        # suffixes in the format "[c]".
        compartment_suffix_re = re.compile(r'\[([a-z])\]$')
        match = re.search(compartment_suffix_re, model.reactions[0])
        if match is None:
            raise ValueError('Unknown compartment suffix format')
            # @todo Could get_kegg_records around this if compartments are set correctly?

    # Set the ID of the cytosol compartment in the model.
    # @todo What if compartments are not set?
    cytosol = None
    for compartment_id in model.compartments:
        if model.compartments[compartment_id].lower() == 'cytosol':
            cytosol = compartment_id

    # Set compartment on all metabolites if it is not already set.
    for index in range(len(model.metabolites)):
        metabolite = model.metabolites[index]
        if metabolite.compartment is None:
            # @todo What if the metabolite does not have a compartment suffix?
            metabolite.compartment = re.search(compartment_suffix_re, metabolite.id).group(1)

    # Look for boundary reactions to find metabolites that are being consumed and produced.
    consumed_metabolites = set()
    produced_metabolites = set()
    boundary_reactions = model.reactions.query(lambda x: x, 'boundary')
    for reaction in boundary_reactions:
        if len(reaction.reactants) > 0:  # @todo Why do I care about this?
            metabolite = reaction.reactants[0]
        else:
            print('Boundary reaction {0}: {1} has product metabolite and was skipped'.
                  format(reaction.id, reaction.reaction))
            print('  Flux is {0}'.format(reaction.x))
            continue

        # A positive flux means the metabolite is consumed after adjusting for metabolite stoichiometry.
        # Reminder: reaction.x is negative for consumed metabolites.
        flux = reaction.x * reaction.metabolites[metabolite]
        if flux > threshold:
            consumed_metabolites.add(metabolite.id)
        elif flux < threshold:
            produced_metabolites.add(metabolite.id)

    # Create a metabolite structure for each metabolite in the model.
    metabolites = dict()
    for metabolite in model.metabolites:
        metabolite_data = dict()
        metabolite_data['id'] = metabolite.id
        if metabolite.name is not None:
            metabolite_data['name'] = metabolite.name
        else:
            metabolite_data['name'] = metabolite.id
        metabolite_data['name'] = re.sub(compartment_suffix_re, '', metabolite_data['name'])
        metabolite_data['formula'] = metabolite.formula
        metabolite_data['compound'] = re.sub(compartment_suffix_re, '', metabolite.id)
        if metabolite.id in consumed_metabolites:
            metabolite_data['media'] = 1
        else:
            metabolite_data['media'] = 0
        if metabolite.id in produced_metabolites:
            metabolite_data['is_output'] = 1
        else:
            metabolite_data['is_output'] = 0
        metabolite_data['compartment'] = metabolite.compartment + '0'
        metabolites[metabolite.id] = metabolite_data

    # Find maximum upper bound.
    max_upper_bound = max(model.reactions.list_attr("upper_bound"))

    # Create the representation of the model for visualization.
    visual_model = dict()
    name_words = model.name.split(' ')  # Used to generate short name
    visual_model['name'] = {'short': name_words[0][0]+name_words[1][0:2], 'long': model.name}

    # Identify the cofactor metabolites based on the source of the model. The values come
    # from the database used to build the model (e.g. ModelSEED or BiGG).
    if len(model.metabolites.query('cpd00001_c', 'id')) > 0:
        visual_model['cofactors'] = [
            'H2O', 'H+', 'Phosphate', 'ATP', 'ADP',
            'GTP', 'GDP', 'CO2', 'NADH', 'NADHP',
            'Oxidizedferredoxin', 'Reducedferredoxin', 'FAD', 'FADH2'
        ]
    elif len(model.metabolites.query('h2o', 'id')) > 0:
        visual_model['cofactors'] = [
            'h2o', 'h', 'pi', 'atp', 'adp',
            'gtp', 'gdp', 'co2', 'nadh', 'nadph',
            'fad', 'fadh2'
        ]

    # The list of nodes describes all of the reactions in the model.
    visual_model['nodes'] = list()
    use_names = True  # @todo Need to detect if names are available in metabolites
    num_active_reactions = 0
    for reaction in model.reactions:
        # Boundary reactions were handled above so we can skip them here.
        if reaction.boundary:
            continue

        # Make a reaction node.
        reaction_data = dict()
        reaction_data['node_type'] = 'reaction'
        reaction_data['id'] = reaction.id
        reaction_data['reaction'] = re.sub(compartment_suffix_re, '', reaction.id)
        if reaction.name != '':
            reaction_data['name'] = reaction.name
        else:
            reaction_data['name'] = 'unknown'
        reaction_data['definition'] = reaction.build_reaction_string(use_metabolite_names=use_names)
        reaction_data['flux'] = reaction.x
        if 'gapfill_data' in reaction.notes:
            reaction_data['gapfilled'] = 'Yes'
        else:
            reaction_data['gapfilled'] = 'No'
        if cytosol in reaction.get_compartments():
            # Any metabolite in cytosol makes the reaction in cytosol.
            reaction_data['compartment'] = 'c0'
        else:
            # Microbes are simple -- anything else is extracellular.
            reaction_data['compartment'] = 'e0'

        # Reaction direction is set based on lower and upper bounds.
        # @todo Does this need to be adjusted for fva?
        if reaction.lower_bound < 0.0 and reaction.upper_bound > 0.0:
            direction = '='
        elif reaction.lower_bound >= 0.0 and reaction.upper_bound > 0.0:
            direction = '>'
        elif reaction.lower_bound < 0.0 and reaction.upper_bound <= 0.0:
            direction = '<'
        else:
            warn('Could not set direction from lower bound {0} and upper bound {1} for reaction {2}'.
                 format(reaction.lower_bound, reaction.upper_bound, reaction.id))
            continue

        # Features is a list of gene IDs.
        # @todo Remove features because not used in visualization?
        reaction_data['features'] = [gene.id for gene in reaction.genes]

        # Add a metabolite data structure for each metabolite in the reaction.
        reactants = list()
        for metabolite in reaction.reactants:
            metabolite_data = metabolites[metabolite.id]
            metabolite_data['stoich'] = reaction.metabolites[metabolite]
            reactants.append(metabolite_data)
        products = list()
        for metabolite in reaction.products:
            metabolite_data = metabolites[metabolite.id]
            metabolite_data['stoich'] = reaction.metabolites[metabolite]
            products.append(metabolite_data)

        # Sometimes the solver gets a solution close to zero. Round the flux to zero when it is very close.
        if reaction_data['flux'] != 0.0 and abs(reaction_data['flux']) < threshold:
            reaction_data['flux'] = 0.0
        if reaction_data['flux'] != 0.0:
            num_active_reactions += 1

        # For unidirectonal reactions, make sure the flux value is consistent.
        if direction == '>' and reaction_data['flux'] < 0.0:
            warn('Negative flux on forward reaction {0}'.format(reaction_data))
        if direction == '<' and reaction_data['flux'] > 0.0:
            warn('Positive flux on reverse reaction {0}'.format(reaction_data))

        # For reverse reactions, flip the reactants and products and the
        # symbol in the definition.
        if direction == '<':
            reaction_data['reactants'] = products
            reaction_data['products'] = reactants
            parts = reaction_data['definition'].split(' <-- ')
            reaction_data['definition'] = '{0} --> {1}'.format(parts[1], parts[0])

        # For bidirectional reactions, flip the reactants and products when the
        # flux value is negative.
        elif direction == '=':
            if reaction_data['flux'] != 0.0:
                if reaction_data['flux'] < 0.0:
                    reaction_data['reactants'] = products
                    reaction_data['products'] = reactants
                    parts = reaction_data['definition'].split(' <=> ')
                    reaction_data['definition'] = '{0} --> {1}'.format(parts[1], parts[0])
                else:
                    reaction_data['reactants'] = reactants
                    reaction_data['products'] = products
                    reaction_data['definition'] = reaction_data['definition'].replace('<=>', '-->')
            else:
                reaction_data['reactants'] = reactants
                reaction_data['products'] = products

        # For forward reactions there are no changes.
        else:
            reaction_data['reactants'] = reactants
            reaction_data['products'] = products

        # Now the flux value can always be a positive value.
        reaction_data['flux'] = abs(reaction_data['flux'])

        # Scale a non-zero flux value to fit in the range 0-100. There is a limitation in
        # MetabolicNetworks that requires the flux range to be fixed.
        # @todo Pick a different transform based on the max upper bound value?
        if reaction_data['flux'] != 0.0 and max_upper_bound > 100.0:
            # COBRA models have a maximum upper bound for reaction flux of 1000.
            # Flux range is transformed from 0-1000 to 0-99.015.
            reaction_data['flux'] = log10(reaction_data['flux']+1.0) * 33.00

        # Add the pathways that the reaction participates in.
        if pathways is not None:
            try:
                reaction_data['pathway'] = pathways['reactions'][reaction.id]
            except KeyError:
                reaction_data['pathway'] = list()

        # Add the reaction node.
        visual_model['nodes'].append(reaction_data)

    # Add the pathway names if needed.
    if pathways is not None:
        visual_model['pathways'] = pathways['pathways']

    # Create the json file with the visual model.
    if pretty:
        dump_opts = {"indent": 4, "separators": (",", ": "), "sort_keys": True}
    else:
        dump_opts = {}
    jsonschema.validate(visual_model, visual_model_json_schema)  # @todo Move to unit test
    json.dump(visual_model, open(file_name, "w"), **dump_opts)

    return


visual_model_json_schema = {
    "$schema": "http://json-schema.org/draft-04/schema#",
    "title": "Visual COBRA",
    "description": "JSON representation of COBRA model for visualization by MetabolicNetworks",
    "type": "object",
    "properties": {
        "name": {
            "type": "object",
            "properties": {
                "short": {"type": "string"},
                "long": {"type": "string"},
            }
        },
        "cofactors": {
            "type": "array",
            "items": {"type": "string"},
        },
        "nodes": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "node_type": {"type": "string"},
                    "id": {"type": "string"},
                    "reaction": {"type": "string"},  # Is this necessary?
                    "name": {"type": "string"},
                    "definition": {"type": "string"},
                    "flux": {"type": "number"},
                    "gapfilled": {"type": "string", "pattern": "Yes|No"},
                    "compartment": {"type": "string", "pattern": "[a-z]0"},
                    "features": {  # Is this necessary?
                        "type": "array",
                        "items": {"type": "string"},
                    },
                    "pathways": {
                        "type": "array",
                        "items": {"type": "string"},
                    },
                    "reactants": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "properties": {
                                "id": {"type": "string"},
                                "name": {"type": "string"},
                                "compound": {"type": "string"},  # Is this necessary?
                                "media": {"type": "integer"},
                                "compartment": {"type": "string", "pattern": "[a-z]0"},
                                "is_output": {"type": "integer"},
                                "formula": {"type": "string"},
                                "stoich": {"type": "number"},  # Is this necessary?
                            },
                        },
                        "required": ["id", "name", "compound", "media", "compartment"],
                    },
                    "products": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "properties": {
                                "id": {"type": "string"},
                                "name": {"type": "string"},
                                "compound": {"type": "string"},  # Is this necessary?
                                "media": {"type": "integer"},
                                "compartment": {"type": "string", "pattern": "[a-z]0"},
                                "is_output": {"type": "integer"},
                                "formula": {"type": "string"},
                                "stoich": {"type": "number"},  # Is this necessary?
                            },
                        },
                        "required": ["id", "name", "compound", "media", "compartment"],
                    },
                },
                "required": ["id", "name", "reaction", "compartment", "definition", "flux", "reactants", "products"],
            },
        },
        "pathways": {
            "type": "object",
            "properties": {
                "id": {"type": "string"},
                "name": {"type": "string"},
            },
        },
        "version": {
            "type": "integer",
            "default": 1,
        },

    },
    "required": ["name", "nodes"],
    "additionalProperties": False,
}
