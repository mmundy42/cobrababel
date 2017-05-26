import requests
from warnings import warn
import re
import logging

from cobra import Model, Metabolite, Reaction, DictList

# Base URL for MetaNetX website
metanetx_url = 'http://www.metanetx.org/cgi-bin/mnxget/mnxref/'

# Regular expression for metabolites in reaction equation
metabolite_re = re.compile(r'(\d*\.\d+|\d+) (MNXM\d+)')

# Logger for this module
LOGGER = logging.getLogger(__name__)


def create_metanetx_universal_model(validate=False, verbose=False):
    """ Create an universal model from MetaNetX universal reactions and metabolites.

    Parameters
    ----------
    validate : bool, optional
        When True, perform validity checks on universal COBRA model
    verbose : bool, optional
        When True, show warning messages

    Returns
    -------
    cobra.Model
        COBRA model object with universal reactions and metabolites
    """

    # Create an empty model.
    universal = Model('metanetx_universal', name='MetaNetX universal model')

    # Download the metabolites file.
    LOGGER.info('Started download of metabolite file')
    response = requests.get('{0}chem_prop.tsv'.format(metanetx_url))
    if response.status_code != requests.codes.OK:
        response.raise_for_status()
    LOGGER.info('Finished download of metabolite file')
    metabolite_list = response.text.split('\n')

    # Map field names to column numbers (hopefully MetaNetX doesn't change this).
    field_names = {
        'MNX_ID': 0,
        'Description': 1,
        'Formula': 2,
        'Charge': 3,
        'Mass': 4,
        'InChI': 5,
        'SMILES': 6,
        'Source': 7
    }

    # Accumulate Metabolite objects separately because it is faster than adding
    # metabolites one at a time to a model.
    metabolites = DictList()

    # Create Metabolite objects for all of the metabolites from the downloaded file.
    # Add all of the universal metabolites from the list.
    LOGGER.info('Started creating Metabolite objects from %d lines in file', len(metabolite_list))
    for index in range(len(metabolite_list)):
        if len(metabolite_list[index]) == 0 or metabolite_list[index][0] == '#':
            continue  # Skip empty lines and comment lines
        fields = metabolite_list[index].split('\t')
        if len(fields) != 8:
            if verbose:
                warn('Skipped metabolite on line {0} with missing fields: {1}'.format(index, metabolite_list[index]))
            continue

        # Create cobra.Metabolite from MetaNetX metabolite.
        metabolite = Metabolite(id=fields[field_names['MNX_ID']],
                                name=fields[field_names['Description']],
                                formula=fields[field_names['Formula']],
                                charge=fields[field_names['Charge']])
        metabolite.notes['mass'] = fields[field_names['Mass']]
        metabolite.notes['InChI'] = fields[field_names['InChI']]
        metabolite.notes['SMILES'] = fields[field_names['SMILES']]
        metabolite.notes['source'] = fields[field_names['Source']]
        metabolites.append(metabolite)
    LOGGER.info('Finished creating %d Metabolite objects', len(metabolites))

    # Add the metabolites to the universal model.
    universal.add_metabolites(metabolites)
    LOGGER.info('Finished adding Metabolite objects to universal model')

    # Download the reactions file.
    LOGGER.info('Started download of reaction file')
    response = requests.get('{0}reac_prop.tsv'.format(metanetx_url))
    if response.status_code != requests.codes.OK:
        response.raise_for_status()
    LOGGER.info('Finished download of reaction file')
    reaction_list = response.text.split('\n')

    # Map field names to column numbers (hopefully MetaNetX doesn't change this).
    field_names = {
        'MNX_ID': 0,
        'Equation': 1,
        'Description': 2,
        'Balance': 3,
        'EC': 4,
        'Source': 5
    }

    # Accumulate Reaction objects separately because it is faster than adding
    # reactions one at a time to a model.
    reactions = DictList()

    # Create Reaction objects for all of the reactions from the downloaded file.
    LOGGER.info('Started creating Reaction objects from %d lines in file', len(reaction_list))
    for index in range(len(reaction_list)):
        if len(reaction_list[index]) == 0 or reaction_list[index][0] == '#':
            continue  # Skip empty lines and comment lines
        fields = reaction_list[index].split('\t')
        if len(fields) != len(field_names):
            if verbose:
                warn('Skipped reaction on line {0} with missing fields: {1}'.format(index, reaction_list[index]))
            continue

        # Create cobra.Reaction from MetaNetX reaction.
        metabolites = _parse_equation(fields[field_names['Equation']], universal)
        if metabolites is None:
            if verbose:
                warn('Could not parse equation for reaction {0} on line {1}: {2}'
                     .format(fields[field_names['MNX_ID']], index, fields[field_names['Equation']]))
            continue
        reaction = Reaction(id=fields[field_names['MNX_ID']],
                            name=fields[field_names['MNX_ID']],
                            lower_bound=-1000.0,
                            upper_bound=1000.0)
        reaction.add_metabolites(metabolites)
        if len(fields[field_names['EC']]) > 0:
            reaction.notes['EC_number'] = fields[field_names['EC']]
        if len(fields[field_names['Source']]) > 1:
            parts = fields[field_names['Source']].split(':')
            if len(parts) == 2:
                reaction.notes['aliases'] = {parts[0]: parts[1]}
            else:
                if verbose:
                    warn('Could not parse source for {0}: {1}'
                         .format(fields[field_names['MNX_ID']], fields[field_names['Source']]))
        reactions.append(reaction)
    LOGGER.info('Finished creating %d Reaction objects', len(reactions))

    # Add the reactions to the universal model.
    universal.add_reactions(reactions)
    LOGGER.info('Finished adding Reaction objects to universal model')

    # If requested, validate the COBRA model.
    if validate:
        warn('Coming soon')

    return universal


def _parse_equation(equation, model):
    """ Parse an equation string into metabolite dictionary for a Reaction object.

        A MetaNetX equation string has reactants and products separated by an '='
        character and metabolites separate by a '+' character. For example:

        2 MNXM2 + 2 MNXM947 = 1 MNXM4 + 2 MNXM470

    Parameters
    ----------
    equation : str
        Equation string from a MetaNetX reaction
    model : cobra.Model
        Model object (metabolites used in reaction must be in model)

    Returns
    -------
    dict or None
        Dictionary of metabolites that can be used as input for Reaction.add_metabolites()
        or None if equation string cannot be parsed
    """

    # Split the equation string into reactant and product strings.
    parts = equation.split(' = ')
    if len(parts) != 2:
        return None

    # Build a dictionary of metabolites that can be used in a cobra.Reaction object.
    metabolites = dict()
    reactants = parts[0].split(' + ')
    for r in reactants:
        match = re.search(metabolite_re, r)
        if match is None:
            return None
        met = model.metabolites.get_by_id(match.group(2))
        metabolites[met] = -1.0 * float(match.group(1))
    products = parts[1].split(' + ')
    for p in products:
        match = re.search(metabolite_re, p)
        if match is None:
            return None
        met = model.metabolites.get_by_id(match.group(2))
        metabolites[met] = float(match.group(1))

    return metabolites
