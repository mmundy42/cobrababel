import requests
from warnings import warn
import re
import logging

from cobra import Model, Metabolite, Reaction, DictList

# Base URL for MetaNetX website
metanetx_url = 'http://www.metanetx.org/cgi-bin/mnxget/mnxref/'

# Supported MetaNetX version (corresponds to MNXref 3.0)
metanetx_version = 'MNXref Version 2017/05/04'

# Regular expression for metabolites in reaction equation
metabolite_re = re.compile(r'(\d*\.\d+|\d+) (MNXM\d+|BIOMASS)@(MNXD[\dX]|BOUNDARY)')

# Logger for this module
LOGGER = logging.getLogger(__name__)


def create_metanetx_universal_model(validate=False, verbose=False):
    """ Create an universal model from MetaNetX universal reactions and metabolites.

    The MetaNetX metabolite list is very large and includes metabolites that are
    not used in any reaction. The returned model only includes metabolites that
    are actually used in a reaction.

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
    metabolite_list = _download_metanetx_file('chem_prop.tsv')

    # Map field names to column numbers (MetaNetX may add fields in future).
    field_names = {
        'MNX_ID': 0,
        'Description': 1,
        'Formula': 2,
        'Charge': 3,
        'Mass': 4,
        'InChI': 5,
        'SMILES': 6,
        'Source': 7,
        'InChIKey': 8
    }

    # Accumulate all available Metabolite objects separately. Later when creating
    # reactions, metabolites are put in a compartment.
    all_metabolites = DictList()

    # Create Metabolite objects for all of the metabolites from the downloaded file.
    # Add all of the universal metabolites from the list.
    LOGGER.info('Started creating Metabolite objects from %d lines in file', len(metabolite_list))
    for index in range(len(metabolite_list)):
        if len(metabolite_list[index]) == 0 or metabolite_list[index][0] == '#':
            continue  # Skip empty lines and comment lines
        fields = metabolite_list[index].split('\t')
        if len(fields) < len(field_names):
            if verbose:
                warn('Skipped metabolite on line {0} with missing fields: {1}'.format(index, metabolite_list[index]))
            continue

        # Create cobra.core.Metabolite from MetaNetX metabolite.
        metabolite = Metabolite(id=fields[field_names['MNX_ID']],
                                name=fields[field_names['Description']],
                                formula=fields[field_names['Formula']])
        charge = fields[field_names['Charge']]
        metabolite.charge = int(charge) if len(charge) > 0 and charge != 'NA' else None
        mass = fields[field_names['Mass']]
        if len(mass) > 0:
            metabolite.notes['mass'] = float(mass)
        metabolite.notes['InChI'] = fields[field_names['InChI']] \
            if len(fields[field_names['InChI']]) > 0 else 'NA'
        metabolite.notes['SMILES'] = fields[field_names['SMILES']] \
            if len(fields[field_names['SMILES']]) > 0 else 'NA'
        metabolite.notes['source'] = fields[field_names['Source']] \
            if len(fields[field_names['Source']]) > 0 else 'NA'
        metabolite.notes['InChIKey'] = fields[field_names['InChIKey']] \
            if len(fields[field_names['InChIKey']]) > 0 else 'NA'
        all_metabolites.append(metabolite)
    LOGGER.info('Finished creating %d Metabolite objects', len(all_metabolites))

    # Download the compartments file.
    compartment_list = _download_metanetx_file('comp_prop.tsv')

    # Map field names to column numbers (MetaNetX may add fields in future).
    field_names = {
        'MNX_ID': 0,
        'Description': 1,
        'Source': 2
    }

    # Add the compartments to the universal model.
    LOGGER.info('Started adding compartments from %d lines in file', len(compartment_list))
    for index in range(len(compartment_list)):
        if len(compartment_list[index]) == 0 or compartment_list[index][0] == '#':
            continue  # Skip empty lines and comment lines
        fields = compartment_list[index].split('\t')
        if len(fields) < len(field_names):
            if verbose:
                warn('Skipped compartment on line {0} with missing fields: {1}'.format(index, compartment_list[index]))
            continue
        universal.compartments[fields[field_names['MNX_ID']]] = fields[field_names['Description']]
    LOGGER.info('Finished adding {0} compartments to universal model'.format(len(universal.compartments)))

    # Download the reactions file.
    reaction_list = _download_metanetx_file('reac_prop.tsv')

    # Map field names to column numbers (hopefully MetaNetX doesn't change this).
    field_names = {
        'MNX_ID': 0,
        'Equation': 1,
        'Description': 2,
        'Balance': 3,
        'EC': 4,
        'Source': 5
    }

    # Accumulate Reaction and Metabolite objects separately because it is faster
    # than adding them one at a time to a model.
    reactions = DictList()
    metabolites = DictList()

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

        # Create cobra.core.Reaction from MetaNetX reaction.
        metabolite_info = _parse_metanetx_equation(fields[field_names['Equation']])
        if metabolite_info is None:
            if verbose:
                warn('Could not parse equation for reaction {0} on line {1}: {2}'
                     .format(fields[field_names['MNX_ID']], index, fields[field_names['Equation']]))
            continue
        rxn_mets = dict()
        for metabolite_id in metabolite_info:
            try:
                rxn_mets[metabolites.get_by_id(metabolite_id)] = metabolite_info[metabolite_id]['coefficient']
            except KeyError:
                metabolite = all_metabolites.get_by_id(metabolite_info[metabolite_id]['mnx_id']).copy()
                metabolite.id = metabolite_id
                metabolite.compartment = metabolite_info[metabolite_id]['compartment']
                metabolites.append(metabolite)
                rxn_mets[metabolite] = metabolite_info[metabolite_id]['coefficient']
        reaction = Reaction(id=fields[field_names['MNX_ID']],
                            name=fields[field_names['MNX_ID']],
                            lower_bound=-1000.0,
                            upper_bound=1000.0)
        reaction.add_metabolites(rxn_mets)
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


def _download_metanetx_file(file_name):
    """ Download and parse a MetaNetX property file.

    Parameters
    ----------
    file_name : str
        Name of property file to download from MetaNetX web site

    Returns
    -------
    list
        List of lines in property file
    """

    LOGGER.info('Started download of %s file', file_name)
    response = requests.get('{0}{1}'.format(metanetx_url, file_name))
    if response.status_code != requests.codes.OK:
        response.raise_for_status()
    LOGGER.info('Finished download of %s file', file_name)
    property_list = response.text.split('\n')
    version = property_list[0].strip('# ')
    if version != metanetx_version:
        warn('MetaNetX version "{0}" in {1} file is not supported version "{2}"'
             .format(version, file_name, metanetx_version))
    return property_list


def _parse_metanetx_equation(equation):
    """ Parse an equation string into a dictionary of metabolite information.

    A MetaNetX equation string has reactants and products separated by an '='
    character and metabolites separate by a '+' character. For example, this
    is the equation for reaction MNXR95362:

    2 MNXM2@MNXD1 + 2 MNXM947@MNXD1 = 1 MNXM4@MNXD1 + 2 MNXM470@MNXD1

    Parameters
    ----------
    equation : str
        Equation string from a MetaNetX reaction

    Returns
    -------
    dict or None
        Dictionary of metabolite information or None if equation string cannot
        be parsed
    """

    # Split the equation string into reactant and product strings.
    parts = equation.split(' = ')
    if len(parts) != 2:
        return None

    # Build a dictionary keyed by metabolite ID with the information needed for
    # setting the metabolites in a cobra.core.Reaction object.
    metabolites = dict()
    reactants = parts[0].split(' + ')
    for r in reactants:
        match = re.search(metabolite_re, r)
        if match is None:
            return None
        met_id = '{0}_{1}'.format(match.group(2), match.group(3))
        metabolites[met_id] = {
            'mnx_id': match.group(2),
            'coefficient': -1.0 * float(match.group(1)),
            'compartment': match.group(3)
        }
    products = parts[1].split(' + ')
    for p in products:
        match = re.search(metabolite_re, p)
        if match is None:
            return None
        met_id = '{0}_{1}'.format(match.group(2), match.group(3))
        metabolites[met_id] = {
            'mnx_id': match.group(2),
            'coefficient': float(match.group(1)),
            'compartment': match.group(3)
        }

    return metabolites
