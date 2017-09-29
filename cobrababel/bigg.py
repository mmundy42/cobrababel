import requests
import io
from warnings import warn
from time import sleep
import logging

from cobra.core import Model, Metabolite, Reaction, DictList
from cobra.io import load_json_model

# Base URL for BiGG website
bigg_url = 'http://bigg.ucsd.edu/api/v2/'

# Pause after this many requests to BiGG data API
PAUSE_COUNT = 500

# Logger for this module
LOGGER = logging.getLogger(__name__)


def create_bigg_universal_model(validate=False, ignore_pseudo_reactions=True):
    """ Create an universal model from BiGG universal reactions and metabolites.

    Parameters
    ----------
    validate : bool, optional
        When True, perform validity checks on universal COBRA model
    ignore_pseudo_reactions : bool, optional
        When True, do not include pseudo reactions

    Returns
    -------
    cobra.core.Model
        COBRA model object with universal reactions and metabolites
    """

    # Get the current version number for the BiGG database.
    response = requests.get(bigg_url + 'database_version')
    if response.status_code != requests.codes.OK:
        response.raise_for_status()
    version = response.json()

    # Create an empty model.
    universal = Model('bigg_universal', name='BiGG universal model {0}'.format(version['bigg_models_version']))
    universal.notes['last_updated'] = version['last_updated']
    universal.notes['source'] = 'BiGG'

    # Get the list of universal metabolites.
    LOGGER.info('Started download of universal metabolite list')
    response = requests.get(bigg_url + 'universal/metabolites')
    if response.status_code != requests.codes.OK:
        response.raise_for_status()
    metabolite_list = response.json()['results']
    LOGGER.info('Finished download of universal metabolite list')

    # Get the details on each universal metabolite. Note there is no bulk download
    # of universal metabolites so there is separate request to the BiGG server
    # for each metabolite which makes this very slow.
    LOGGER.info('Started download of %d metabolites', len(metabolite_list))
    metabolites = list()
    for index in range(len(metabolite_list)):
        if index % PAUSE_COUNT == 0:
            LOGGER.info('Paused at metabolite number %s', index)
            sleep(1)  # Be nice to BiGG data API and take a breath
        bigg_metabolite = get_bigg_metabolite(metabolite_list[index]['bigg_id'])

        # Add a metabolite for each unique compartment.
        compartment_list = set(x['bigg_id'] for x in bigg_metabolite['compartments_in_models'])
        for compartment in compartment_list:
            bigg_metabolite['compartment_bigg_id'] = compartment
            metabolites.append(bigg_metabolite.copy())
    LOGGER.info('Finished download of metabolites')

    # Add the metabolites to the universal model.
    LOGGER.info('Started adding %d Metabolite objects to universal model', len(metabolites))
    add_bigg_metabolites(metabolites, universal)
    LOGGER.info('Finished adding Metabolite objects to universal model')

    # Get the list of universal reactions.
    LOGGER.info('Started download of universal reaction list')
    response = requests.get(bigg_url + 'universal/reactions')
    if response.status_code != requests.codes.OK:
        response.raise_for_status()
    reaction_list = response.json()['results']
    LOGGER.info('Finished download of universal reaction list')

    # Get the details on each universal reaction. Remember there is no bulk download
    # of universal reactions which makes this very slow.
    LOGGER.info('Started creating Reaction objects for %d reactions', len(reaction_list))
    reactions = list()
    for index in range(len(reaction_list)):
        if index % PAUSE_COUNT == 0:
            LOGGER.info('Paused at reaction number %s', index)
            sleep(1)  # Be nice to BiGG data API and take a breath
        reactions.append(get_bigg_reaction(reaction_list[index]['bigg_id']))
    LOGGER.info('Finished creating %d reaction objects', len(reactions))

    # Add the reactions to the universal model.
    LOGGER.info('Started adding %d Reaction objects to universal model', len(reactions))
    add_bigg_reactions(reactions, universal, ignore_pseudo_reactions)
    LOGGER.info('Finished adding Reaction objects to universal model')

    # If requested, validate the COBRA model.
    if validate:
        warn('Coming soon')

    return universal


def create_cobra_model_from_bigg_model(bigg_id, validate=False):
    """ Create a COBRA model from a BiGG model.

    Parameters
    ----------
    bigg_id : str
        ID of BiGG model
    validate : bool, optional
        When True, perform validity checks on COBRA model

    Returns
    -------
    cobra.core.Model
        COBRA model created from JSON representation of BiGG model
    """

    # Download the JSON representation and details of the model from BiGG.
    LOGGER.info('Started download of %s model', bigg_id)
    response = requests.get('{0}models/{1}'.format(bigg_url, bigg_id))
    if response.status_code != requests.codes.OK:
        response.raise_for_status()
    details = response.json()

    response = requests.get('{0}models/{1}/download'.format(bigg_url, bigg_id))
    if response.status_code != requests.codes.OK:
        response.raise_for_status()
    LOGGER.info('Finished download of %s model', bigg_id)

    # Convert to a cobra.Model object.
    with io.StringIO(response.text) as f:
        model = load_json_model(f)

    # Add some details to the Model object.
    model.name = details['organism']
    model.notes['genome_name'] = details['genome_name']
    model.notes['reference_type'] = details['reference_type']
    model.notes['reference_id'] = details['reference_id']
    model.notes['source'] = 'BiGG'

    # Confirm a few basics.
    if len(model.reactions) != details['reaction_count']:
        warn('{0} reactions in model does not equal {1} in model details'
             .format(len(model.reactions), details['reaction_count']))
    if len(model.metabolites) != details['metabolite_count']:
        warn('{0} metabolites in model does not equal {1} in model details'
             .format(len(model.metabolites), details['metabolite_count']))
    if len(model.genes) != details['gene_count']:
        warn('{0} genes in model does not equal {1} in model details'
             .format(len(model.genes), details['gene_count']))

    # If requested, validate the COBRA model.
    if validate:
        warn('Coming soon')

    return model


def get_bigg_model_list():
    """ Get the list of models available from BiGG website.

    Each entry in the list is a dictionary with keys that include 'bigg_id',
    'organism', 'metabolite_count', 'reaction_count', and 'gene_count'.

    Returns
    -------
    list of dict
        List of models available from BiGG website
    """

    LOGGER.info('Started download of model list')
    response = requests.get(bigg_url + 'models')
    if response.status_code != requests.codes.OK:
        response.raise_for_status()
    LOGGER.info('Finished download of model list')
    return response.json()['results']


def get_bigg_metabolite(bigg_id, model_bigg_id='universal'):
    """ Get a metabolite from the BiGG database.

    Parameters
    ----------
    bigg_id : str
        ID of BiGG metabolite
    model_bigg_id : str, optional
        ID of model containing metabolite

    Returns
    -------
    dict
        Dictionary with metabolite data
    """

    if model_bigg_id == 'universal':
        url = '{0}universal/metabolites/{1}'.format(bigg_url, bigg_id)
    else:
        url = '{0}models/{1}/metabolites/{2}'.format(bigg_url, model_bigg_id, bigg_id)
    LOGGER.debug('Started download of %s metabolite from %s model', bigg_id, model_bigg_id)
    response = requests.get(url)
    if response.status_code != requests.codes.OK:
        response.raise_for_status()
    LOGGER.debug('Finished download of %s metabolite from %s model', bigg_id, model_bigg_id)
    return response.json()


def add_bigg_metabolites(bigg_list, model):
    """ Create a COBRA metabolite from a BiGG metabolite.

    Parameters
    ----------
    bigg_list : list of dict
        List of dictionaries with BiGG metabolite data
    model : cobra.core.Model
        Model to add metabolites to
    """

    # Create a Metabolite object for each BiGG metabolite.
    metabolites = DictList()
    for bigg_metabolite in bigg_list:
        # Available data is different for a metabolite from an organism model versus
        # a metabolite from the universal model.
        if 'compartment_bigg_id' in bigg_metabolite:
            compartment = bigg_metabolite['compartment_bigg_id']
        elif 'compartments_in_models' in bigg_metabolite:
            compartment = bigg_metabolite['compartments_in_models'][0]['bigg_id']
        else:
            raise ValueError('BiGG metabolite {0} does not have a compartment'.format(bigg_metabolite['bigg_id']))
        metabolite = Metabolite(id='{0}_{1}'.format(bigg_metabolite['bigg_id'], compartment),
                                name=bigg_metabolite['name'],
                                compartment=compartment)
        try:
            metabolite.formula = bigg_metabolite['formula']
        except KeyError:
            try:
                if len(bigg_metabolite['formulae']) > 0:
                    metabolite.formula = bigg_metabolite['formulae'][0]
            except KeyError:
                pass
        try:
            metabolite.charge = bigg_metabolite['charge']
        except KeyError:
            try:
                if len(bigg_metabolite['charges']) > 0:
                    metabolite.charge = bigg_metabolite['charges'][0]
            except KeyError:
                pass
        metabolite.notes['aliases'] = bigg_metabolite['database_links']
        metabolites.append(metabolite)

        if compartment not in model.compartments:
            try:
                model.compartments[compartment] = bigg_metabolite['compartment_name']
            except KeyError:
                model.compartments[compartment] = 'unknown'

    # Add all of the metabolites to the model.
    model.add_metabolites(metabolites)
    return


def get_bigg_reaction(bigg_id, model_bigg_id='universal'):
    """ Get a reaction from the BiGG database.

    Parameters
    ----------
    bigg_id : str
        ID of BiGG reaction
    model_bigg_id : str, optional
        ID of model containing reaction

    Returns
    -------
    dict
        Dictionary with reaction data
    """

    if model_bigg_id == 'universal':
        url = '{0}universal/reactions/{1}'.format(bigg_url, bigg_id)
    else:
        url = '{0}models/{1}/reactions/{2}'.format(bigg_url, model_bigg_id, bigg_id)
    LOGGER.debug('Started download of reaction %s from model %s', bigg_id, model_bigg_id)
    response = requests.get(url)
    if response.status_code != requests.codes.OK:
        response.raise_for_status()
    LOGGER.debug('Finished download of reaction %s from model %s', bigg_id, model_bigg_id)
    return response.json()


def add_bigg_reactions(bigg_list, model, ignore_pseudo_reactions=True):
    """ Create a COBRA reaction from a BiGG reaction.

    Parameters
    ----------
    bigg_list : list of dict
        List of dictionaries with BiGG reaction data
    model : cobra.core.Model
        Model to add reactions to
    ignore_pseudo_reactions : bool, optional
        When True, do not include pseudo reactions
    """

    # Create a Reaction object for each BiGG reaction.
    reactions = DictList()
    for bigg_reaction in bigg_list:
        if bigg_reaction['pseudoreaction'] and ignore_pseudo_reactions:
            continue
        reaction = Reaction(id=bigg_reaction['bigg_id'], name=bigg_reaction['name'])
        reaction.notes['aliases'] = bigg_reaction['database_links']
        metabolites = dict()
        for met in bigg_reaction['metabolites']:
            metabolite = model.metabolites.get_by_id('{0}_{1}'.format(met['bigg_id'], met['compartment_bigg_id']))
            metabolites[metabolite] = met['stoichiometry']
        reaction.add_metabolites(metabolites)
        try:
            reaction.bounds = (bigg_reaction['results'][0]['lower_bound'], bigg_reaction['results'][0]['upper_bound'])
        except KeyError:
            if '&#8652' in bigg_reaction['reaction_string']:
                reaction.bounds = (-1000.0, 1000.0)
            else:
                warn('Unknown direction symbol in reaction string {0}'.format(bigg_reaction['reaction_string']))
        reactions.append(reaction)

    # Add all of the reactions to the model.
    model.add_reactions(reactions)
    return


def create_bigg_xref(model, to_namespace, reaction_xref_file_name, metabolite_xref_file_name,
                     reaction_alias_name=None, metabolite_alias_name=None):
    """ Create cross reference files using a model created from BiGG.

    Parameters
    ----------
    model : cobra.core.Model
        COBRA model object created from BiGG database
    to_namespace : str
        Namespace of IDs to cross reference to
    reaction_xref_file_name : str
        Path to cross reference file with ID mapping for reactions
    metabolite_xref_file_name : str
        Path to cross reference file with ID mapping for metabolites
    reaction_alias_name : str, optional
        Name of alias for reactions
    metabolite_alias_name : str, optional
        Name of alias for metabolites
    """

    if model.notes['source'] != 'BiGG':
        warn('Model {0} ({1}) is not a BiGG model'.format(model.id, model.name))
    if reaction_alias_name is None:
        reaction_alias_name = to_namespace
    if metabolite_alias_name is None:
        metabolite_alias_name = to_namespace

    # Build the reaction cross reference from the reaction note for aliases.
    with open(reaction_xref_file_name, 'w') as handle:
        handle.write('bigg\t{0}\n'.format(to_namespace))
        xref_reactions = model.reactions.query(lambda x: reaction_alias_name in x['aliases'], 'notes')
        if len(xref_reactions) == 0:
            raise ValueError('Model {0} ({1}) does not have any reactions with alias name {2}'
                             .format(model.id, model.name, reaction_alias_name))
        for reaction in xref_reactions:
            for alias in reaction.notes['aliases'][reaction_alias_name]:
                handle.write('{0}\t{1}\n'.format(reaction.id, alias['id']))

    # Build the metabolite cross reference from the metabolite note for aliases.
    with open(metabolite_xref_file_name, 'w') as handle:
        handle.write('bigg\t{0}\n'.format(to_namespace))
        xref_metabolites = model.metabolites.query(lambda x: metabolite_alias_name in x['aliases'], 'notes')
        if len(xref_metabolites) == 0:
            raise ValueError('Model {0} ({1}) does not have any metabolites with alias name {2}'
                             .format(model.id, model.name, metabolite_alias_name))
        for metabolite in xref_metabolites:
            for alias in metabolite.notes['aliases'][metabolite_alias_name]:
                handle.write('{0}\t{1}\n'.format(metabolite.id, alias['id']))

    return


def get_bigg_alias_names(model):
    """ Get the set of alias names in a model created from BiGG.

    Parameters
    ----------
    model : cobra.core.Model
        COBRA model object created from BiGG database

    Returns
    -------
    set
        Set of alias names from notes attribute in reactions and metabolites
    """

    names = set()
    for rxn in model.reactions:
        try:
            names.update(set(rxn.notes['aliases'].keys()))
        except KeyError:
            pass
    for met in model.metabolites:
        try:
            names.update(set(met.notes['aliases'].keys()))
        except KeyError:
            pass
    return names
