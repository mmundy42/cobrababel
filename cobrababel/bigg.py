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


def create_bigg_universal_model(validate=False):
    """ Create an universal model from BiGG universal reactions and metabolites.

    Parameters
    ----------
    validate : bool, optional
        When True, perform validity checks on universal COBRA model

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
    add_bigg_reactions(reactions, universal)
    LOGGER.info('Finished adding Reaction objects to universal model')

    # If requested, validate the COBRA model.
    if validate:
        warn('Coming soon')

    return universal


def create_cobra_model_from_bigg_model(bigg_id, validate=False):
    """ Create a COBRA model from a BiGG model.

    Parameters
    ----------
    bigg_id: str
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
    bigg_id: str
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
    bigg_list: list of dict
        List of dictionaries with BiGG metabolite data
    model: cobra.core.Model
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
        if len(bigg_metabolite['database_links']) > 0:
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
    bigg_id: str
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


def add_bigg_reactions(bigg_list, model):
    """ Create a COBRA reaction from a BiGG reaction.

    Parameters
    ----------
    bigg_list: list of dict
        List of dictionaries with BiGG reaction data
    model: cobra.core.Model
        Model to add reactions to
    """

    # Create a Reaction object for each BiGG reaction.
    reactions = DictList()
    for bigg_reaction in bigg_list:
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
