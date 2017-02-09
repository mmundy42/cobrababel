import requests
import io
from warnings import warn
from time import sleep
from cobra import Model, Metabolite, Reaction
from cobra.io import load_json_model

# Base URL for BiGG website
bigg_url = 'http://bigg.ucsd.edu/api/v2/'

# Pause after this many requests to BiGG data API
PAUSE_COUNT = 250


def create_bigg_universal_model(validate=False):
    """ Create an universal model from BiGG universal reactions and metabolites.

    Parameters
    ----------
    validate : bool, optional
        When True, perform validity checks on universal COBRA model

    Returns
    -------
    cobra.Model
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
    response = requests.get(bigg_url + 'universal/metabolites')
    if response.status_code != requests.codes.OK:
        response.raise_for_status()
    metabolite_list = response.json()['results']

    # Keep track of unique compartments across all universal metabolites (hoping every
    # model defines the compartments the same way).
    all_compartments = dict()

    # Get the details on each universal metabolite and add to the model.
    for index in range(len(metabolite_list)):
        if index % PAUSE_COUNT == 0:
            sleep(1)  # Be nice to BiGG data API and take a breath
        bigg_metabolite = get_bigg_metabolite(metabolite_list[index]['bigg_id'])

        # Check for any new compartment IDs included with this metabolite.
        for cindex in range(len(bigg_metabolite['compartments_in_models'])):
            compartment = bigg_metabolite['compartments_in_models'][cindex]['bigg_id']
            if compartment not in all_compartments:
                met = get_bigg_metabolite('{0}_{1}'.format(bigg_metabolite['bigg_id'], compartment),
                                          bigg_metabolite['compartments_in_models'][cindex]['model_bigg_id'])
                all_compartments[compartment] = met['compartment_name']

        # Create a cobra.Metabolite object for each unique compartment.
        for compartment in set(x['bigg_id'] for x in bigg_metabolite['compartments_in_models']):
            bigg_metabolite['compartment_bigg_id'] = compartment
            bigg_metabolite['compartment_name'] = all_compartments[compartment]
            if len(bigg_metabolite['formulae']) > 0:
                bigg_metabolite['formula'] = bigg_metabolite['formulae'][0]
            if len(bigg_metabolite['charges']) > 0:
                bigg_metabolite['charge'] = bigg_metabolite['charges'][0]
            metabolite = add_bigg_metabolite(bigg_metabolite, universal)

            # @todo Still need to figure out how to handle multiple formulae
            if len(bigg_metabolite['formulae']) > 1:
                metabolite.notes['formulae'] = bigg_metabolite['formulae']
            # @todo Still need to decide how to handle multiple charges
            if len(bigg_metabolite['charges']) > 1:
                metabolite.notes['charges'] = bigg_metabolite['charges']

    # Get the list of universal reactions.
    response = requests.get(bigg_url + 'universal/reactions')
    if response.status_code != requests.codes.OK:
        response.raise_for_status()
    reaction_list = response.json()['results']

    # Get the details on each universal reaction and add to the model.
    for index in range(len(reaction_list)):
        if index % PAUSE_COUNT == 0:
            sleep(1)  # Be nice to BiGG data API and take a breath
        add_bigg_reaction(get_bigg_reaction(reaction_list[index]['bigg_id']), universal)

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
    cobra.Model
        COBRA model created from JSON representation of BiGG model
    """

    # Download the JSON representation and details of the model from BiGG.
    response = requests.get('{0}models/{1}'.format(bigg_url, bigg_id))
    if response.status_code != requests.codes.OK:
        response.raise_for_status()
    details = response.json()

    response = requests.get('{0}models/{1}/download'.format(bigg_url, bigg_id))
    if response.status_code != requests.codes.OK:
        response.raise_for_status()

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
    list of dictionary
        List of models available from BiGG website
    """

    response = requests.get(bigg_url + 'models')
    if response.status_code != requests.codes.OK:
        response.raise_for_status()
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
    response = requests.get(url)
    if response.status_code != requests.codes.OK:
        response.raise_for_status()
    return response.json()


def add_bigg_metabolite(bigg_metabolite, model):
    """ Create a COBRA metabolite from a BiGG metabolite and add it to COBRA model.

    Parameters
    ----------
    bigg_metabolite: dict
        Dictionary with BiGG metabolite data
    model: cobra.Model
        Model to add metabolite to

    Returns
    -------
    cobra.Metabolite
        Metabolite object created from BiGG metabolite
    """

    metabolite = Metabolite(id='{0}_{1}'.format(bigg_metabolite['bigg_id'], bigg_metabolite['compartment_bigg_id']),
                            name=bigg_metabolite['name'])
    try:
        metabolite.formula = bigg_metabolite['formula']
        metabolite.charge = bigg_metabolite['charge']
    except KeyError:
        pass
    if len(bigg_metabolite['database_links']) > 0:
        metabolite.notes['aliases'] = bigg_metabolite['database_links']
    model.add_metabolites([metabolite])
    if bigg_metabolite['compartment_bigg_id'] not in model.compartments:
        model.compartments[bigg_metabolite['compartment_bigg_id']] = bigg_metabolite['compartment_name']
    return metabolite


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
    response = requests.get(url)
    if response.status_code != requests.codes.OK:
        response.raise_for_status()
    return response.json()


def add_bigg_reaction(bigg_reaction, model):
    """ Create a COBRA reaction from a BiGG reaction and add it to COBRA model.

    Parameters
    ----------
    bigg_reaction: dict
        Dictionary with BiGG reaction data
    model: cobra.Model
        Model to add reaction to

    Returns
    -------
    cobra.Reaction
        Reaction object created from BiGG reaction
    """

    reaction = Reaction(id=bigg_reaction['bigg_id'], name=bigg_reaction['name'])
    reaction.notes['aliases'] = bigg_reaction['database_links']
    metabolites = dict()
    for m in bigg_reaction['metabolites']:
        metabolite = model.metabolites.get_by_id('{0}_{1}'.format(m['bigg_id'], m['compartment_bigg_id']))
        metabolites[metabolite] = m['stoichiometry']
    reaction.add_metabolites(metabolites)
    try:
        reaction.bounds = (bigg_reaction['results'][0]['lower_bound'], bigg_reaction['results'][0]['upper_bound'])
    except KeyError:
        if '&#8652' in bigg_reaction['reaction_string']:
            reaction.bounds = (-1000.0, 1000.0)
        else:
            warn('Unknown direction symbol in reaction string {0}'.format(bigg_reaction['reaction_string']))
    model.add_reaction(reaction)
    return reaction
