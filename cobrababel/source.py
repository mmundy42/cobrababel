from os.path import splitext
from warnings import warn
import re

from cobra.io import load_matlab_model, load_json_model, read_sbml_model, save_json_model
from cobra.core import Model

# Regular expression for compartment suffix on ModelSEED IDs
modelseed_suffix_re = re.compile(r'_([ce])$')

# Regular expression for compartment suffix on BiGG IDs
bigg_suffix_re = re.compile(r'\[([cefghilmnpsruvx])\]$')


def _id_type(object_id):
    """ Figure out the ID type for an object ID.

    Parameters
    ----------
    object_id : str

    Returns
    -------
    str
        ID type ('bigg', 'modelseed', or 'unknown')
    """

    if re.search(bigg_suffix_re, object_id) is not None:
        return 'bigg'
    if re.search(modelseed_suffix_re, object_id) is not None:
        return 'modelseed'
    return 'unknown'


def load_model_from_file(filename):
    """ Load a model from a file based on the extension of the file name.

    Parameters
    ----------
    filename : str
        Path to model file

    Returns
    -------
    cobra.Model
        Model object loaded from file

    Raises
    ------
    Exception
        If model file extension is not supported.
    """

    (root, ext) = splitext(filename)
    if ext == '.mat':
        model = load_matlab_model(filename)
    elif ext == '.xml' or ext == '.sbml':
        model = read_sbml_model(filename)
    elif ext == '.json':
        model = load_json_model(filename)
    else:
        raise IOError('Model file extension not supported for {0}'.format(filename))

    return model


def create_universal_model_from_source(source_models, file_name=None, validate=False):
    """ Create an universal model from a list of source models.

    Parameters
    ----------
    source_models : list of str
        List of path names to source model files
    file_name : str, optional
        Path to file for saving universal COBRA model in JSON format
    validate : bool, optional
        When True, perform validity checks on universal COBRA model

    Returns
    -------
    cobra.Model
        COBRA model object with universal reactions and metabolites
    """

    # Create a new COBRApy Model object.
    universal = Model('universal', name='Universal')
    universal.notes['sources'] = list()

    # Add all of the source models to the universal model.
    for model_filename in source_models:
        # Load the model from a file.
        model = load_model_from_file(model_filename)

        # All metabolites need to have a compartment suffix.
        for metabolite in model.metabolites:
            metabolite.notes['type'] = _id_type(metabolite.id)
        unknown = model.metabolites.query(lambda x: 'unknown' in x['type'], 'notes')
        if len(unknown) > 0:
            raise Exception('Unknown compartment suffixes found in metabolites for {0}'.format(model_filename))

        # Remove any duplicate reactions from model which will leave a model with just the new reactions.
        duplicates = model.reactions.query(lambda x: universal.reactions.has_id(x), 'id')
        model.remove_reactions(duplicates, remove_orphans=True)

        # Remove objectives and genes from model.
        model.objective = {}
        for reaction in model.reactions:
            reaction.gene_reaction_rule = ''

        # Add new metabolites and reactions from the model to the universal model.
        universal += model

        # Add attributes not automatically added above.
        universal.notes['sources'].append({'id': model.id, 'filename': model_filename})
        for compartment in model.compartments:
            if compartment in universal.compartments:
                if model.compartments[compartment] != universal.compartments[compartment]:
                    warn('Model {0} has a different meaning for compartment {1} {2}'
                         .format(model.id, compartment, model.compartments[compartment]))
            else:
                universal.compartments[compartment] = model.compartments[compartment]

    # If requested, validate the COBRA model.
    if validate:
        warn('Coming soon')

    # If requested, save the COBRA model.
    if file_name is not None:
        save_json_model(universal, file_name)

    return universal
