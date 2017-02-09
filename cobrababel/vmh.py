import requests
from warnings import warn
from tempfile import gettempdir
from os.path import join
from os import unlink
import zipfile
import io

from cobra.io import read_sbml_model, load_matlab_model

# Base URL for Virtual Metabolic Human website
vmh_url = 'https://webdav-r3lab.uni.lu/public/msp/'

# File name of current Recon2 model.
recon2_file_name = 'Recon2.v04.mat'


def create_cobra_model_from_agora_model(agora_name, validate=False):
    """ Create a COBRA model from an AGORA model.

    Parameters
    ----------
    agora_name: str
        Name of AGORA model
    validate : bool, optional
        When True, perform validity checks on COBRA model

    Returns
    -------
    cobra.Model
        COBRA model created from SBML representation of AGORA model
    """

    # Download the SBML file.
    response = requests.get('{0}AGORA/sbml/{1}.xml'.format(vmh_url, agora_name))
    if response.status_code != requests.codes.OK:
        response.raise_for_status()

    # Convert to a cobra.Model object.
    with io.BytesIO(response.content) as f:
        model = read_sbml_model(f)

    # If requested, validate the COBRA model.
    if validate:
        warn('Coming soon')

    return model


def create_cobra_model_from_vmh_recon2(validate=False):
    """ Create a COBRA model from the current Recon2 model.

    Parameters
    ----------
    validate : bool, optional
        When True, perform validity checks on COBRA model

    Returns
    -------
    cobra.Model
        COBRA model created from Matlab representation of Recon2 model
    """

    # Download the zip file that contains the Matlab file and extract it.
    response = requests.get('{0}{1}_.zip'.format(vmh_url, recon2_file_name), stream=True)
    if response.status_code != requests.codes.OK:
        response.raise_for_status()
    zipfile.ZipFile(io.BytesIO(response.content), 'r').extract(recon2_file_name, gettempdir())

    # Convert to a cobra.Model object.
    temp_file = join(gettempdir(), recon2_file_name)
    model = load_matlab_model(temp_file)
    unlink(temp_file)

    # If requested, validate the COBRA model.
    if validate:
        warn('Coming soon')

    return model
