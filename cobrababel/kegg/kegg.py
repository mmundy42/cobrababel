import requests
from warnings import warn

from .KeggReaction import KeggReaction
from .KeggEnzyme import KeggEnzyme

# Base URL for KEGG website
kegg_url = 'http://rest.kegg.jp/'


class QueryError(Exception):
    """ Exception raised when there is an error with a KEGG query """
    pass


def _kegg_request(url, result):
    """ Send a request to the KEGG web service.

    Parameters
    ----------
    url : str
        KEGG API request URL
    result : list of str
        List of lines from query result

    Returns
    -------
    list of str
        Updated list of lines with results from this query
    """

    # Send the request and receive the response from the web service.
    response = requests.get(url)
    if response.status_code != requests.codes.OK:
        response.raise_for_status()

    # When the status code indicates a failure, raise an exception with the details.
    # @todo will this even run?
    if response.status_code != 200:
        if response.status_code == 400:
            reason = 'bad request (syntax error, wrong database name, etc.)'
        elif response.status_code == 404:
            reason = 'not found'
        else:
            reason = 'unknown ({0})'.format(response.status_code)
        raise QueryError('Query {0} failed: {1}'.format(url, reason))

    # Add the response to the current results.
    result.extend(response.text.split('\n')[:-1])  # Remove the empty last line
    return result


def get_kegg_records(db_name, id_list, option=None):
    """ Get records from a KEGG database.

    Parameters
    ----------
    db_name : str
        Name of database to get_kegg_records records from
    id_list : list of str
        List of IDs of records to get_kegg_records
    option : str
        Option value for request

    Returns
    -------
    list of str
         List of lines from records for specified IDs
    """

    # The web service limits the size of the query to 10 items.
    increment = 10
    start = 0
    if len(id_list) > increment:
        end = start + increment
    else:
        end = len(id_list)
    counter = len(id_list)

    # Run through the list of IDs until all have been processed.
    result = list()
    while counter > 0:
        # Build the query string.
        query = ''
        for index in range(start, end):
            query += db_name + ':' + id_list[index] + '+'

        # Adjust for the next query.
        start += increment
        end += increment
        if end > len(id_list):
            end = len(id_list)
        counter -= increment

        # Send the request.
        url = '{0}get/{1}'.format(kegg_url, query[:-1])  # Take off trailing +
        if option is not None:
            url = url + '/' + option
        result = _kegg_request(url, result)

    return result


def list_kegg_ids(db_name):
    """ Return a list of entry identifiers and associated definitions for a given database.

    Parameters
    ----------
    db_name : str
        Name of database to query

    Returns
    -------
    list of str
         List of lines from query result
    """

    # Retrieve the list.
    return _kegg_request('{0}list/{1}'.format(kegg_url, db_name), list())


def get_kegg_metabolites(id_list):

    warn('Getting metabolites (or compounds) from KEGG is not supported yet')
    return


def get_kegg_reactions(id_list):
    """ Get reactions from the KEGG reaction database.

    Parameters
    ----------
    id_list : list of str
        List of KEGG reaction IDs

    Returns
    -------
    list of KeggReaction objects
         Objects for returned reactions
    """

    # Get the records from the reaction database.
    record_list = get_kegg_records('rn', id_list)

    # Convert the returned reaction records to KeggReaction objects.
    reaction_list = list()
    start = 0
    for index in range(len(record_list)):
        if record_list[index] == '///':
            reaction_list.append(KeggReaction(record_list[start:index + 1]))
            start = index + 1

    return reaction_list


def get_kegg_enzymes(id_list):
    """ Get enzymes from the enzyme database.

    Parameters
    ----------
    id_list : list of str
        List of KEGG enzyme IDs

    Returns
    -------
    list of KEGGEnyzme objects
         Objects for returned enzymes
    """

    # Get the records from the enzyme database.
    record_list = get_kegg_records('ec', id_list)

    # Convert the returned enzyme records to KeggEnzyme objects.
    enzyme_list = list()
    start = 0
    for index in range(len(record_list)):
        if record_list[index] == '///':
            enzyme_list.append(KeggEnzyme(record_list[start:index + 1]))
            start = index + 1

    return enzyme_list


def get_kegg_amino_acid_seq(code, gene_list):
    """ Get amino acid sequences for a list of genes from an organism.

    Parameters
    ----------
    code : str
        Organism code (3 or 4 character string)
    gene_list : list of str
        List of KEGG gene IDs in organism

    Returns
    -------
    list of str
        List of lines in FASTA format with amino acid sequence for specified genes
    """

    # @todo Maybe the sequences should be returned in a better data structure
    return get_kegg_records(code, gene_list, option='aaseq')


def get_kegg_dna_seq(code, gene_list):
    """ Get DNA sequences for a list of genes from an organism.

    Parameters
    ----------
    code : str
        Organism code (3 or 4 character string)
    gene_list : list of str
        List of KEGG gene IDs in organism

    Returns
    -------
    list of str
        List of lines in FASTA format with DNA sequence for specified genes
    """

    return get_kegg_records(code, gene_list, option='ntseq')
