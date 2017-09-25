def format_long_string(string, max_length):
    """ Format a string so it fits in column of a specific width.

    Parameters
    ----------
    string : str
        String to format
    max_length : int
        Maximum length of returned string

    Returns
    -------
    str
        Formatted string
    """

    if len(string) > max_length:
        string = string[:max_length - 3]
        string += '...'
    return string


def get_models_in_folder(source_folder):
    """ Get a list of model file path names based on supported extensions.

    Note that the contents of a file with one of the supported extensions might
    not be a model so this should only be used for folders where all of the files
    are model files.

    Parameters
    ----------
    source_folder : str
        Path to folder with source model files

    Returns
    -------
    list of str
        List of path names to source model files
    """

    source_models = list()
    for filename in listdir(source_folder):
        if filename.endswith('.mat') or filename.endswith('.xml') or \
                filename.endswith('.sbml') or filename.endswith('.json'):
            source_models.append(join(source_folder, filename))
    return source_models
