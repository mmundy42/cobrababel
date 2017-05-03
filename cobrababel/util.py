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
