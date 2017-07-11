
def translate(model, reaction_xref_file_name, metabolite_xref_file_name, from_namespace, to_namespace):
    """ Translate IDs in a model from one namespace to another namespace.

    Parameters
    ----------
    model : cobra.core.Model
        Model to translate
    reaction_xref_file_name : str
        Path to cross reference file with ID mapping for reactions
    metabolite_xref_file_name : str
        Path to cross reference file with ID mapping for metabolites
    from_namespace : str
        Namespace of IDs in input model
    to_namespace : str
        Namespace of IDs in output model

    Returns
    -------
    cobra.core.Model
        New model with translated IDs
    """

    return
