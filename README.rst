CobraBabel: COBRA Model Translator
==================================

Metabolic models are available from many sources and different source systems use different
ID namespaces. That makes combining and comparing models from different sources
difficult. CobraBabel supports the following source systems:

1. `ModelSEED <http://modelseed.org>`_
2. `Virtual Metabolic Human <http://vmh.uni.lu/>`_ (VMH)
3. `Biochemical, Genetic and Genomic knowledge base <http://bigg.ucsd.edu/>`_ (BiGG)
4. `MetaNetX <http://www.metanetx.org/>`_
5. `Kyoto Encyclopedia of Genes and Genomes <http://www.kegg.jp>`_ (KEGG)

Features include:

* Create a universal model from BiGG, MetaNetX, VMH
* Get a list of available models from BiGG
* Create an organism model from BiGG, VMH
* Get a specific reaction or metabolite from BiGG, KEGG

Installation
------------

Use pip to install CobraBabel from
`PyPI <https://pypi.python.org/pypi/cobrababel>`_ (we recommend doing this
inside a `virtual environment
<http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_)::

    pip install cobrababel

Source systems
--------------

CobraBabel uses the data access methods provided by the source systems. The current
version of the source system supported by CobraBabel is shown below. If the source
systems change the interface it is possible for CobraBabel to return an error or
incorrect data. Additional details on the source systems are in the notes below.

Remember that when accessing the source systems performance is variable.

VMH notes
^^^^^^^^^

CobraBabel uses https://vmh.uni.lu/#downloadview

VMH provides a SEED2VMH_translation.csv file which is not the same as the translation table
in Supplementary Table 15 in the AGORA paper. Supplementary Table 15 is provided as an
Excel workbook which was parsed into two files: (1) vmh_metabolite_xref.tsv and
(2) vmh_reaction_xref.tsv which are provided in the "data" folder.

BiGG Notes
^^^^^^^^^^

CobraBabel uses `API version 2 <http://bigg.ucsd.edu/data_access>`_.

Universal metabolites can have more than one formula if different models use a
different formula for a metabolite with the same ID (same situation can occur
with charge). Currently, CobraBabel blindly picks the first formula in the list
of formulae for a universal metabolite.

Universal metabolites have a list of all the compartments where the metabolite
is used. Currently, CobraBabel assumes that all compartment IDs are defined the
same across all models. For example, compartment ID "c" means "cytosol" in all
BiGG models. The first instance of the compartment in a metabolite is the one
that sets the name.

Note that creating a universal model takes a long time.

MetaNetX Notes
^^^^^^^^^^^^^^

CobraBabel uses `MNXref Version 2015/09/03 <http://www.metanetx.org/mnxdoc/mnxref.html>`_.

Some reactions are defined with unspecified stoichiometry coefficients, for example:
`(2n) MNXM1471 + 1 MNXM3341 = (2n) MNXM1 + 1 MNXM4074 + (2n) MNXM537`. Currently,
CobraBabel does not include these reactions when creating an universal model.

There are some reactions where there is a BIOMASS metabolite in the reaction
definition. But the BIOMASS metabolite is not defined so any reactions with the
BIOMASS metabolite are not included in the universal model.

If a reaction has a value in the Source field, there is only one source in the
format: `source:id`. Sources include Rhea (rhea), KEGG (kegg), MetaCyc (metacyc),
UniPathway (upa), The Seed (seed), BiGG (bigg), BioPath (biopath), and Reactome
(reactome). Set the `verbose` parameter when calling `create_metanetx_universal_model()`
to show a warning for reactions with an invalid format in the Source field.

There are no compartments so metabolites are defined without a compartment.

All reactions are defined as bi-directional so the lower and upper bounds are
set to the default values `(-1000.0, 1000.0)` for all reactions.

KEGG Notes
^^^^^^^^^^

There are some limitations on the `KEGG API <http://www.kegg.jp/kegg/rest/>`_,
most importantly:

**Restrictions:** KEGG API is provided for academic use by academic users
belonging to academic institutions. This service should not be used for bulk
data downloads. Please obtain KEGG FTP academic subscription for downloading
KEGG data.

When setting OTU representatives in the organism database, there can be variability
in which organism is selected when there are multiple equally good matches. For
example, in the KEGG database there are two Mycobacterium tuberculosis H37Rv
organisms, one with ID "T00015" and one with ID "T02178". It is random as to which
one gets picked as the OTU representative.

A sample OTU representative source file is provided in the "data/otu-reps.tsv" file
which was retrieved from the KBase Central Data Model using the all_entities_OTU
command.

References
----------

1. `COBRApy: COnstraints-Based Reconstruction and Analysis for Python <http://dx.doi.org/doi:10.1186/1752-0509-7-74>`_
2. `BiGG Models: A platform for integrating, standardizing, and sharing genome-scale models <http://dx.doi.org/doi:10.1093/nar/gkv1049>`_
3. `KEGG: Kyoto Encyclopedia of Genes and Genomes <http://www.kegg.jp>`_
4. `MetaNetX.org: a website and repository for accessing, analyzing, and manipulating metabolic networks <http://dx.doi.org/doi:10.1093/bioinformatics/btt036>`_
5. ``_
6. `AGORA`_
7. `Systems Biology Knowledgebase <http://kbase.us>`_
