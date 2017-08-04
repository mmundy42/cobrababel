CobraBabel: COBRA Model Translator
==================================

Metabolic models are available from many sources and different source systems use different
ID namespaces. That makes combining and comparing models from different sources
difficult. CobraBabel supports the following source systems:

1. `Virtual Metabolic Human <http://vmh.uni.lu/>`_ (VMH)
2. `Biochemical, Genetic and Genomic knowledge base <http://bigg.ucsd.edu/>`_ (BiGG)
3. `MetaNetX <http://www.metanetx.org/>`_
4. `Kyoto Encyclopedia of Genes and Genomes <http://www.kegg.jp>`_ (KEGG)

Features include:

* Create a universal model from VMH, BiGG, MetaNetX
* Get a list of available models from BiGG
* Create an organism model from  VMH, BiGG
* Get a specific reaction or metabolite from BiGG, KEGG
* Get a specific enzyme, DNA sequence, or amino acid sequence from KEGG

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
version of the source system supported by CobraBabel is shown below. If a source
system changes its interface, it is possible for CobraBabel to return an error or
incorrect data. Additional details on the source systems are in the notes below.

CobraBabel uses web services provided by other organizations which can be offline,
the interface can change, or the URL can change. CobraBabel uses these default URLs:

* VMH download service at https://webdav-r3lab.uni.lu/public/msp
* BiGG web service at http://bigg.ucsd.edu/api/v2
* MetaNetX download service at http://www.metanetx.org/cgi-bin/mnxget/mnxref
* KEGG web service at http://rest.kegg.jp

VMH notes
^^^^^^^^^

CobraBabel uses the most recent versions of `Recon2 and AGORA
<https://vmh.uni.lu/#downloadview>`_.

VMH provides a SEED2VMH_translation.csv file which is not the same as the translation table
in Supplementary Table 15 in the AGORA paper. Supplementary Table 15 is provided as an
Excel workbook which was parsed into two files: (1) vmh_metabolite_xref.tsv and
(2) vmh_reaction_xref.tsv which are provided in the "data" folder.

BiGG Notes
^^^^^^^^^^

CobraBabel uses `BiGG API version 2 <http://bigg.ucsd.edu/data_access>`_.

Universal metabolites can have more than one formula if different models use a
different formula for a metabolite with the same ID (same situation can occur
with charge). Currently, CobraBabel blindly picks the first formula in the list
of formulae for a universal metabolite.

Universal metabolites have a list of all the compartments where the metabolite
is used but a compartment name is not provided. Since different models may use
different names for the same compartment ID, when a universal metabolite is
added to a model the compartment name is not set.

There is no bulk download of universal metabolites and reactions so creating a
universal model is very slow because each metabolite and reaction is downloaded
separately.

MetaNetX Notes
^^^^^^^^^^^^^^

CobraBabel uses `MNXref Version 2017/05/04 <http://www.metanetx.org/mnxdoc/mnxref.html>`_.

The list of metabolites is very large and includes metabolites that are not used
in any reaction. CobraBabel only includes metabolites that are used in a reaction
when creating a universal model.

Some reactions are defined with unspecified stoichiometry coefficients. For example,
reaction MNXR109485 has the equation:
``(2n) MNXM17@MNXD1 + 1 MNXM18575@MNXD1 = 1 MNXM1@MNXD1 + (n) MNXM47@MNXD1 + (n) MNXM87@MNXD1``.
Currently, CobraBabel does not include these reactions when creating a universal model.

A reaction equation is defined with metabolites in a generic compartment using
the syntax ``MNXM17@MNXD1`` to identify the metabolite and compartment. CobraBabel
replaces the "@" character with an "_" character in the metabolite ID when creating
a universal model.

If a reaction has a value in the Source field, there is only one source in the
format: ``source:id``. Sources include Rhea (rhea), KEGG (kegg), MetaCyc (metacyc),
UniPathway (upa), The Seed (seed), BiGG (bigg), BioPath (biopath), and Reactome
(reactome). Set the ``verbose`` parameter when calling ``create_metanetx_universal_model()``
to show a warning for reactions with an invalid format in the Source field.

There are no names for reactions so the reaction name is set to the reaction ID.

All reactions are defined as bi-directional so the lower and upper bounds are
set to the default values ``(-1000.0, 1000.0)`` for all reactions.

KEGG Notes
^^^^^^^^^^

CobraBabel uses the current version of the `KEGG API <http://www.kegg.jp/kegg/rest/>`_.

Note, there are some limitations on using the KEGG API, most importantly:

    **Restrictions:** KEGG API is provided for academic use by academic users
    belonging to academic institutions. This service should not be used for bulk
    data downloads. Please obtain KEGG FTP academic subscription for downloading
    KEGG data.

    -- www.kegg.jp

When setting OTU representatives in the organism database, there can be variability
in which organism is selected when there are multiple equally good matches. For
example, in the KEGG database there are two Mycobacterium tuberculosis H37Rv
organisms, one with ID "T00015" and one with ID "T02178". It is random as to which
one gets picked as the OTU representative.

A sample OTU representative source file is provided in the "data/otu-reps.tsv" file
which was retrieved from the KBase Central Data Model using the all_entities_OTU
command.

Release Notes
-------------

Version 0.1.3 (August 4, 2017)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Updated create_metanetx_universal_model() to support MetaNetX/MNXref Version 3.0
  including new InChIKey value for metabolites and fully compartmentalized reactions
  (see `#1 <https://github.com/mmundy42/cobrababel/issues/1>`_)
* Fixed setting numeric fields in MetaNetX metabolites (see
  `#2 <https://github.com/mmundy42/cobrababel/issues/2>`_)
* Fixed failing tests due to missing comma (see
  `#3 <https://github.com/mmundy42/cobrababel/issues/3>`_)
* Added support for DBLINKS field in KEGG Reaction database entry

Version 0.1.2 (June 13, 2017)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Updated compare_model() with better comparison of reaction definition and better
  output formatting

Version 0.1.1 (June 7, 2017)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Added more details on how models are compared

Version 0.1.0 (May 30, 2017)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Initial version

References
----------

1. `COBRApy: COnstraints-Based Reconstruction and Analysis for Python <http://dx.doi.org/doi:10.1186/1752-0509-7-74>`_
2. `A community-driven global reconstruction of human metabolism <http://dx.doi.org/doi:10.1038/nbt.2488>`_
3. `Generation of genome-scale metabolic reconstructions for 773 members of the human gut microbiota <http://dx.doi.org/doi:doi:10.1038/nbt.3703>`_
4. `BiGG Models: A platform for integrating, standardizing, and sharing genome-scale models <http://dx.doi.org/doi:10.1093/nar/gkv1049>`_
5. `MetaNetX.org: a website and repository for accessing, analyzing, and manipulating metabolic networks <http://dx.doi.org/doi:10.1093/bioinformatics/btt036>`_
6. `KEGG: Kyoto Encyclopedia of Genes and Genomes <http://www.kegg.jp>`_
7. `Systems Biology Knowledgebase <http://kbase.us>`_
