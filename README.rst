CobraBabel: COBRA Model Translator
==================================

Metabolic models are available from many sources and different sources use different ID namespaces. That makes
combining and comparing models from different sources difficult. CobraBabel provides functions to translate
between ID namespaces and convert metabolic models to a common ID namespace. CobraBabel supports the
following namespaces:

1. `ModelSEED <http://modelseed.org>`_
2. `Virtual Metabolic Human <http://vmh.uni.lu/>`_ (VMH)
3. `Biochemical, Genetic and Genomic knowledge base <http://bigg.ucsd.edu/>`_ (BiGG)
4. `MetaNetX <http://www.metanetx.org/>`_
5. `Kyoto Encyclopedia of Genes and Genomes <http://www.kegg.jp>`_ (KEGG)

Features include:

1. Take a list of models, open all of them, create a universal model with all unique reactions and metabolites.
2. Download specific files, can get the whole BiGG database, VMH mapping file, MetaNetX translation file.
3. Function to convert compartment suffixes.
4. Import from a Excel file or csv file.
5. ID match without special characters

VMH Notes
^^^^^^^^^

VMH provides a SEED2VMH_translation.csv file which is not current with the translation table
in Supplementary Table 15 in the AGORA paper. Supplementary Table 15 is provided as an
Excel workbook which was parsed into two files: (1) vmh_metabolite_xref.tsv and
(2) vmh_reaction_xref.tsv which are provided in the "data" folder.

BiGG Notes
^^^^^^^^^^

Universal metabolites can have more than one formula if different models use a different formula
for a metabolite with the same ID (same situation can occur with charge). Currently, CobraBabel blindly
picks the first formula in the list of formulae for a universal metabolite.

Universal metabolites have a list of all the compartments where the metabolite is used. Currently,
CobraBabel assumes that all compartment IDs are defined the same across all models. For example, compartment
ID "c" means "cytosol" in all BiGG models. The first instance of the compartment in a metabolite is the one
that sets the name.

MetaNetX Notes
^^^^^^^^^^^^^^

Some reactions are defined with unspecified stoichiometry coefficients, for example:
`(2n) MNXM1471 + 1 MNXM3341 = (2n) MNXM1 + 1 MNXM4074 + (2n) MNXM537`. Currently, CobraBabel
does not include these reactions when creating an universal model.

There are some reactions where there is a BIOMASS metabolite in the reaction definition.
But the BIOMASS metabolite is not defined so any reactions with the BIOMASS metabolite are
not included in the universal model.

If a reaction has a value in the Source field, there is only one source in the format:
`source:id`. Sources include Rhea (rhea), KEGG (kegg), MetaCyc (metacyc), UniPathway (upa),
The Seed (seed), BiGG (bigg), BioPath (biopath), and Reactome (reactome). Set the `verbose`
parameter when calling `create_metanetx_universal_model()` to show a warning for reactions with
an invalid format in the Source field.

There are no compartments so metabolites are defined without a compartment. Still not sure
how exchange and transport reactions are defined.

All reactions are defined as bi-directional so the lower and upper bounds are set to the
default values `(-1000.0, 1000.0)` for all reactions.

KEGG Notes
^^^^^^^^^^

There are some limitations on the `KEGG API <http://www.kegg.jp/kegg/rest/>`_, most importantly:

**Restrictions:** KEGG API is provided for academic use by academic users belonging to academic
institutions. This service should not be used for bulk data downloads. Please obtain KEGG FTP
academic subscription for downloading KEGG data.

When setting OTU representatives in the organism database, there can be variability in which
organism is selected when there are multiple equally good matches. For example, in the KEGG database
there are two Mycobacterium tuberculosis H37Rv organisms, one with ID "T00015" and one with ID
"T02178". It is random as to which one gets picked as the OTU representative.

A sample OTU representative source file is provided in the "data/otu-reps.tsv" file which was
retrieved from the KBase Central Data Model using the all_entities_OTU command. No idea how to
do this with the KBase narrative interface. Maybe need to use the data_api service?


References
----------

BiGG

MetaNetX.org: a website and repository for accessing, analyzing and manipulating metabolic networks.
Bioinformatics (2013)

VMH

Systems Biology Knowledgebase http://kbase.us