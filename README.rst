cobrababel: COBRA Model Translator
==================================

Metabolic models are available from many sources and different sources use different ID namespaces. That makes
combining and comparing models from different sources difficult. cobrababel provides functions to translate
between ID namespaces and convert metabolic models to a common ID namespace. cobrababel supports the
following namespaces:

1. `ModelSEED <http://modelseed.org>`_
2. `Virtual Metabolic Human <http://vmh.uni.lu/>`_ (VMH)
3. `Biochemical, Genetic and Genomic knowledge base <http://bigg.ucsd.edu/>`_ (BiGG)

Features include:

1. Take a list of models, open all of them, create a universal model with all unique reactions and metabolites.
2. Download specific files, can get the whole BiGG database, VMH mapping file, MetaNetX translation file.
3. Function to convert compartment suffixes.
4. Import from a Excel file or csv file.
5. ID match without special characters

BiGG Notes
~~~~~~~~~~

Universal metabolites can have more than one formula if different models use a different formula
for a metabolite with the same ID (same situation can occur with charge). Currently, cobrababel blindly
picks the first formula in the list of formulae for a universal metabolite.

Universal metabolites have a list of all the compartments where the metabolite is used. Currently,
cobrababel assumes that all compartment IDs are defined the same across all models. For example, compartment
ID "c" means "cytosol" in all BiGG models. The first instance of the compartment in a metabolite is the one
that sets the name.

MetaNetX Notes
~~~~~~~~~~~~~~

Some reactions are defined with unspecified stoichiometry coefficients, for example:
`(2n) MNXM1471 + 1 MNXM3341 = (2n) MNXM1 + 1 MNXM4074 + (2n) MNXM537`. Currently, cobrababel
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
