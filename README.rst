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

Notes
~~~~~

1. In BiGG, universal metabolites can have more than one formula if different models use a different formula
   for a metabolite with the same ID (same situation can occur with charge). Currently, cobrababel blindly
   picks the first formula in the list of formulae for a universal metabolite.
2. In BiGG, universal metabolites have a list of all the compartments where the metabolite is used. Currently,
   cobrababel assumes that all compartment IDs are defined the same across all models. For example, compartment
   ID "c" means "cytosol" in all BiGG models.
