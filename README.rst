cobra-babel
===========

Metabolic models are available from many sources but each source uses a different ID namespace. That makes
combining and comparing models from different sources difficult. cobra-babel provides functions to translate
between ID namespaces and convert metabolic models to a common ID namespace. cobra-babel supports the
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

BiGG has limited curation so some metabolites have more than one formula and/or charge. One model with really
weird formulas (hundreds of Carbon atoms), others with different number of hydrogen (but not always a
different charge), some with unknown formula (R)