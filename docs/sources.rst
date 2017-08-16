
Create models from source systems
---------------------------------

CobraBabel supports creating COBRA models from multiple source systems.

VMH
~~~

Virtual Metabolic Human includes a model for human metabolism (Recon2)
and models for human gut microbes (AGORA).

.. code:: ipython3

    btheta = cobrababel.create()

BiGG
~~~~

Biochemical, Genetic, and Genomic knowledge base (BiGG) includes
manually curated models for microbes along with universal databases of
reactions, metabolites, and compartments.

MetaNetX
~~~~~~~~

MetaNetX includes universal databases of reactions, metabolites, and
compartments along with a reaction and a metabolite cross reference to
map between IDs between source systems.

KEGG
~~~~

Kyoto Encyclopedia of Genes and Genomes (KEGG) includes a universal
database of reactions and metabolites.
