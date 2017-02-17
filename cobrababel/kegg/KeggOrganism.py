import re

# Regular expression to remove text in parenthesis at end of string.
# @todo Remove leading space?
paren_re = re.compile(r' *\(.*\)$')


class KeggOrganism:
    
    """ A KEGG organism record contains the following fields:

        ID      Each entry is identified by the unique identifier called the T number ('T' followed
                by five-digit number).
        Code    The code for the organism is a three or four character string based on the name
                and is used for identifying genes in complete genomes.
        Name    The scientific name of the organism (common name in parenthesis).
        Taxonomy    The taxonomy of the organism as a string with levels separated by semicolons.
        OtuRep  Flag indicating if the organism is representative of the OTU.

        Each record is a single line with fields separated by tabs (which is different from other record types).
    """
    
    def __init__(self, record=None):
        """ Initialize object.

        Parameters
        ----------
        record : list of str, optional
            List of lines with enzyme record from flat file database
        """
        self.id = None
        self.code = None
        self.name = None
        self.taxonomy = tuple()
        self.otu_representative = 0
        self.search_name = None

        if record is not None:
            self.parse(record)
        return
    
    def as_dict(self):
        """ Return the Organism object as a dictionary.

        Returns
        -------
        dict
            Dictionary representing Organism object
        """

        organism = dict()
        if self.id is not None:
            organism['id'] = self.id
        if self.code is not None:
            organism['code'] = self.code
        if self.name is not None:
            organism['name'] = self.name
        organism['taxonomy'] = self.taxonomy
        organism['otu_representative'] = self.otu_representative
        return organism
    
    def parse(self, record):
        """ Parse a record from the flat file database to complete this object.

        Parameters
        ----------
        record : str
            Line with organism record from flat file database
        """

        # Each record is a single line with the fields separated by tabs. 
        fields = record.split('\t')
        self.id = fields[0]
        self.code = fields[1]
        self.name = fields[2]
        self.taxonomy = tuple(fields[3].split(';'))  # Taxonomy levels are separated by semicolons
        self.otu_representative = int(fields[4])
        self.search_name = re.sub(paren_re, '', self.name)  # Remove text in parenthesis at end of string
        return

    def make_record(self):
        """ Make a record for the flat file database.
        
        Returns
        -------
        str
            Line with organism record
        """
        
        taxonomy = ';'.join(self.taxonomy)
        record = '\t'.join([self.id, self.code, self.name, taxonomy, str(self.otu_representative)])
        return [record]

    # Make these methods properties, then should be able to use with query() method on DictList
    @property
    def is_prokaryote(self):
        """ Check if the organism is a prokaryote.

        Returns
        -------
        bool
            True if the organism is a prokaryote
        """

        return self.taxonomy[0] == 'Prokaryotes'

    @property
    def is_eukaryote(self):
        """ Check if the organism is a eukaryote.

        Returns
        -------
        bool
            True if the organism is a eukaryote
        """
        
        return self.taxonomy[0] == 'Eukaryotes'

    @property
    def is_bacteria(self):
        """ Check if the organism is a bacteria.

        Returns
        -------
        bool
            True if the organism is a bacteria
        """
        
        return self.taxonomy[1] == 'Bacteria'

    @property
    def is_archaea(self):
        """ Check if the organism is an archaea.

        Returns
        -------
        bool
           True if the organism is an archaea
        """
        
        return self.taxonomy[1] == 'Archaea'

    @property
    def is_representative(self):
        """ Check if the organism is a representative of the OTU.

        Returns
        -------
        bool
            True if organism is an OTU representative
        """

        return self.otu_representative
