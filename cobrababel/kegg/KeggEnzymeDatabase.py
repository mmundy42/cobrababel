import re
from warnings import warn

from .KeggDatabase import KeggDatabase
from .KeggEnzyme import KeggEnzyme


class KeggEnzymeDatabase(KeggDatabase):
    """ Manage a KEGG enzyme flat file database. """

    def __init__(self, filename):
        """ Initialize object.
        
        Parameters
        ----------
        filename : str
            Path to enzyme database file
        """

        super(KeggEnzymeDatabase, self).__init__(filename)
        return
    
    def load(self):
        """ Load the enzyme database from a flat file. """
        
        with open(self.filename, 'r') as handle:
            enzymes = [KeggEnzyme(record) for record in self.get_record(handle)]
        self.records += enzymes
        return
    
    def store(self, filename=None):
        """ Save the enzyme database to a flat file.
        
        Parameters
        ----------
        filename : str, optional
            Path to enzyme database file
        """
        
        if filename is None:
            filename = self.filename
        
        # Convert all of the KeggEnzyme objects to flat file database records
        # and write the records to the file.
        self.records.sort()
        with open(filename, 'w') as handle:
            for index in range(len(self.records)):
                for line in self.records[index].make_record():
                    handle.write(line+'\n')
        return

    def update(self, enzyme):
        """ Update an enzyme in the database (add new or replace existing enzyme).
        
        Parameters
        ----------
        enzyme : KeggEnzyme
            Enzyme to add or replace
        """

        # Replace the current KeggEnzyme object if it already exists in the database.
        if self.records.has_id(enzyme.id):
            self.records._replace_on_id(enzyme)
            return
            
        # Add the new KeggEnzyme object to the database.
        self.records += [enzyme]
        return

    def validate(self):
        """ Validate that obsolete enzymes are linked to valid enzymes in the database.

        Returns
        -------
        dict
            Dictionary with counts and list of IDs of missing enzymes
        """

        # Keep track of what is found.
        obsolete = dict()
        obsolete['num_obsolete'] = 0
        obsolete['num_transferred'] = 0
        obsolete['num_deleted'] = 0
        obsolete['missing_enzymes'] = list()
        ec_number_re = re.compile(r'([\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+)')

        # Check each enzyme in the database to see it is obsolete.
        for enzyme in self.records:
            # enzyme = self.records.get_by_id(id)
            if enzyme.obsolete:
                obsolete['num_obsolete'] += 1
                if len(enzyme.name) > 1:
                    warn('Obsolete enzyme {0} has more than one name: {1}'.format(id, enzyme.name))

                # When an obsolete enzyme has been transferred to another enzyme, the name field
                # has the reference to the replacement enzyme.
                if enzyme.name[0].startswith('Transferred to '):
                    obsolete['num_transferred'] += 1
                    matches = re.findall(ec_number_re, enzyme.name[0])
                    for ec_number in matches:
                        if not self.records.has_id(ec_number):
                            obsolete['missing_enzymes'].append(ec_number)

                # When an obsolete enzyme has been deleted, the comment field may have a
                # reference to a replacement enzyme.
                elif enzyme.comment[0].startswith('Deleted entry: '):
                    obsolete['num_deleted'] += 1
                    matches = re.findall(ec_number_re, enzyme.comment[0])
                    for ec_number in matches:
                        if not self.records.has_id(ec_number):
                            obsolete['missing_enzymes'].append(ec_number)
                else:
                    warn('Obsolete enzyme {0} has not been transferred or deleted'.format(id))

        return obsolete
