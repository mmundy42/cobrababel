from .KeggDatabase import KeggDatabase
from .KeggReaction import KeggReaction


class KeggReactionDatabase(KeggDatabase):
    """ Manage a KEGG reaction flat file database. """

    def __init__(self, filename):
        """ Initialize object.

        Parameters
        ----------
        filename : str
            Path to reaction database file
        """

        super(KeggReactionDatabase, self).__init__(filename)
        return
    
    def load(self):
        """ Load the reaction database from a flat file. """
        
        with open(self.filename, 'r') as handle:
            reactions = [KeggReaction(record) for record in self.get_record(handle)]
        self.records += reactions
        return
    
    def store(self, file_name=None):
        """ Save the reaction database to a flat file.

        Parameters
        ----------
        file_name : str, optional
            Path to reaction database file
        """
        
        if file_name is None:
            file_name = self.filename
            
        # Convert all of the KeggReaction objects to flat file database records
        # and write the records to the file.    
        self.records.sort()
        with open(file_name, 'w') as handle:
            for index in range(len(self.records)):
                for line in self.records[index].make_record():
                    handle.write(line + '\n')
        return

    def update(self, reaction):
        """ Update a reaction in the database (add new or replace existing reaction).

        Parameters
        ----------
        reaction : kegg.KeggReaction
            Reaction to add or replace
        """

        # Replace the current KeggReaction object if it already exists in the database.
        if self.records.has_id(reaction.id):
            self.records._replace_on_id(reaction)
            return

        # Add the new KeggReaction object to the database.
        self.records += [reaction]
        return
