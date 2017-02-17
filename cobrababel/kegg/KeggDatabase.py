from cobra.core import DictList


class RecordNotFound(Exception):
    """ Exception raised when a record is not found. """
    pass


class DatabaseError(Exception):
    """ Exception raised when there is a problem with a database. """
    pass


class KeggDatabase(object):
    """ Base class for managing a KEGG flat file database. """

    def __init__(self, filename):
        """ Initialize object.

        Parameters
        ----------
        filename : str
            Path to database file
        """

        self.filename = filename
        self.records = DictList()
        return
    
    def get_record(self, handle):
        """ Get a record from a database file.
        
        Parameters
        ----------
        handle : file handle
            File handle of database file

        Returns
        -------
        list of str
             List of lines in record
        """
        
        record = list()
        for line in handle:
            record.append(line.strip('\n'))
            if line[:3] == '///':
                yield record
                record = list()
                continue

    def store(self, filename=None):
        """ Save the database to a flat file.

        Parameters
        ----------
        filename : str, optional
            Path to database file
\       """

        if filename is None:
            filename = self.filename

        # Convert all of the record objects to flat file database records and write to the file.
        self.records.sort()
        with open(filename, 'w') as handle:
            for index in range(len(self.records)):
                for line in self.records[index].make_record():
                    handle.write(line + '\n')
        return

    def update(self, new_object):
        """ Update a record in the database (add new or replace existing record).

        Parameters
        ----------
        new_object : object
            Record object to add or replace
        """

        # Replace the current object if it already exists in the database.
        if self.records.has_id(new_object.id):
            self.records._replace_on_id(new_object)
            return

        # Add the new object to the database.
        self.records += [new_object]
        return

    def size(self):
        """ Get the number of records in the database.

        Returns
        -------
        int
            Number of records in database
        """

        return len(self.records)

    def has_id(self, id):
        """ Check if an ID exists in the database.

        Parameters
        ----------
        id : str
            ID to check

        Returns
        -------
        bool
            True when ID exists, otherwise False
        """

        return self.records.has_id(id)

    def get_by_id(self, id):
        """ Get an record with the specified ID.

        Parameters
        ----------
        id : str
            ID of record to return

        Returns
        -------
        object
            Object with specified ID
        """

        return self.records.get_by_id(id)
