from fuzzywuzzy import fuzz
from pandas import DataFrame, Series

from .KeggDatabase import KeggDatabase, DatabaseError
from .KeggOrganism import KeggOrganism
from .kegg import list_kegg_ids

# Cutoff for identifying a good match when comparing names.
MATCH_CUTOFF = 97

# Column names for reporting results of setting OTU representatives.
otu_report_columns = ['OTU_REP_NAME', 'OTU_REP_ID', 'NAME', 'ID', 'RATIO',
                      'PARTIAL_RATIO', 'TOKEN_SET_RATIO', 'RESULT']


class KeggOrganismDatabase(KeggDatabase):
    """ Manage a KEGG organism flat file database. """

    def __init__(self, filename):
        """ Initialize object.
        
        Parameters
        ----------
        filename : str
            Path to organism database file
        """

        super(KeggOrganismDatabase, self).__init__(filename)
        self.code_to_id = dict()  # Keyed by organism code
        return

    def load(self):
        """ Load the organism database from a flat file. """

        with open(self.filename, 'r') as handle:
            # Organism database is different and each record is one tab-delimited line.
            organisms = [KeggOrganism(record) for record in handle]
        self.records += organisms

        # Build a dictionary to lookup organisms by code.
        for index in range(len(self.records)):
            self.code_to_id[self.records[index].code] = self.records[index].id
        return

    def download(self):
        """ Download the organism database from KEGG web service.

            The only way to get_kegg_records taxonomic info for organisms is with the "list" operation
            from the KEGG web service.
        """

        # Make sure the database is empty.
        if len(self.records) > 0:
            raise DatabaseError('Organism database must be empty before downloading from web service')

        # Get the current organism database from KEGG.
        organism_list = list_kegg_ids('organism')

        # Add all of the organisms to the database.
        organisms = list()
        for index in range(len(organism_list)):
            record = organism_list[index]
            record += '\t0'  # Mark as not a representative by default
            organisms.append(KeggOrganism(record))
        self.records += organisms

        # Build a dictionary to lookup organisms by code.
        for index in range(len(self.records)):
            self.code_to_id[self.records[index].code] = self.records[index].id

        return

    def set_representatives(self, otu_file):
        """ Set the otu_representative flag for organisms from source file.

        The OTU file lists the organisms that are members of an OTU and which organism
        is the OTU representative. Each line in the file has four fields: (1) OTU
        number, (2) boolean flag that identifies OTU representative ('1' means true),
        (3) genome ID, and (4) scientific name of organism.

        A DataFrame is returned with a report that gives details on the analysis done
        to pick the OTU representative organisms. Each row in the DataFrame contains
        the following columns: (1) 'OTU_REP_NAME' name of organism that is OTU representative
        in source file, (2) 'OTU_REP_ID' ID of organism that is OTU representative in source
        file, (3) 'NAME' name of organism that is a possible match, (4) 'ID' ID of organism
        that is a possible match, (5) 'RATIO' value of match ratio, (6) 'PARTIAL_RATIO'
        value of partial match ratio, (7) 'TOKEN_SET_RATIO' value of token set match ratio,
        (8) 'RESULT' result of analysis where 'selected' means the organism was selected as
        the OTU representative, 'alternate' means the organism was a good match but was not
        selected, 'skipped' means the organism was not selected, 'not found' means a good
        match was not found, and 'no match' means the organism was not considered.

        Note that if multiple organisms are equally good matches, it is random which one
        is picked as the representative.

        Parameters
        ----------
        otu_file : str
            Path to file that identifies OTU representatives in TSV format

        Returns
        -------
        pandas.DataFrame
            Report with results of analysis
        """

        if len(self.records) == 0:
            raise DatabaseError('Organism database is empty')

        # Reset all organisms to not be a OTU representative.
        for organism in self.records:
            organism.otu_representative = 0

        # Build a dictionary where the key is name of OTU representative organism and
        # value is organism ID.
        representatives = dict()
        with open(otu_file, 'r') as handle:
            for line in handle:
                fields = line.strip().split('\t')
                if fields[1] == '1':
                    # Remove 'substr' and 'str' from scientific name for better matching to KEGG names.
                    name = fields[3].replace('substr.', '').replace('str.', '')
                    representatives[name] = fields[2]

        # Build a dictionary for finding prokaryote organisms by name where the key
        # is the genus and species portion of the organism name and the value is a
        # list of KeggOrganism objects that match the genus and species. This reduces
        # the number of organisms that need to have a match ratio calculated.
        lookup = dict()
        for organism in self.records:
            if organism.is_prokaryote:
                num_words = 2
                if organism.name.startswith('Candidatus'):
                    num_words = 3  # Need first three words for candidate names
                parts = organism.name.split()
                lookup_name = ' '.join(parts[0:num_words])
                try:
                    lookup[lookup_name].append(organism)
                except KeyError:
                    lookup[lookup_name] = [organism]

        # For each OTU representative from OTU file, see if there is a matching organism
        # from KEGG. If there is a match, mark the organism as a OTU representative.
        report = DataFrame(columns=otu_report_columns)
        for name in representatives:
            num_words = 2
            if name.startswith('Candidatus'):
                num_words = 3  # Need first three words for candidate names
            parts = name.split()
            lookup_name = ' '.join(parts[0:num_words])

            # See if there are potential matches in database.
            if lookup_name in lookup:
                # Calculate the match ratio for all potential matches using the search name.
                matches = list()
                for organism in lookup[lookup_name]:
                    match = dict()
                    match['organism'] = organism
                    match['ratio'] = fuzz.ratio(name, organism.search_name)
                    match['partial_ratio'] = 0
                    match['token_set_ratio'] = 0
                    match['possible'] = False
                    matches.append(match)

                # Sort the matches by the ratio so best matches are first.
                matches = sorted(matches, key=lambda k: k['ratio'], reverse=True)

                # See if there are any good matches and select the best matches.
                possible = list()
                max_ratio = matches[0]['ratio']
                if max_ratio >= MATCH_CUTOFF:
                    for match in matches:
                        if match['ratio'] == max_ratio:
                            possible.append(match)
                            match['possible'] = True

                # If there are no good matches yet, calculate the partial match ratio
                # for all potential matches using the search name.
                if len(possible) == 0:
                    for match in matches:
                        match['partial_ratio'] = fuzz.partial_ratio(name, match['organism'].search_name)
                        if match['partial_ratio'] >= MATCH_CUTOFF:
                            possible.append(match)
                            match['possible'] = True

                # If there are no good matches yet, calculate the token set ratio
                # for all potential matches using the name.
                if len(possible) == 0:
                    for match in matches:
                        match['token_set_ratio'] = fuzz.token_set_ratio(name, match['organism'].name)
                        if match['token_set_ratio'] == 100:
                            possible.append(match)
                            match['possible'] = True

                # Set the OTU representative based on the best match.
                none_selected = False
                if len(possible) > 0:
                    possible[0]['organism'].otu_representative = 1
                else:
                    # There is only possible match that is very close.
                    if len(lookup[lookup_name]) == 1:
                        lookup[lookup_name][0].otu_representative = 1
                    else:
                        none_selected = True

                # Add the results to the report.
                for match in matches:
                    organism = match['organism']
                    if none_selected:
                        selected_value = 'not found'
                    elif len(possible) > 1:
                        if organism.is_representative:
                            selected_value = 'selected'
                        elif match['possible']:
                            selected_value = 'alternate'
                        else:
                            selected_value = 'skipped'
                    else:
                        if organism.is_representative:
                            selected_value = 'selected'
                        else:
                            selected_value = 'skipped'
                    result = Series([name, representatives[name], organism.name, organism.id,
                                    match['ratio'], match['partial_ratio'], match['token_set_ratio'],
                                    selected_value], index=otu_report_columns)
                    report = report.append(result, ignore_index=True)

            else:
                result = Series([name, representatives[name], 'no match', 'no match',
                                0, 0, 0, 'no match'], index=otu_report_columns)
                report = report.append(result, ignore_index=True)

        return report

    def num_otu_representatives(self):
        """ Get the number of OTU representative organisms in the database.

        Returns
        -------
        int
            Number of OTU representative organisms
        """

        return len(self.records.query(lambda x: x, 'is_representative'))

    def num_prokaryotes(self):
        """ Get the number of prokaryote organisms in the database.

        Returns
        -------
        int
            Number of prokaryote organisms
        """

        return len(self.records.query(lambda x: x, 'is_prokaryote'))

    def num_eukaryotes(self):
        """ Get the number of eukaryote organisms in the database.

        Returns
        -------
        int
            Number of eukaryote organisms
        """

        return len(self.records.query(lambda x: x, 'is_eukaryote'))

    def num_bacteria(self):
        """ Get the number of bacteria organisms in the database.

        Returns
        -------
        int
            Number of bacteria organisms
        """

        return len(self.records.query(lambda x: x, 'is_bacteria'))

    def num_archaea(self):
        """ Get the number of archaea organisms in the database.

        Returns
        -------
        int
            Number of archaea organisms
        """

        return len(self.records.query(lambda x: x, 'is_archaea'))

    def has_code(self, code):
        """ Check if a code exists in the organism database.

        Parameters
        ----------
        code : str
            Organism code to check

        Returns
        -------
        bool
            True when code exists, otherwise False
        """

        return code in self.code_to_id

    def get_by_code(self, code):
        """ Get an organism with the specified code.
        
        Parameters
        ----------
        code : str
            Code of record to return

        Returns
        -------
        KeggOrganism
            Organism object with specified code
        """

        return self.records.get_by_id(self.code_to_id[code])
