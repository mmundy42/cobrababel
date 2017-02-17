from warnings import warn


class KeggReaction:
    
    """ A KEGG reaction record contains the following fields:

        Entry   The KEGG REACTION database is a manually curated collection of biochemical reactions,
                mostly enzymatic reactions but including some spontaneous reactions, which are derived
                from KEGG ENZYME (Enzyme Nomenclature) or KEGG PATHWAY (KEGG metabolic pathway maps).
                Each entry is identified by the unique identifier called the R number ('R' followed
                by five-digit number).
        Name    The name of the reaction, usually the systematic name of the enzyme that catalyzes the
                reaction (can be not specified).
        Definition    The chemical reaction in the form of an equation. The reaction is assumed to be
                reversible and reactants (substrates and products) are separated by '<=>'. Each compound
                in the left or the right side is separated by ' + '. There may be a coefficient before
                the compound name.
        Equation    The C number representation of the reaction equation with links to the COMPOUND
                database entries.
        Remark  The same reaction entries, if any, are given.
        Comment Text information commenting on the reaction.
        RPair   Links to the corresponding KEGG RPAIR database entries, which contain chemical structure
                transformation patterns of reactant pairs (substrate-product pairs).
        Reaction class  Links to the corresponding KEGG RCLASS database entries, which contain reaction
                class information defined by chemical structure transformation patterns of substrate-
                product pairs.
        Enzyme  Links to the corresponding KEGG ENZYME database entries.
        Pathway Links to the KEGG pathway maps, where the corresponding reaction is marked in red.
        Orthology    Links to the corresponding KEGG ORTHOLOGY (KO) database entries.
        Module  Links to the corresponding KEGG MODULE database entries.
        Reference    References for the reaction.
    """
    
    def __init__(self, record=None):
        """ Initialize object.

        Parameters
        ----------
        record : list of str, optional
            List of lines with reaction record from flat file database
        """

        self.id = None
        self.name = list()
        self.definition = None
        self.equation = None
        self.remark = None
        self.comment = list()
        self.rpair = list()
        self.rclass = list()
        self.enzyme = list()
        self.pathway = list()
        self.orthology = list()
        self.reference = list()
        self.module = list()

        if record is not None:
            self.parse(record)
        return
    
    def as_dict(self):
        """ Return the Reaction object as a dictionary.

        Returns
        -------
        dict
            Dictionary representing Reaction object
        """

        reaction = dict()
        if self.id is not None:
            reaction['id'] = self.id
        if len(self.name) > 0:
            reaction['name'] = self.name
        if self.definition is not None:
            reaction['definition'] = self.definition
        if self.equation is not None:
            reaction['equation'] = self.equation
        if self.remark is not None:
            reaction['remark'] = self.remark
        if len(self.comment) > 0:
            reaction['comment'] = self.comment
        if len(self.rpair) > 0:
            reaction['rpair'] = self.rpair
        if len(self.rclass) > 0:
            reaction['rclass'] = self.rclass
        if len(self.enzyme) > 0:
            reaction['enzyme'] = self.enzyme
        if len(self.pathway) > 0:
            reaction['pathway'] = self.pathway
        if len(self.orthology) > 0:
            reaction['orthology'] = self.orthology
        if len(self.reference) > 0:
            reaction['reference'] = self.reference
        if len(self.module) > 0:
            reaction['module'] = self.module

        return reaction
    
    def parse(self, record):
        """ Parse a record from the flat file database to complete this object.

        Parameters
        ----------
        record : list of str
            List of lines with reaction record from flat file database
        """

        for line in record:
            if line[:3] == '///':  # End of record delimiter
                return
            
            # A field in the record has the name in the first 12 characters of the line.
            # If the 12 character prefix is blank, this line is a part of the current field.
            if line[:12] != '            ':
                field_name = line[:12].strip()
            value = line[12:].strip()
            
            # Process each field in the record.
            if field_name == 'ENTRY':
                parts = value.split()
                self.id = parts[0]
            elif field_name == 'NAME':
                # A continuation of the field is delimited by a semicolon at the end of the line.
                self.name.append(value.strip(';'))
            elif field_name == 'DEFINITION':
                self.definition = value
            elif field_name == 'EQUATION':
                self.equation = value
            elif field_name == 'REMARK':
                self.remark = value
            elif field_name == 'COMMENT':
                self.comment.append(value)
            elif field_name == 'RPAIR':
                parts = value.split()
                self.rpair.append(tuple(parts))
            elif field_name == 'RCLASS':
                parts = value.split()
                self.rclass.append(tuple(parts))
            elif field_name == 'ENZYME':
                # Multiple enzymes are all listed on the same line.
                self.enzyme.extend(value.split())
            elif field_name == 'PATHWAY':
                self.pathway.append([value[:7], value[8:].strip()])
            elif field_name == 'ORTHOLOGY':
                self.orthology.append([value[:6], value[7:].strip()])
            elif field_name == 'REFERENCE':
                self.reference.append(value)
            elif field_name == 'MODULE':
                self.module.append([value[:7], value[9:].strip()])
            else:
                warn('Skipping field {0} with value {1}'.format(field_name, value))
        return

    def make_record(self):
        """ Make an reaction record for the flat file database.

        Returns
        -------
        list of str
            List of lines with reaction record
        """

        record = list()
        record.append('ENTRY       {0}                      Reaction'.format(self.id))
        if len(self.name) > 0:
            line = 'NAME        {0}'.format(self.name[0])
            if len(self.name) == 1:
                record.append(line)
            else:
                record.append(line+';')
                for index in range(1, len(self.name)-1):
                    record.append('            {0};'.format(self.name[index]))
                record.append('            {0}'.format(self.name[-1]))
        if self.definition is not None:
            record.append('DEFINITION  {0}'.format(self.definition))
        if self.equation is not None:
            record.append('EQUATION    {0}'.format(self.equation))
        if self.remark is not None:
            record.append('REMARK      {0}'.format(self.remark))
        if len(self.comment) > 0:
            record.append('COMMENT     {0}'.format(self.comment[0]))
            for index in range(1, len(self.comment)):
                record.append('            {0}'.format(self.comment[index]))
        if len(self.rpair) > 0:
            line = 'RPAIR       {0}  '.format(self.rpair[0][0])
            line += ' '.join(self.rpair[0][1:])
            record.append(line)
            for index in range(1, len(self.rpair)):
                line = '            {0}  '.format(self.rpair[index][0])
                line += ' '.join(self.rpair[index][1:])
                record.append(line)
        if len(self.rclass) > 0:
            line = 'RCLASS      {0}  '.format(self.rclass[0][0])
            line += ' '.join(self.rclass[0][1:])
            record.append(line)
            for index in range(1, len(self.rclass)):
                line = '            {0}  '.format(self.rclass[index][0])
                line += ' '.join(self.rclass[index][1:])
                record.append(line)
        if len(self.enzyme) > 0:
            enzymes = '       '.join(self.enzyme)
            record.append('ENZYME      '+enzymes)
        if len(self.pathway) > 0:
            record.append('PATHWAY     {0}  {1}'.format(self.pathway[0][0], self.pathway[0][1]))
            for index in range(1, len(self.pathway)):
                record.append('            {0}  {1}'.format(self.pathway[index][0], self.pathway[index][1]))
        if len(self.orthology) > 0:
            record.append('ORTHOLOGY   {0}  {1}'.format(self.orthology[0][0], self.orthology[0][1]))
            for index in range(1, len(self.orthology)):
                record.append('            {0}  {1}'.format(self.orthology[index][0], self.orthology[index][1]))
        if len(self.reference) > 0:
            record.append('REFERENCE   {0}'.format(self.reference[0]))
            for index in range(1, len(self.reference)):
                record.append('            {0}'.format(self.reference[index]))
        if len(self.module) > 0:
            record.append('MODULE      {0}  {1}'.format(self.module[0][0], self.module[0][1]))
            for index in range(1, len(self.module)):
                record.append('            {0}  {1}'.format(self.module[index][0], self.module[index][1]))
        record.append('///')
        return record
