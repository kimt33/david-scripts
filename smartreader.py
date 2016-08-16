""" The aim is to find the value associated with some regular expression which is specific to a sequence of regular expressions
i.e. We want to find some value in the text: Though we know the regular expression to find that value in that line, there may be
many lines that satisfies the regular expression. We can differentiate between them by some flag, a certain pattern that occurs
before the desired regular expression.

example:
Header0
sometext: item
sometext: item
Header1
sometext: item <-- we want to find this value: simply using r'^sometext: (item)' wont do
"""

class Header:
    """ landmarks to an input that would help find the right item

    Attributes:
        regobj: compiled regular expression object
        header: instance of Header that is the Header of self
        headerflag: boolean object that describe weather the header has been found

    Note: a header can have a header which can have a header ...
    """
    def __init__(self, regobj, header=None, reg_ender=None):
        self.regobj = regobj
        self.header = header
        if self.header:
            self.headerflag = False
        else:
            self.headerflag = True
        self.reg_ender = reg_ender

    def update_header(self, line):
        """ updates the appropriate headerflag if any of the headers (all headers of headers...) is found

        Arg:
            line: a string
        """
        # if header is not found yet and header exists
        if not self.headerflag and self.header is not None:
            # if header's header exists
            if self.header.headerflag:
                match = self.header.regobj.search(line)
                # if header's header has been found
                if match:
                    self.headerflag = True
            #if header's header doesn't exist
            else:
                # update to see if header's header exists
                self.header.update_header(line)

    def end_header(self, line):
        if self.reg_ender and self.reg_ender.search(line):
            self.headerflag=False
        if self.header and self.headerflag==False:
            self.header.end_header(line)

class Item(Header):
    """
    Attributes:
        regobj: compiled regular expression object
        header: instance of Header that is the Header of self
        headerflag: boolean object that describe weather the header has been found
        value: list of list of values that have been matched (one for each header)
        max_match: integer that sets the maximum number of matches allowed in header (set to -1 if limitless)
        num_match: counter for number of matches occured in header
        head_count: counter for header
        reg_ender: compiled regular expression object that would disable header

    Note:
        value is list of list because multiple header business
    """

    def __init__(self, regobj, header=None, max_match=1, reg_ender=None):
        self.regobj = regobj
        self.value = []
        self.header = header
        self.reg_ender = reg_ender
        if self.header:
            self.headerflag = False
        else:
            self.headerflag = True
        self.max_match = max_match
        self.num_match = 0
        self.head_count = -1

    def __call__(self, line):
        """ finds the right item, append the matched group to the value.

        Note: headerflag is set to False after matching certain number of times (unless self.max_match is set to anything that is not positive int).
        Note: if you are going to use an item as a header, i.e. you want to store something from the header, then you should
        call the last item first, then move up (because of the headerflag reset)
        """
        # INITIALIZE  (i.e. reset value, head_count, num_match, headerflag)
        # if header exists and was found in current line
        # OR
        # if it doesn't exist and item was found
        if (self.header and self.header.regobj.search(line)) or \
           (self.header is None and self.regobj.search(line)):
            # NOTE: in case of multiple identical headers, each header would initiate value and items are added into it
            self.value.append([])
            self.head_count += 1
            self.num_match = 0
            self.headerflag = True

        # FIND (i.e. find item or header)
        # if header was found already
        if self.headerflag:
            match = self.regobj.search(line)
            # if item was found
            if match:
                self.num_match += 1
                # if item was found max_match number of times (in total, NOT within given header)
                if self.num_match == self.max_match:
                    self.headerflag = False
                    self.num_match = 0
                self.value[self.head_count].append(match.groups())
            # ends header only if items have been found
            if len(self.value[-1]) > 0:
                self.end_header(line)
        # if header was not found
        else:
            self.update_header(line)

    def clear(self):
        """ clears the attributes

        try not to use: ugly.
        """
        self.value = []
        if self.header:
            self.headerflag = False
        else:
            self.headerflag = True
        self.num_match = 0
        self.head_count = -1
