class DiachromaticInteraction:
    """
    Class to represent an interaction between two different parts of the genome with status and count.
    Each instance of this class represents one line of the Diachromatic output file.
    See https://github.com/TheJacksonLaboratory/diachromatic
    """

    def __init__(self, chrA: str, fromA: int, toA: int, statusA:str, chrB: str, fromB: int, toB: int, statusB: str, simple: int, twisted: int):
        """
        @param F an array with 9 elements that represent an interaction
        """

        # map status flag pair to int
        if statusA == 'I' and statusB == 'I':
            self.s_flag_int = 0
        elif statusA == 'I' and statusB == 'A':
            self.s_flag_int = 1
        elif statusA == 'A' and statusB == 'I':
            self.s_flag_int = 2
        elif statusA == 'A' and statusB == 'A':
            self.s_flag_int = 3
        else:
            raise ValueError("Bad status values: A-%s, B-%s" % (statusA, statusB))

        if chrA.startswith("chr"):
            self._chrA = chrA[3:]
        else:
            self._chrA = chrA
        self._fromA = fromA
        self._toA = toA

        if chrB.startswith("chr"):
            self._chrB = chrB[3:]
        else:
            self._chrB = chrB
        self._fromB = fromB
        self._toB = toB

        self._simple = simple
        self._twisted = twisted

        self._rep_num = 1
        

    def append_interaction_data(self, simple, twisted):
        """
        This method is called if we already have one or more observations
        for this interaction. We have thus compared the key already.
        We can just add the interaction counts.
        """
        self._rep_num += 1
        self._simple += simple
        self._twisted += twisted

    def has_data_for_required_replicate_num(self, k):
        return k <= self._rep_num

    @property
    def chrA(self):
        return "chr%s" % self._chrA

    @property
    def chrB(self):
        return "chr%s" % self._chrB

    @property
    def simple(self):
        return self._simple

    @property
    def twisted(self):
        return self._twisted

    @property
    def total_readpairs(self):
        return self._simple + self._twisted

    @property
    def status(self):
        if self.s_flag_int == 0:
            return 'II'
        elif self.s_flag_int == 1:
            return 'IA'
        elif self.s_flag_int == 2:
            return 'AI'
        else:
            return 'AA'


    @property
    def key(self):
        """
        Unique key to refer to this pair of Digests
        """
        return "%s%d%s%d" % (self._chrA, self._fromA, self._chrB, self._fromB)


    def output_summary(self):
        """
        This function essentially recreates the original line from Diachromatic (for a single interaction)
        or creates an analogous line for the recreated interactions.
        """
        if self.s_flag_int == 0:
            statusA = 'I'
            statusB = 'I'
        elif self.s_flag_int == 1:
            statusA = 'I'
            statusB = 'A'
        elif self.s_flag_int == 2:
            statusA = 'A'
            statusB = 'I'
        else:
            statusA = 'A'
            statusB = 'A'
        return "chr{}\t{}\t{}\t{}\tchr{}\t{}\t{}\t{}\t{}:{}".format(self._chrA,
                                                                self._fromA,
                                                                self._toA,
                                                                statusA,
                                                                self._chrB,
                                                                self._fromB,
                                                                self._toB,
                                                                statusB,
                                                                self._simple,

                                                                self._twisted)
    def get_diachromatic_interaction_line(self):
        """
        :return: A string that represents this interaction in Diachromatic's interaction format
        """

        # Determine enrichment state of the two digests
        if self.s_flag_int == 0:
            statusA = 'I'
            statusB = 'I'
        elif self.s_flag_int == 1:
            statusA = 'I'
            statusB = 'A'
        elif self.s_flag_int == 2:
            statusA = 'A'
            statusB = 'I'
        else:
            statusA = 'A'
            statusB = 'A'

        di_line = \
            self.chrA + "\t" + \
            str(self._fromA) + "\t" + \
            str(self._toA) + "\t" + \
            statusA + "\t" + \
            self.chrB + "\t" + \
            str(self._fromB) + "\t" + \
            str(self._toB) + "\t" + \
            statusB + "\t" + \
            str(self._simple)+ ":" + \
            str(self._twisted)

        return di_line

    def __hash__(self):
        """
        Note that it is sufficient to know the chromosome and the start position to 
        unambiguously identify the digest. Therefore, the has function only needs
        to use the following four items.
        """
        return hash((self._chrA, self._fromA, self._chrB, self._fromB))

    def __eq__(self, other):
        return self._fromA == other._fromA and self._fromB == other._fromB and self._chrA == other._chrA and self._chrB == other._chrB
