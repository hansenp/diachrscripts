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

        # Map status flag pair to int
        if statusA == 'N' and statusB == 'N':
            self.s_flag_int = 0
        elif statusA == 'N' and statusB == 'E':
            self.s_flag_int = 1
        elif statusA == 'E' and statusB == 'N':
            self.s_flag_int = 2
        elif statusA == 'E' and statusB == 'E':
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
    def fromA(self):
        return self._fromA

    @property
    def toA(self):
        return self._toA

    @property
    def chrB(self):
        return "chr%s" % self._chrB

    @property
    def fromB(self):
        return self._fromB

    @property
    def toB(self):
        return self._toB

    @property
    def n_simple(self):
        return self._simple

    @property
    def n_twisted(self):
        return self._twisted

    @property
    def rp_total(self):
        return self._simple + self._twisted

    @property
    def enrichment_status_tag_pair(self):
        if self.s_flag_int == 0:
            return 'NN'
        elif self.s_flag_int == 1:
            return 'NE'
        elif self.s_flag_int == 2:
            return 'EN'
        else:
            return 'EE'

    @property
    def key(self):
        """
        Unique key to refer to this pair of Digests
        """
        return "%s:%d:%s:%d" % (self._chrA, self._fromA, self._chrB, self._fromB)

    def get_diachromatic_interaction_line(self):
        """
        :return: A string that represents this interaction in Diachromatic's interaction format
        """

        # Determine enrichment state of the two digests
        if self.s_flag_int == 0:
            statusA = 'N'
            statusB = 'N'
        elif self.s_flag_int == 1:
            statusA = 'N'
            statusB = 'E'
        elif self.s_flag_int == 2:
            statusA = 'E'
            statusB = 'N'
        else:
            statusA = 'E'
            statusB = 'E'

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

        if type(self).__name__ == "DiachromaticInteraction11":
            di_line = di_line + "\t" + str(self._nln_pval) + "\t" + self.get_category()

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


class DiachromaticInteraction11(DiachromaticInteraction):
    """
    This is an extension of the class DiachromaticInteraction that can store also a P-value and a category.
    """

    def __init__(self, chrA: str, fromA: int, toA: int, statusA:str, chrB: str, fromB: int, toB: int, statusB: str, simple: int, twisted: int, nln_pval: float):
        super().__init__(chrA, fromA, toA, statusA, chrB, fromB, toB, statusB, simple, twisted)

        # Negative of the natural logarithm of P-values
        self._nln_pval = nln_pval

        # Interaction category
        self._category = None

    def set_category(self, category: str):
        """
        Sets the category of an interaction.
        :param category: A string, either UI, DI or UIR
        """

        if category == "UI":
            self._category = 0
        elif category == "DI":
            self._category = 1
        elif category == "UIR":
            self._category = 2
        else:
            raise TypeError("Invalid tag for interaction category: " + category + ". Must be either 'UI', 'DI' or 'UIR'")

    def get_category(self):
        """
        Returns the category of an interaction.
        :param category: A string, either UI, DI or UIR
        """

        if self._category == 0:
            return "UI"
        elif self._category == 1:
            return "DI"
        elif self._category == 2:
            return "UIR"
        else:
            raise NameError("Interaction category not yet defined.")
