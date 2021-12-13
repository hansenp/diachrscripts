from numpy import array, argsort


class DiachromaticInteraction:
    """
    Class to represent an interaction between two different parts of the genome with status and count.
    Each instance of this class represents one line of the Diachromatic output file.
    See https://github.com/TheJacksonLaboratory/diachromatic
    """

    def __init__(self, chrA: str, fromA: int, toA: int, statusA: str, chrB: str, fromB: int, toB: int, statusB: str,
                 simple_1: int, simple_2: int, twisted_1: int, twisted_2: int):
        """
        @param F an array with 9 elements that represent an interaction
        """

        # Map status flag pair to int
        if (statusA == 'N' and statusB == 'N') or (statusA == 'I' and statusB == 'I'):
            self.s_flag_int = 0
        elif (statusA == 'N' and statusB == 'E') or (statusA == 'I' and statusB == 'A'):
            self.s_flag_int = 1
        elif (statusA == 'E' and statusB == 'N') or (statusA == 'A' and statusB == 'I'):
            self.s_flag_int = 2
        elif (statusA == 'E' and statusB == 'E') or (statusA == 'A' and statusB == 'A'):
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

        self._simple_1 = simple_1
        self._simple_2 = simple_2
        self._twisted_1 = twisted_1
        self._twisted_2 = twisted_2

        self._rep_num = 1

    def append_interaction_data(self, simple_1, simple_2, twisted_1, twisted_2):
        """
        This method is called if we already have one or more observations
        for this interaction. We have thus compared the key already.
        We can just add the interaction counts.
        """
        self._rep_num += 1
        self._simple_1 += simple_1
        self._simple_2 += simple_2
        self._twisted_1 += twisted_1
        self._twisted_2 += twisted_2

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
        return self._simple_1 + self._simple_2

    @property
    def n_twisted(self):
        return self._twisted_1 + self._twisted_2

    @property
    def n_heaviest_two(self):
        rp_counts = sorted([self._simple_1, self._simple_2, self._twisted_1, self._twisted_2], reverse=True)
        return rp_counts[0] + rp_counts[1]

    @property
    def n_lightest_two(self):
        rp_counts = sorted([self._simple_1, self._simple_2, self._twisted_1, self._twisted_2], reverse=True)
        return rp_counts[2] + rp_counts[3]

    @property
    def rp_total(self):
        return self._simple_1 + self._simple_2 + self._twisted_1 + self._twisted_2

    @property
    def i_dist(self):
        return self._fromB - self._toA

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
            str(self._simple_1) + ":" + \
            str(self._simple_2) + ":" + \
            str(self._twisted_1) + ":" + \
            str(self._twisted_2)

        if type(self).__name__ == "DiachromaticInteraction11":
            log10_pval = round(self._log10_pval, 2)
            log10_pval += 0.  # Get rid of negative sign on numbers close to zero
            di_line = di_line + "\t" + "{:.2f}".format(log10_pval) + "\t" + self.get_category()

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

    def __init__(self, chrA: str, fromA: int, toA: int, statusA: str, chrB: str, fromB: int, toB: int, statusB: str,
                 simple_1: int, simple_2: int, twisted_1: int, twisted_2: int, log10_pval: float):
        super().__init__(chrA, fromA, toA, statusA, chrB, fromB, toB, statusB, simple_1, simple_2, twisted_1, twisted_2)

        # Negative of the decadic logarithm of P-values
        self._log10_pval = log10_pval

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
        elif category == "DIX":
            self._category = 3
        else:
            raise TypeError(
                "Invalid tag for interaction category: " + category + ". Must be either 'UI', 'DI' or 'UIR'")

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
        elif self._category == 3:
            return "DIX"
        else:
            raise NameError("Interaction category not yet defined.")

    def get_pval(self):
        return 10 ** -self._log10_pval

    def get_ht_tag(self):
        """
        Returns one of the tags '01', '02', '03', '12', '13', '23' depending on which of the four read pair counts are
        the largest. For instance, for '10:29:2:7' the tag '01' would be returned because the first two read pair
        counts are the largest, whereas for '6:0:3:21' the tag '03' would be returned because the first and the last
        read pair counts are the largest.
        """

        rpc_list = [self._simple_1, self._simple_2, self._twisted_2, self._twisted_1]
        rp_total = sum(rpc_list)
        if 0.75 < self._simple_1/rp_total:
            return '0X'
        if 0.75 < self._simple_2/rp_total:
            return '1X'
        if 0.75 < self._twisted_2/rp_total:
            return '2X'
        if 0.75 < self._twisted_1/rp_total:
            return '3X'
        rpcs = -1 * array(rpc_list)
        rpcs_sort_idx = argsort(rpcs)
        ht_tag = "".join(map(str, sorted(rpcs_sort_idx[:2])))
        return ht_tag
