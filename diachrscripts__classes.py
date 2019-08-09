class Digest:

    # Class to represent a genomic region that corresponds to a restriction digest.

    # Class Attribute
    __active = False

    # Initializer / Instance Attributes
    def __init__(self, chromosome, start, end):
        self.chromosome = chromosome
        self.sta = start
        self.end = end

    def set_active(self):
        self.__active = True

    def is_active(self):
        return self.__active


class Interaction:

    # Class to represent a genomic interaction between two restriction digests.

    # Class Attribute
    __digest_distance = None  # Distance between the two interacting digests
    __n_simple = 0            # Number of simple read pairs
    __n_twisted = 0           # Number of twisted read pairs
    __cis = None
    __simple = None

    # Initializer / Instance Attributes
    def __init__(self, digets_1, digets_2):
        self.digets_1 = digets_1
        self.digets_2 = digets_2
        if digets_1.chromosome == digets_2.chromosome:
            self.__cis = True
        else:
            self.__cis = False

    # Instance method
    def get_digest_distance(self):
        if self.__digest_distance == None:
            self.__digest_distance = self.digets_2.end - self.digets_1.end
        return self.__digest_distance

    def is_cis(self):
        return self.__cis

    def get_interaction_category(self):
        if self.digets_1.is_active():
            state_1 = 'A'
        else:
            state_1 = 'I'
        if self.digets_2.is_active():
            state_2 = 'A'
        else:
            state_2 = 'I'
        category = state_1 + state_2
        return sorted(category)[0]+sorted(category)[1]

    def get_binomial_p_value(self):
        pass