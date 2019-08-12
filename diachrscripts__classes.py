from scipy.stats import binom

class Digest:

    # Class to represent a genomic region that corresponds to a restriction digest.

    # Class Attribute
    __active = False

    # Initializer / Instance Attributes
    def __init__(self, chromosome, start, end):
        self.chromosome = chromosome
        self.start = start
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
    def __init__(self, digets_1, digets_2, n_simple, n_twisted):
        self.digest_1 = digets_1
        self.digest_2 = digets_2
        if digets_1.chromosome == digets_2.chromosome:
            self.__cis = True
        else:
            self.__cis = False
        self.__n_simple = n_simple
        self.__n_twisted = n_twisted

    # Instance method
    def get_digest_distance(self):
        if self.__digest_distance == None:
            self.__digest_distance = self.digest_2.start - self.digest_1.end # distance between digest ends
        return self.__digest_distance

    def is_cis(self):
        return self.__cis

    def get_interaction_category(self):
        if self.digest_1.is_active():
            state_1 = 'A'
        else:
            state_1 = 'I'
        if self.digest_2.is_active():
            state_2 = 'A'
        else:
            state_2 = 'I'
        category = state_1 + state_2
        return sorted(category)[0]+sorted(category)[1]

    def get_binomial_p_value(self): # check this function for small simple and twisted read counts
        if self.__n_simple < self.__n_twisted:
            return 1 - binom.cdf(self.__n_twisted-1, self.__n_simple + self.__n_twisted, 0.5)
        else:
            return 1 - binom.cdf(self.__n_simple-1, self.__n_simple + self.__n_twisted, 0.5)
