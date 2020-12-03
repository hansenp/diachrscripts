from scipy.stats import binom
from numpy import log, logaddexp
import numpy as np
import warnings

class BinomialModel:
    """
    In this class, interactions are statistically evaluated in terms of directionality. The main task of this class
    is the calculation of P-values. The distribution required for this is instantiated only once for an object of this
    class and is then used for all calculations. In addition, to save time, P-values that have already been calculated
    for given pairs of simple and twisted read pair counts are saved in a dictionary to save time.
    """
    def __init__(self, prob_simple:float, n_max:int) -> None:
        self._prob_simple = prob_simple
        self._n_max = n_max


        self._pval_dict = {}

        print("n_max: " + str(self._n_max))
        #self._precalculate_p_values_using_pmf()


        #self._pre_calculate_p_values_using_pmf_and_log()
        self.compare_different_ways_of_p_value_calculation(N=10,prob_simple=0.5)
        #self._pre_calculate_p_values_using_logcdf_and_logsf()

    def get_binomial_logsf_p_value(self, n_simple, n_twisted):
        """
        This function returns already calculated P-values from the dictionary or,
        if the P-value has not yet been calculated, a funnction is called to
        calculate the P-value. This P-value is then added to the dictionary and returned.

        :return: P-value.
        """
        # Create key from simple and twisted read pair counts
        key = "{}-{}".format(n_simple, n_twisted)

        if key not in self._pval_dict:
            self._pval_dict[key] = self._calculate_binomial_logsf_p_value(n_simple, n_twisted)

        return self._pval_dict[key]


    def _calculate_binomial_logsf_p_value(self, n_simple, n_twisted): # (natural) logsf
        """
        So far, we have used this function to calculate the P-values in script '04' and we will
        continue to uses.
        This function assumes a background frequency of 0.5 for simple and twisted read pairs.
        Further down in this class an attempt is made to allow different background frequencies
        for simple and twisted read pairs.

        :param n_simple:
        :param n_twisted:
        :return: Since the distribution is symmetrical because of p=0.5,
        we return twice the P-value from a one-sided test as the P-value for a two-sided test.
        For P-values that are too small to be represented, we return -np.inf.
        """
        try:
            if n_simple < n_twisted:
                p_value = binom.logsf(n_twisted - 1, n_simple + n_twisted, 0.5)
                return 2*p_value
            else:
                p_value = binom.logsf(n_simple - 1, n_simple + n_twisted, 0.5)
                return 2*p_value
        except RuntimeWarning: # Underflow: P-value is too small
            return -np.inf
            #return log(sys.float_info.min * sys.float_info.epsilon) # return natural log of smallest possible float

    def compare_different_ways_of_p_value_calculation(self, N=5, prob_simple=0.5):
        """
        In this function, various methods for calculating the the P-values are tried out
        and compared. First, the P-values are calculated by adding up the pmf.
        Then the P-values are calculated using the cdf.
        These P-values are then logarithmized and compared with those calculated with logcdf and logsf.

        :param N: Maximum n
        :param prob_simple: Background frequency of simple read pairs
        """

        prob_twisted = 1.0 - prob_simple
        n_list = list(range(1, N + 1))
        print("n_list: " + str(n_list))
        print("prob_simple: " + str(prob_simple))
        print("prob_twisted: " + str(prob_twisted))

        for n in n_list:

            print("--------")
            print("n: " + str(n))
            k_list = list(range(0, n + 1))
            print("k_list: " + str(k_list))

            print()

            try:
                # Get PMF for simple and twisted
                pmf_list_simple = binom.pmf(k_list, n, prob_simple)
                print("pmf_list_simple: " + str(pmf_list_simple))
                pmf_list_twisted = binom.pmf(k_list, n, prob_twisted)
                print("pmf_list_twisted: " + str(pmf_list_twisted))

                print()

                # Calculate one-sided P-values for simple by adding up the densities
                p_values_one_sided_simple = []
                d_sum = 0
                for d in pmf_list_simple:
                    d_sum = d_sum + d
                    p_values_one_sided_simple.append(d_sum)
                print("p_values_one_sided_simple: " + str(p_values_one_sided_simple))

                # Calculate one-sided P-values for twisted by adding up the densities
                p_values_one_sided_twisted = []
                d_sum = 0
                for d in pmf_list_twisted:
                    d_sum = d_sum + d
                    p_values_one_sided_twisted.append(d_sum)
                print("p_values_one_sided_twisted: " + str(p_values_one_sided_twisted))

                d_sum = 0
                p_values_two_sided = []
                for k in k_list[:int(n/2)]:
                    d_sum = d_sum + pmf_list_simple[k] + pmf_list_twisted[k]
                    p_values_two_sided.append(d_sum)
                print("p_values_two_sided: " + str(p_values_two_sided))


            except RuntimeWarning:
                print("Underflow at 'pmf' and n = " + str(n) + "!")
                exit(1)

            print()

            # Get CDF for simple and twisted
            try:
                cdf_list_simple = binom.cdf(k_list, n, prob_simple)
                print("cdf_list_simple: " + str(cdf_list_simple))
                cdf_list_twisted = binom.cdf(k_list, n, prob_twisted)
                print("cdf_list_twisted: " + str(cdf_list_twisted))

                p_values_two_sided_cdf = []
                for k in k_list[:int(n/2)]:
                    p_values_two_sided_cdf.append(cdf_list_simple[k] + cdf_list_twisted[k])
                print("p_values_two_sided_cdf: " + str(p_values_two_sided_cdf))

            except RuntimeWarning:
                print("Underflow at 'cdf' and n = " + str(n) + "!")
                exit(1)

            print()

            # Get natural logarithm of CDF for simple and twisted
            try:
                print("-ln(cdf_list_simple): " + str(-log(cdf_list_simple)))
                print("-ln(cdf_list_twisted): " + str(-log(cdf_list_twisted)))
                print("-ln(p_values_two_sided_cdf): " + str(-log(p_values_two_sided_cdf)))
            except RuntimeWarning:
                print("Underflow at '-ln' and n = " + str(n) + "!")
                exit(1)

            print()

            # Get logcdf for simple and twisted
            try:
                logcdf_list_simple = binom.logcdf(k_list, n, prob_simple)
                print("logcdf_list_simple: " + str(logcdf_list_simple))
                logsf_list_simple = binom.logsf([n - x - 1 for x in k_list], n, prob_simple)
                print("logsf_list_simple: " + str(logsf_list_simple))

                p_values_two_sided_logcdfsf = []
                for k in k_list[:int(n/2)]:
                    p_values_two_sided_logcdfsf.append(-logaddexp(logcdf_list_simple[k],logsf_list_simple[k]))
                print("p_values_two_sided_logcdfsf: " + str(p_values_two_sided_logcdfsf))

            except RuntimeWarning:
                print("Underflow at 'logcdf' and n = " + str(n) + "!")
                exit(1)

    def _pre_calculate_p_values_using_logcdf_and_logsf(self):
        """
        This function is an attempt to implement a two-sided binomial test for directionality,
        whereby the probability of success (prob_simple) does not necessarily have to be 0.5.
        The functions uses 'logcdf', 'logsf' and 'logaddexp' from numpy are used, which allow
        particularly small P-values to be calculated.
        For a given n, the P-values for k = 0, ..., n are calculated at once and saved in a
        dictionary, which is much faster than calculating the P-values one by one.
        The function appears to be working for the left tail of the distribution,
        but for the right tail (simple counts greater than n/2) a case distinction would have
        to be maade and it is not clear how.

        :return: PVAL_DICT dictionary with n as key and P-values for k = 0, ..., n as values.
        """

        # Dictionary with pre-calculated P-values to be returned
        PVAL_DICT = {}

        # List from 1 to n
        n_list = list(range(1, self._n_max + 1))

        # In one iteraation of this loop, we calculate all P-values for a given n.
        for n in n_list:

            print(n)

            # List k = 0, ..., n for which P-values are to be calculated
            k_list = list(range(0, n + 1))

            # Maybe we shoul work with two lists that we put together afterwards
            k_list_a = list(range(0, int(n/2) + 1))
            k_list_b = list(range(int(n/2)+1, n + 1))

            # Init array that will be stored under the key n
            PVAL_DICT[n] = []

            try:

                logcdf_list_simple = binom.logcdf(k_list, n, self._prob_simple)
                logsf_list_simple = binom.logsf([n - k for k in k_list], n, self._prob_simple)

                PVAL_DICT[n] = -logaddexp(logcdf_list_simple[k_list], logsf_list_simple[k_list])

                print("p_values_two_sided_logcdfsf n=" + str(n) + ": " + str(PVAL_DICT[n]))

            except RuntimeWarning:
                print("Underflow at n = " + str(n) + "!")
                exit(1)

        return PVAL_DICT