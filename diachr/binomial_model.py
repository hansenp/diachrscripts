from scipy.stats import binom
from numpy import exp, inf, log, logaddexp
import time


class BinomialModel:
    """
    In this class, interactions are statistically evaluated in terms of directionality. The main task of this class
    is the calculation and storage of P-values. To save time, P-values that have already been calculated for given
    pairs of simple and twisted read pair counts are saved in a dictionary to save time.

    There is a function under development further down in this class. This function is fast and could therefore be
    used to calculate all P-values in advance. In addition, this function does not require a probability of success
    of 0.5, but any probabilities can be selected.
    """

    def __init__(self, two_sided: bool = True) -> None:

        # Perform two-sided or one-sided tests
        self._two_sided = two_sided

        # Background frequency of simple read pairs
        self._prob_simple = 0.5

        # Dictionary for already calculated P-values
        self._pval_dict = {}

    def _get_pval_dict_size(self):
        print("Size of dictionary for P-values: " + str(len(self._pval_dict)))
        return len(self._pval_dict)

    def get_binomial_log10_p_value(self, n_simple, n_twisted):
        """
        This function returns already calculated P-values from the dictionary or,
        if the P-value has not yet been calculated, a function is called to
        calculate the P-value. This P-value is then added to the dictionary and returned.

        :return: The negative of the decadic logarithm of the P-value
        """

        # Create key from simple and twisted read pair counts
        key = "{}-{}".format(n_simple, n_twisted)

        if key not in self._pval_dict:
            self._pval_dict[key] = self._calculate_binomial_log10_p_value(n_simple, n_twisted)

        return self._pval_dict[key]

    def get_binomial_p_value(self, n_simple, n_twisted):
        """
        :param n_simple: Number of simple read pairs
        :param n_twisted: Number of simple read pairs
        :return: P-value of a two-sided test for directionality
        """
        return 10 ** -self.get_binomial_log10_p_value(n_simple, n_twisted)

    def _calculate_binomial_log10_p_value(self, n_simple, n_twisted):
        """
        This function calculates the decadic logarithm of a P-value from a two-sided test and returns the negative
        of it. We assume a background frequency of 0.5 for simple and twisted read pairs. Therefore, the distribution
        is symmetrical and we return twice the P-value from a one-sided test as the P-value for a two-sided test.
        For P-values that are too small to be represented, we return -inf.

        :param n_simple: Number of simple read pairs
        :param n_twisted: Number of twisted read pairs
        :return: The negative of the decadic logarithm of the P-value
        """
        try:
            if n_simple == n_twisted:
                return 0  # All other cases are more extreme, so the P-value must be 1.
            if n_simple < n_twisted:
                ln_p_value = binom.logcdf(n_simple, n_simple + n_twisted, 0.5)
            else:
                ln_p_value = binom.logcdf(n_twisted, n_simple + n_twisted, 0.5)

            if self._two_sided:
                return -logaddexp(ln_p_value, ln_p_value) / log(10) # logcdf returns the natural logarithm
            else:
                return -ln_p_value / log(10)

        except RuntimeWarning:  # Underflow: P-value is too small
            return inf
            # or return log of smallest possible float
            # return log(sys.float_info.min * sys.float_info.epsilon)/log(10)

    def find_smallest_significant_n(self, p_val_thresh: float = 0.01, verbose: bool = False):
        """
        This function finds the smallest number of read pairs that gives a significant P-value at a chosen threshold.
        A tuple consisting of this smallest number and the associated largest significant P-value is returned.

        :param p_val_thresh: P-value threshold
        :param verbose: If true, messages about progress will be written to the screen
        :return: (Smallest number of read pairs, Associated P-value):
        """

        if verbose:
            print("[INFO] Looking for smallest number of read pairs with a significant P-value at the chosen "
                  "threshold of " + str(p_val_thresh) + ".")

        for rp_num in range(1, 1000):
            p_val = self.get_binomial_p_value(rp_num, 0)
            if p_val < p_val_thresh:
                if verbose:
                    print("\t[INFO] Smallest number of read pairs: " + str(rp_num) + " read pairs (" + str(p_val) + ")")
                    print("[INFO] ... done.")
                return rp_num, p_val

    def measure_time_savings_due_to_pval_dict(self):
        """
        This function determines time saved by using the dictionary for P-values.
        For this purpose, the time is measured which is required to to calculate approx. 500,000 P-values.
        Then the time is measured that is required to query these P-values from the dictionary.
        """
        start_time = time.time()
        n_max = 1000
        for n in range(1, n_max + 1):
            for k in range(0, n + 1):
                self.get_binomial_log10_p_value(k, n - k)
        end_time = time.time()
        pval_num = self._get_pval_dict_size()

        print("It took " + "{:.2f}".format(end_time - start_time) + " seconds to calculate " + str(
            pval_num) + " P-values.")

        start_time = time.time()
        n_max = 1000
        for n in range(1, n_max + 1):
            for k in range(0, n + 1):
                self.get_binomial_log10_p_value(k, n - k)
        end_time = time.time()
        pval_num = self._get_pval_dict_size()

        print("It took " + "{:.2f}".format(end_time - start_time) + " seconds to get " + str(
            pval_num) + " P-values from the dictionary.")
