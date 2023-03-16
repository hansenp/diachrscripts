from scipy.stats import binom
from numpy import inf, log, logaddexp
import time


class BinomialModel:
    """
    In this class, interactions are statistically evaluated with respect to imbalances in their read pair counts. To
    save time, P-values for given pairs of counts (the two highest vs. the two others) are stored in a dictionary so
    they only need to be computed once.
    """

    def __init__(self, two_sided: bool = True) -> None:

        # Perform two-sided or one-sided tests
        self._two_sided = two_sided

        # Dictionary for already calculated P-values
        self._pval_dict = {}

    def _get_pval_dict_size(self):
        print("Size of dictionary for P-values: " + str(len(self._pval_dict)))
        return len(self._pval_dict)

    def get_binomial_log10_p_value(self, n_simple, n_twisted):
        """
        Given a pair of counts, this function returns the negative decadic logarithm of the binomial P-value. If an
        identical pair of counts has been seen before, then the P-value is retrieved from the dictionary. Otherwise,
        the P-value is calculated, added to the dictionary, and returned.

        :return: Negative decadic logarithm of the binomial P-value
        """

        # Create key for pair of counts
        key = "{}-{}".format(n_simple, n_twisted)

        if key not in self._pval_dict:
            self._pval_dict[key] = self._calculate_binomial_log10_p_value(n_simple, n_twisted)

        return self._pval_dict[key]

    def get_binomial_p_value(self, n_simple, n_twisted):
        """
        :param n_simple: Number of simple read pairs (or the two highest counts)
        :param n_twisted: Number of simple read pairs (or the two counts that are not the highest)
        :return: P-value
        """
        return 10 ** -self.get_binomial_log10_p_value(n_simple, n_twisted)

    def _calculate_binomial_log10_p_value(self, n_simple, n_twisted):
        """
        For a given pair of counts, this function calculates the negative decadic logarithm of the P-value from a one-
        or two-sided binomial test. For P-values that are too small to be represented, 'inf' is returned.

        :param n_simple: Number of simple read pairs (or the two highest counts)
        :param n_twisted: Number of twisted read pairs (or the two counts that are not the highest)
        :return: The negative of the decadic logarithm of the binomial P-value
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
        This function finds the smallest total number of read pairs that gives a significant P-value at a chosen
        threshold (classifiable interactions).

        :param p_val_thresh: P-value threshold
        :param verbose: If true, messages about progress will be written to the screen
        :return: (Smallest total number of read pairs, Associated P-value):
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
