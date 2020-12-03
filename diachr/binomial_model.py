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
        #self.compare_different_ways_of_p_value_calculation(N=1000,prob_simple=0.5)

        self._pre_calculate_p_values_using_pmf_and_log()
        #self._pre_calculate_p_values_using_logcdf_and_logsf()

    def _pre_calculate_p_values_using_pmf_and_log(self):

        PVAL_DICT = {}
        n_list = list(range(1, self._n_max + 1))
        print(n_list)
        print("#####")

        for n in n_list:
            print("============================= " + str(n))
            k_list = list(range(0, n + 1))
            print(binom.pmf(k_list, n, self._prob_simple))

            for k in k_list:
                if k < n-k-1:
                    d_sum_left = 0
                    d_sum_right = 0
                    for i in range(0, k + 1):
                        print(str(i) + "\t" + str(n-i))
                        d_sum_left = d_sum_left + binom.pmf(i, n, self._prob_simple)
                        d_sum_right = d_sum_right + binom.pmf(n-i, n, self._prob_simple)
                else:
                    d_sum_left = 0
                    d_sum_right = 0
                    for i in range(k, n + 1):
                        print(str(i) + "\t" + str(n-i))
                        d_sum_left = d_sum_left + binom.pmf(i, n, self._prob_simple)
                        d_sum_right = d_sum_right + binom.pmf(n-i, n, self._prob_simple)

                print(str(k) + "\t" + str(n-k) + "\t" + str(d_sum_left) + "\t" + str(d_sum_right) + "\t" + str(d_sum_left + d_sum_right))

        return PVAL_DICT


    def _pre_calculate_p_values_using_logcdf_and_logsf(self):

        PVAL_DICT = {}

        n_list = list(range(1, self._n_max + 1))

        for n in n_list:

            print(n)
            k_list = list(range(0, n + 1))
            PVAL_DICT[n] = []

            try:
                logcdf_list_simple = binom.logcdf(k_list, n, self._prob_simple)
                logsf_list_simple = binom.logsf([n - k for k in k_list], n, self._prob_simple)
                print(k_list)
                print([n - k for k in k_list])
                print(logcdf_list_simple)
                print(logsf_list_simple)
                print(-logaddexp(logcdf_list_simple[k_list], logsf_list_simple[k_list]))

                print(n)
                for k in k_list[:int(n+1)]:
                    print(k)
                    PVAL_DICT[n].append(-logaddexp(logcdf_list_simple[k], logsf_list_simple[k]))

                print("p_values_two_sided_logcdfsf " + str(n) + ": " + str(PVAL_DICT[n]))

            except RuntimeWarning:
                print("Underflow at n = " + str(n) + "!")
                exit(1)

        return PVAL_DICT

    def compare_different_ways_of_p_value_calculation(self, N=5, prob_simple=0.5):

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

    def calculate_binomial_logsf_p_value(n_simple, n_twisted): # (natural) logsf
        try:
            if n_simple < n_twisted:
                p_value = binom.logsf(n_twisted - 1, n_simple + n_twisted, 0.5)
                return p_value
            else:
                p_value = binom.logsf(n_simple - 1, n_simple + n_twisted, 0.5)
                return p_value
        except RuntimeWarning: # Underflow: P-value is too small
            return -np.inf
            #return log(sys.float_info.min * sys.float_info.epsilon) # return natural log of smallest possible float