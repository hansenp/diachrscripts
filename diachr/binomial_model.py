from scipy.stats import binom
from numpy import log, exp, logaddexp
import warnings

class BinomialModel:
    """
    In this class, interactions are statistically evaluated in terms of directionality. The main task of this class
    is the calculation of P-values. The distribution required for this is instantiated only once for an object of this
    class and is then used for all calculations. In addition, to save time, P-values that have already been calculated
    for given pairs of simple and twisted read pair counts are saved in a dictionary to save time.
    """
    def __init__(self, probabilty_of_simple:float, n_max:int) -> None:
        self._probabilty_of_simple = probabilty_of_simple
        self._probabilty_of_twisted = 1.0 - probabilty_of_simple
        self._n_max = n_max


        self._pval_dict = {}

        print("n_max: " + str(self._n_max))
        #self._precalculate_p_values_using_pmf()
        self.compare_different_ways_of_p_value_calculation(N=10,prob_simple=0.5)

    def _pre_calculate_p_values_using_logcdf(self):

        n_list = list(range(1, self._n_max + 1))

        for n in n_list:
            pass

        cdf_list_simple = -binom.logcdf(k_list, n, prob_simple)
        print("logcdf_list_simple: " + str(cdf_list_simple))
        cdf_list_twisted = -binom.logcdf(k_list, n, prob_twisted)
        print("logcdf_list_twisted: " + str(cdf_list_twisted))

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
                logcdf_list_simple = -binom.logcdf(k_list, n, prob_simple)
                print("logcdf_list_simple: " + str(logcdf_list_simple))
                logcdf_list_twisted = -binom.logsf(k_list, n, prob_twisted)
                print("logcdf_list_twisted: " + str(logcdf_list_twisted))

                p_values_two_sided_logcdf = []
                for k in k_list[:int(n/2)]:
                    p_values_two_sided_logcdf.append(logaddexp(logcdf_list_simple[k],logcdf_list_twisted[k]))
                print("p_values_two_sided_logcdf: " + str(log(p_values_two_sided_logcdf)))

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