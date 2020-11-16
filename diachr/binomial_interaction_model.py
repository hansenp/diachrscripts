from typing import Dict, List, Tuple
from collections import defaultdict
from scipy.stats import binom
import numpy as np
from .diachr_util import calculate_binomial_p_value
from .enhanced_interaction_parser import EnhancedInteractionParser
import scipy, scipy.stats, numpy
import matplotlib.pyplot as plt

class BinomialInteractionModel:
    """
    This class is to investigate the relationship between signal strength, measured as the total number of read pairs,
    and the binomial model for directionality of interactions.

    For this purpose, a specified number of interactions is simulated, whereby the the total numbers of read pairs are drawn
    from a uniform distribution. For individual interactions with n read pairs, the number of simple read pairs is drawn
    from a binomial distribution with p = 0.5 (null model) and the number of twisted read pairs is set to n minus the number
    of simple read pairs. Subsequently, the simulated interactions are evaluated for statistical significance given the null
    model and a specified significance threshold. Finally, a plot for the number of significant simulated interactions for
    each n is generated and, in addition, the corresponding numbers are written to a tab separated text file.

    If a enhanced interaction (EI) file is specified, the numbers of significant empirical interactions for each n will also
    be determined and written to file.


    """

    def __init__(
        self,
        p_value_cutoff: float=0.05,
        out_prefix: str = "OUT_PREFIX") -> None:

        # Threshold used to define significant interactions
        self._p_value_cutoff = p_value_cutoff

        # Prefix for the names of generated files
        self._out_prefix = out_prefix

        # Dictionary that keeps track of already calculated P-values (saves time)
        #    key - a string like 2-7
        #    value - our corresponding binomial p-value
        self._pval_memo = defaultdict(float)

        # Maximum n for simulated interactions
        self.n_max = None

        # List containing 1, ..., n_max
        self.n_list = None

        # List containing numbers of simulated interactions for each n
        self.n_sim_list = None

        # List containing numbers of simulated significant interactions for each n
        self.n_sig_list = None

        # Define colors for directed, undirected reference and all undirected interactions
        self._di_color = (255/255, 163/255, 0/255, 1)
        self._uir_color = (171/255,215/255,230/255,1)
        self._ui_color = (210/255,210/255,210/255,1)

        # Output parameters
        #self._print_params()

    def _print_params(self):
        print("[INFO] " + "Parameters")
        print("\t[INFO] _out_prefix: " + self._out_prefix)
        print("\t[INFO] _p_value_cutoff: " + str(self._p_value_cutoff))


    def simulate_interactions(self, n_max: int=500,  i_num: int=100000, pvt: float=0.05) -> Tuple[List, List, List]:
        """
        Simulate interactions and return lists with counts of simulated and significant interactions with n read pairs.

        Parameters
        ----------------------------
        n_max: int = 500,
            Simulate interactions with 0 to n_max read pairs.
        i_num: int = 100000,
            Number of simulated interactions.
        pvt: float = 0.05,
            P-value threshold.

        Returns
        ----------------------------
        A tuple of a lists with counts of simulated and significant interactions with n read pairs.
        """

        print("[INFO] " + "Running simulation ...")

        # List of n
        n_list = list(range(n_max + 1))

        # Lists containing numbers of simulated and significant interactions for each n
        n_sim_list = [int(i_num / (n_max + 1)) for _ in range(n_max + 1)]
        n_sig_list = [0 for _ in range(n_max + 1)]

        print("\t[INFO] " + "Counting significant interactions for each n ...")

        # Iterate dictionary with numbers of interactions for each read pair number n
        for n in range(0, (n_max + 1)):
            if n_sim_list[n]==0:
                continue

            # Generate random simple read pair counts for current n
            simple_count_list = list(binom.rvs(n, p = 0.5, size = n_sim_list[n]))

            # Determine P-values and count significant interactions
            for simple_count in simple_count_list:
                twisted_count = n - simple_count
                pv = self.binomial_p_value(simple_count=simple_count, twisted_count=twisted_count)
                # Our test of directionality is two-sided (either simple or twisted).
                # Therefore, we have to divide the threshold by 2.
                if pv <= pvt/2.0:
                    n_sig_list[n] += 1

        print("[INFO] " + "...done.")

        return n_list, n_sim_list, n_sig_list

    def simulate_interactions_plot(self, n_max: int=None,  i_num: int=None, pvt: float=None,
                                   n_list: List=None, n_sim_list: List=None,  n_sig_list: List=None, CREATE_PDF=False):
        fig, ax1 = plt.subplots()
        plt.title("n_max: " + str(n_max) + ", i_num: " + str(i_num) + ", pvt: " + str(pvt))
        color = 'gray'
        ax1.set_ylim(0, max(n_sim_list) + int(max(n_sim_list) / 25))
        ax1.set_xlabel('Number of read pairs per interaction (n)')
        ax1.set_ylabel('Simulated interactions', color=color)
        ax1.scatter(n_list, n_sim_list, color=color, alpha=0.5)
        ax1.tick_params(axis='y', labelcolor=color)
        ax2 = ax1.twinx()
        color = self._di_color
        ax2.set_ylabel('Significant interactions', color=color)
        ax2.plot(n_list, n_sig_list, color=color, linewidth=0.5)
        ax2.hlines(pvt * (i_num / (n_max + 1)), 0, n_max, color='red', linestyle='dashed')
        ax2.tick_params(axis='y', labelcolor=color)
        fig.tight_layout()
        if CREATE_PDF:
            plt.savefig(self._out_prefix + "_simulated_interactions_n_max_" + str(n_max)  + "_i_num_" + str(i_num)  + "_pvt_" + str(pvt) + ".pdf", format="pdf")


    def analyze_N_binomial_distributions_with_fixed_p_thresh(self, N=40, pvt=0.05, CREATE_DIST_PLOTS=True, CREATE_PDF=False):

        # Intit lists for the final plot
        n_list = []
        d_sum_prev_list = []
        k_pvt_list = []

        # Create binomial distributions for n = 1, ..., N
        for n in range(1, N + 1):

            # Get list of densities for current n
            x = list(scipy.arange(0, n + 1, 1.0))
            pmf = scipy.stats.binom.pmf(x, n, 0.5)

            # Find largest k for which the sum of densities from 0 up to k is below the P-value threshold
            d_sum = 0.0
            d_sum_prev = 0.0
            k_pvt = 0
            for density in pmf:
                d_sum = d_sum + density
                if pvt / 2 < d_sum:
                    break
                else:
                    k_pvt += 1
                    d_sum_prev = d_sum

            # Collect determined values for the final plot
            n_list.append(n)
            d_sum_prev_list.append(d_sum_prev)
            k_pvt_list.append(k_pvt)

            # Create distribution plot for current n
            if CREATE_DIST_PLOTS:
                if d_sum_prev > 0.0:
                    color = 'blue'
                else:
                    color = 'red'
                plt.title('P_VAL_THRESH=' + str(pvt) + \
                          ', n=' + str(n) + \
                          ', $k_{max}$=' + str(k_pvt) + \
                          ', P-val=' + "{:.4f}".format(d_sum_prev))
                plt.xlabel('n')
                plt.ylabel('Density')
                plt.xlim(-1, N)
                plt.ylim(0, 0.5)
                plt.bar(x, pmf, color=color)
                plt.vlines(0, 0, N)
                plt.vlines(n, 0, N)
                plt.axvspan(0, k_pvt, facecolor='b', alpha=0.2)
                plt.axvspan(n - k_pvt, n, facecolor='b', alpha=0.2)
                plt.show()

        # Calculate trendline
        spl = scipy.interpolate.UnivariateSpline(n_list, d_sum_prev_list)
        spl.set_smoothing_factor(0.5)

        # Create final plot with two y-axes, one for d_sum and one for k_pvt
        fig, ax1 = plt.subplots()
        plt.title("N: " + str(N) + ", pvt: " + str(pvt))
        color = 'cornflowerblue'
        ax1.set_xlabel('n')
        ax1.set_ylabel('P-value', color=color)
        ax1.plot(n_list, d_sum_prev_list, color=color)
        ax1.plot(n_list, spl(n_list), color='blue')
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.hlines(pvt / 2, 0, N, color='red', linestyle='dashed')

        ax2 = ax1.twinx()

        color = 'black'
        ax2.set_ylabel('$k_{max}$', color=color)
        ax2.scatter(n_list, k_pvt_list, color=color)
        ax2.tick_params(axis='y', labelcolor=color)

        fig.tight_layout()
        if CREATE_PDF:
            plt.savefig(self._out_prefix + "_N_binom_dist_with_fixed_pvt_N_" + str(N) + "_pvt_" + str(pvt) + ".pdf", format="pdf")

    def count_di_uir_and_ui_for_each_n(self, ei_file: str, n_max=2000) -> Tuple[List, List, List, List]:
        """
        Get the counts of DI, UIR and UI interactions for each read pair number.

        Parameters
        ----------------------------
        iefile: str,
            path to a enhanced interaction file.
        n_max: int,
            maximum read pair number for which the interaction numbers are determined.

        Returns
        ----------------------------
        n_list: List,
            contains the index which corresponds to read pair numbers.
        n_di_list: List,
            contains the DI counts for each read pair number.
        n_uir_list: List,
            contains the UIR counts for each read pair number.
        n_ui_list: List,
            contains the UI counts for each read pair number.
        """

        # Count significant interactions in empirical data for each n
        n_dict_di = {}  # Dictionary that stores the numbers of DI interactions with n read pairs
        n_dict_uir = {}  # Dictionary that stores the numbers of UIR interactions with n read pair
        n_dict_ui = {}  # Dictionary that stores the numbers of UI interactions with n read pairs

        # Get list of EI objects
        parser = EnhancedInteractionParser(ei_file)
        ei_list = parser.parse()

        # Iterate list with EI objects
        n_progress = 0
        for ei in ei_list:
            n_progress += 1
            if n_progress % 1000000 == 0:
                print("\t[INFO] Processed " + str(n_progress) + " interactions ...")

            n  = ei.rp_total
            if ei._interaction_category == 'DI':
                if n in n_dict_di:
                    n_dict_di[n] += 1
                else:
                    n_dict_di[n] = 1
            elif ei._interaction_category == 'UIR':
                if n in n_dict_uir:
                    n_dict_uir[n] += 1
                else:
                    n_dict_uir[n] = 1
            elif ei._interaction_category == 'UI':
                if n in n_dict_ui:
                    n_dict_ui[n] += 1
                else:
                    n_dict_ui[n] = 1
            else:
                print("[ERROR] Invalid interaction category tag: " + ei._interaction_category)

        # List of n
        n_list = list(range(n_max + 1))

        # Lists of interaction counts for n = 0, ..., n_max
        n_di_list = []
        n_uir_list = []
        n_ui_list = []

        # Get ordered list of interaction counts for different n
        for n in n_list:
            n_di = n_dict_di.get(n,0)
            n_di_list.append(n_di)
            n_uir = n_dict_uir.get(n,0)
            n_uir_list.append(n_uir)
            n_ui = n_dict_ui.get(n, 0)
            n_ui_list.append(n_ui)

        return n_list, n_di_list, n_uir_list, n_ui_list

    def count_di_uir_and_ui_for_each_n_plot(self, n_list: List=None, n_di_list: List=None, n_uir_list: List=None, n_ui_list: List=None):

        # Get interaction sums and relative frequencies
        max_rel = 0
        di_rel = []
        di_sum = sum(n_di_list)
        for i in n_di_list:
            di_rel.append(i / di_sum)
            if max_rel < i / di_sum:
                max_rel = i / di_sum
        uir_rel = []
        uir_sum = sum(n_uir_list)
        for i in n_uir_list:
            uir_rel.append(i / uir_sum)
            if max_rel < i / uir_sum:
                max_rel = i / uir_sum
        ui_rel = []
        ui_sum = sum(n_ui_list)
        for i in n_ui_list:
            ui_rel.append(i / ui_sum)
            if max_rel < i / ui_sum:
                max_rel = i / ui_sum

        # Plot UI vs. DI
        plt.title('UI vs. DI')
        plt.xlabel('Number of read pairs per interaction (n)')
        plt.ylabel('Relative frequency of interactions')
        plt.xlim(0, 400)
        plt.ylim(0, max_rel)
        plt.scatter(n_list, ui_rel, alpha=1, color=self._ui_color,
                    label='UI: ' + str("{:,d}".format(ui_sum)))
        plt.scatter(n_list, di_rel, alpha=0.5, color=self._di_color,
                    label='DI: ' + str("{:,d}".format(di_sum)))
        plt.legend()
        plt.savefig(self._out_prefix + "_emp_ui_vs_di" + ".pdf", format="pdf")
        plt.close()

        # Plot UI vs. DI
        plt.title('UI vs. UIR')
        plt.xlabel('Number of read pairs per interaction (n)')
        plt.ylabel('Relative frequency of interactions')
        plt.xlim(0, 400)
        plt.ylim(0, max_rel)
        plt.scatter(n_list, ui_rel, alpha=1, color=self._ui_color,
                    label='UI: ' + str("{:,d}".format(ui_sum)))
        plt.scatter(n_list, uir_rel, alpha=0.5, color=self._di_color,
                    label='uir: ' + str("{:,d}".format(uir_sum)))
        plt.legend()
        plt.savefig(self._out_prefix + "_emp_ui_vs_uir" + ".pdf", format="pdf")
        plt.close()

    def binomial_p_value(self, simple_count: int, twisted_count: int):
        """
        Locally defined method for the calculation of the binomial P-value that uses a dictionary that keeps track of
        P-values that already have been calculated.

        Parameters
        ----------------------------
        simple_count: int,
            Number of simple read pairs.
        twisted_count: int,
            Number of twisted read pairs.

        Returns
        ----------------------------
        Binomial P-value.
        """
        # Create key from simple and twisted read pair counts
        key = "{}-{}".format(simple_count, twisted_count)
        # Check whether a P-value for this combination of simple and twisted counts has been calculated already
        if key in self._pval_memo :
            return self._pval_memo [key]
        else:
            # Calculate P-value and add to dictionary
            p_value = calculate_binomial_p_value(simple_count, twisted_count)
            self._pval_memo [key] = p_value
            return p_value
