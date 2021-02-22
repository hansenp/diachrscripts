from typing import List, Tuple
from collections import defaultdict
from scipy import arange, interpolate
from scipy.stats import binom
import matplotlib.pyplot as plt
from diachr.binomial_model import BinomialModel
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet


class BinomialInteractionModel:
    """
    This class is to review our implementation of the P-value calculation and to investigate the consequences of the
    different null distributions used for the test due to different total read pair numbers of interactions.

    A more detailed description and application examples can be found in the following Jupyter notebook:

           jupyter_notebooks/binomialModel.ipynb

    In addition, this class is used in the following script to create all plots that are also created in the
    Jupyter notebook:

           00_explore_binomial_model.py
    """

    def __init__(
            self,
            p_value_cutoff: float = 0.05,
            out_prefix: str = "OUT_PREFIX") -> None:

        # Threshold used to define significant interactions
        self._p_value_cutoff = p_value_cutoff

        # Prefix for the names of generated files
        self._out_prefix = out_prefix

        # Dictionary that keeps track of already calculated P-values (saves time)
        #    key - a string like 2-7
        #    value - our corresponding binomial p-value
        self._pval_memo = defaultdict(float)

        # Object for calculating P-values
        self._p_values = BinomialModel()

        # Maximum n for simulated interactions
        self.n_max = None

        # List containing 1, ..., n_max
        self.n_list = None

        # List containing numbers of simulated interactions for each n
        self.n_sim_list = None

        # List containing numbers of simulated significant interactions for each n
        self.n_sig_list = None

        # Define colors for directed, undirected reference and all undirected interactions
        self._di_color = (255 / 255, 163 / 255, 0 / 255, 1)
        self._uir_color = (171 / 255, 215 / 255, 230 / 255, 1)
        self._ui_color = (210 / 255, 210 / 255, 210 / 255, 1)

    def simulate_interactions(self, n_max: int = 500, i_num: int = 100000, pvt: float = 0.05) -> Tuple[
        List, List, List]:
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
            if n_sim_list[n] == 0:
                continue

            # Generate random simple read pair counts for current n
            simple_count_list = list(binom.rvs(n, p=0.5, size=n_sim_list[n]))

            # Determine P-values and count significant interactions
            for simple_count in simple_count_list:
                twisted_count = n - simple_count
                pv = self._p_values.get_binomial_p_value(simple_count, twisted_count)
                if pv <= pvt:
                    n_sig_list[n] += 1

        print("[INFO] " + "...done.")

        return n_list, n_sim_list, n_sig_list

    def simulate_interactions_plot(self, n_max: int = None, i_num: int = None, pvt: float = None,
                                   n_list: List = None, n_sim_list: List = None, n_sig_list: List = None,
                                   CREATE_PDF=False):
        """
        Plot all simulated and simulated significant interactions.

        Parameters
        ----------------------------
        n_max: int = None,
            Shown in the title and file name and used to calculate the expected number of significant interactions.
        i_num: int = None,
            Shown in the title and file name and used to calculate the expected number of significant interactions.
        pvt: float = None,
            Shown in the title and file name and used to calculate the expected number of significant interactions.
        n_list: List = None,
            List of indices (x-axis).
        n_sim_list: List = None,
            List of simulated interaction numbers for n=1,...,n_max (left y-axis).
        n_sig_list: List = None,
            List of simulated significant interaction numbers for n=1,...,n_max (right y-axis).
        CREATE_PDF: bool = None,
            If 'True', a PDF will be created.

        Returns
        ----------------------------
        Nothing. Only the plot is generated.
        """
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
            plt.savefig(self._out_prefix + "_simulated_interactions_n_max_" + str(n_max) + "_i_num_" + str(
                i_num) + "_pvt_" + str(pvt) + ".pdf", format="pdf")

    def analyze_N_binomial_distributions_with_fixed_p_thresh(self, N=40, pvt=0.05, CREATE_DIST_PLOTS=True,
                                                             CREATE_PDF=False):
        """
        Plot PMF of binomial distributions for n=1,...,N and create a combined plot with maximal P-values and
        associated k.

        Parameters
        ----------------------------
        N: int = 40,
            Number of PMF that will be analyzed.
        pvt: float = 0.05,
            P-value threshold.
        CREATE_DIST_PLOTS: bool = None,
            If 'True', the PMF will be plotted for each n=1,...,N.
        CREATE_PDF: bool = None,
            If 'True', a PDF with the combined plot will be created.

        Returns
        ----------------------------
        Nothing. Only the plots are generated.
        """
        n_list = []
        d_sum_prev_list = []
        k_pvt_list = []

        # Create binomial distributions for n = 1, ..., N
        for n in range(1, N + 1):

            # Get list of densities for current n
            x = list(arange(0, n + 1, 1.0))
            pmf = binom.pmf(x, n, 0.5)

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
                plt.close()

        # Calculate trendline
        spl = interpolate.UnivariateSpline(n_list, d_sum_prev_list)
        spl.set_smoothing_factor(0.5)

        # Set ticks for x-axis
        xticks = list(range(0, N + 1, int(N / 10)))

        # Set ticks for right y-axis (k_max)
        max_k = max(k_pvt_list)
        if max_k < 11:
            yticks_k = list(range(0, max_k + 1, 1))
        elif max_k < 21:
            yticks_k = list(range(0, max_k + 1, 2))
        elif max_k < 51:
            yticks_k = list(range(0, max_k + 1, 5))
        elif max_k < 101:
            yticks_k = list(range(0, max_k + 1, 10))
        else:
            yticks_k = list(range(0, max_k + 1, 20))

        # Create final plot with two y-axes, one for d_sum and one for k_pvt
        fig, ax1 = plt.subplots()
        plt.title("N: " + str(N) + ", pvt: " + str(pvt))
        plt.xticks(xticks)
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
        ax2.set_yticks(yticks_k)

        fig.tight_layout()

        if CREATE_PDF:
            plt.savefig(self._out_prefix + "_N_binom_dist_with_fixed_pvt_N_" + str(N) + "_pvt_" + str(pvt) + ".pdf",
                        format="pdf")

        plt.show()
        plt.close()

    def count_di_uir_and_ui_for_each_n(self, di11_file: str, n_max=2000) -> Tuple[List, List, List, List]:
        """
        Get the counts of DI, UIR and UI interactions for each read pair number.

        Parameters
        ----------------------------
        iefile: str,
            Path to a enhanced interaction file.
        n_max: int,
            Maximum read pair number for which the interaction numbers are determined.

        Returns
        ----------------------------
        n_list: List,
            Contains the index which corresponds to read pair numbers (n).
        n_di_list: List,
            Contains the DI counts for each read pair number.
        n_uir_list: List,
            Contains the UIR counts for each read pair number.
        n_ui_list: List,
            Contains the UI counts for each read pair number.
        """

        # Count significant interactions in empirical data for each n
        n_dict_di = {}  # Dictionary that stores the numbers of DI interactions with n read pairs
        n_dict_uir = {}  # Dictionary that stores the numbers of UIR interactions with n read pair
        n_dict_ui = {}  # Dictionary that stores the numbers of UI interactions with n read pairs

        # Load interactions
        interaction_set = DiachromaticInteractionSet()
        interaction_set.parse_file(di11_file, verbose=True)

        # Iterate interaction set
        n_progress = 0
        for di_11_inter in interaction_set._inter_dict.values():
            n_progress += 1
            if n_progress % 1000000 == 0:
                print("\t[INFO] Processed " + str(n_progress) + " interactions ...")

            n = di_11_inter.rp_total
            if di_11_inter.get_category() == 'DI':
                if n in n_dict_di:
                    n_dict_di[n] += 1
                else:
                    n_dict_di[n] = 1
            elif di_11_inter.get_category() == 'UIR':
                if n in n_dict_uir:
                    n_dict_uir[n] += 1
                else:
                    n_dict_uir[n] = 1
            elif di_11_inter.get_category() == 'UI':
                if n in n_dict_ui:
                    n_dict_ui[n] += 1
                else:
                    n_dict_ui[n] = 1
            else:
                print("[ERROR] Invalid interaction category tag: " + d_inter.get_category())

        # List of n
        n_list = list(range(n_max + 1))

        # Lists of interaction counts for n = 0, ..., n_max
        n_di_list = []
        n_uir_list = []
        n_ui_list = []

        # Get ordered list of interaction counts for different n
        for n in n_list:
            n_di = n_dict_di.get(n, 0)
            n_di_list.append(n_di)
            n_uir = n_dict_uir.get(n, 0)
            n_uir_list.append(n_uir)
            n_ui = n_dict_ui.get(n, 0)
            n_ui_list.append(n_ui)

        return n_list, n_di_list, n_uir_list, n_ui_list

    def count_di_uir_and_ui_for_each_n_plot(self, n_list: List = None, n_di_list: List = None, n_uir_list: List = None,
                                            n_ui_list: List = None, x_max: int = None, y_max: int = None,
                                            l_wd: float = 0.0):
        """
        Create plots for empirical data.

        Parameters
        ----------------------------
        n_list: List = None,
            List of indices (x-axis).
        n_di_list: List = None,
            List of DI counts.
        n_uir_list: List = None,
            List of UIR counts.
        n_ui_list: List = None,
            List of UI counts.
        x_max: int = None,
            Limit for x-axis.
        y_max: int = None,
            Limit for y-axis.
        l_wd: float = 0.0,
            If 0.0<l_wd, neighboring points in the scatterplots are connected by lines.

        Returns
        ----------------------------
        Nothing. Only the plots are generated.
        """

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

        # Set common limit of y-axes
        if y_max == None:
            y_max = max_rel

        # Set ticks for x-axes
        xticks = list(range(0, x_max + 1, int(x_max / 10)))

        # Plot UI vs. DI
        plt.title('UI vs. DI')
        plt.xlabel('Number of read pairs per interaction (n)')
        plt.ylabel('Relative frequency of interactions')
        plt.xlim(0, x_max)
        plt.ylim(0, y_max)
        plt.xticks(xticks)

        plt.scatter(n_list, ui_rel, alpha=1, color=self._ui_color,
                    label='UI: ' + str("{:,d}".format(ui_sum)))

        plt.plot(n_list, di_rel, color=self._di_color, linewidth=l_wd)
        plt.scatter(n_list, di_rel, alpha=0.5, color=self._di_color,
                    label='DI: ' + str("{:,d}".format(di_sum)))
        plt.legend()
        plt.savefig(self._out_prefix + "_x_max_" + str(x_max) + "_emp_ui_vs_di" + ".pdf", format="pdf")
        plt.show()
        plt.close()

        # Plot UI vs. UIR
        plt.title('UI vs. UIR')
        plt.xlabel('Number of read pairs per interaction (n)')
        plt.ylabel('Relative frequency of interactions')
        plt.xlim(0, x_max)
        plt.ylim(0, y_max)
        plt.xticks(xticks)

        plt.scatter(n_list, ui_rel, alpha=1, color=self._ui_color,
                    label='UI: ' + str("{:,d}".format(ui_sum)))

        plt.plot(n_list, uir_rel, color=self._uir_color, linewidth=l_wd)
        plt.scatter(n_list, uir_rel, alpha=0.5, color=self._uir_color,
                    label='UIR: ' + str("{:,d}".format(uir_sum)))
        plt.legend()
        plt.savefig(self._out_prefix + "_x_max_" + str(x_max) + "_emp_ui_vs_uir" + ".pdf", format="pdf")
        plt.show()
        plt.close()

        # Plot DI vs. UIR
        plt.title('DI vs. UIR')
        plt.xlabel('Number of read pairs per interaction (n)')
        plt.ylabel('Relative frequency of interactions')
        plt.xlim(0, x_max)
        plt.ylim(0, y_max)
        plt.xticks(xticks)

        plt.plot(n_list, di_rel, color=self._di_color, linewidth=l_wd)
        plt.scatter(n_list, di_rel, alpha=1, color=self._di_color,
                    label='DI: ' + str("{:,d}".format(di_sum)))

        plt.plot(n_list, uir_rel, color=self._uir_color, linewidth=l_wd)
        plt.scatter(n_list, uir_rel, alpha=0.5, color=self._uir_color,
                    label='UIR: ' + str("{:,d}".format(uir_sum)))
        plt.legend()
        plt.savefig(self._out_prefix + "_x_max_" + str(x_max) + "_emp_di_vs_uir" + ".pdf", format="pdf")
        plt.show()
        plt.close()
