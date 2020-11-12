from typing import Dict, List, Tuple
from collections import defaultdict
from scipy.stats import binom
import numpy as np
from .diachr_util import calculate_binomial_p_value, find_indefinable_n, get_n_dict, get_n_dict_definable
from .enhanced_interaction_parser import EnhancedInteractionParser

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

        # Find smallest number of read pairs (n) that yields a significant P-value
        self._n_indef, self._pv_indef = find_indefinable_n(self._p_value_cutoff)

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

        # Output parameters
        self._print_params()

    def _print_params(self):
        print("[INFO] " + "Parameters")
        print("\t[INFO] _out_prefix: " + self._out_prefix)
        print("\t[INFO] _p_value_cutoff: " + str(self._p_value_cutoff))
        print("\t[INFO] _n_indef: " + str(self._n_indef))
        print("\t[INFO] _pv_indef: " + str(self._pv_indef))


    def simulate_interactions(self, n_max: int=500,  i_num: int=100000) -> Tuple[List, List, List]:
        """
        Simulate interactions and return lists with counts of simulated and significant interactions with n read pairs.

        Parameters
        ----------------------------
        n_max: int = 500,
            Simulate interactions with _n_indef to n_max read pairs.
        i_num: int = 100000,
            Number of simulated interactions.

        Returns
        ----------------------------
        A tuple of a lists with counts of simulated and significant interactions with n read pairs.
        """

        # List of n
        n_list = list(range(n_max + 1))

        # Lists containing numbers of simulated and significant interactions for each n
        n_sim_list = [0 for _ in range(n_max + 1)]
        n_sig_list = [0 for _ in range(n_max + 1)]

        print("[INFO] " + "Generating random numbers of simple and twisted read pairs ...")

        # Get random numbers from uniform distribution
        random_n_vec = np.random.randint(low = self._n_indef, high = n_max  + 1, size = i_num)

        # Count interactions that have the same number of read pairs for each n
        N_DICT_SIM = defaultdict(int)
        for n in random_n_vec:
            N_DICT_SIM[n] += 1

        # Transform dictionary to list (index=n)
        for i in range(0, (n_max + 1)):
            n_sim_list[i] = N_DICT_SIM[i]

        print("[INFO] " + "Counting significant interactions for each n ...")

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
                if pv <= self._p_value_cutoff/2:
                    n_sig_list[n] += 1

        self.n_max = n_max
        self.n_list = n_list
        self.n_sim_list = n_sim_list
        self.n_sig_list = n_sig_list
        return n_list, n_sim_list, n_sig_list


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

    def write_simulated_interaction_counts(self):
        """
        Get the counts of significant interactions in simulated data, and write to file
        Write a text file with significant empirical interactions for each n
        If no value for n can be found, write 0.
        Uses count_simulated_interactions()
        Writes to file "{OUTPREFIX}_sig_interactions_vs_uniform_n.tab"
        Parameters
        ----------------------------
        signum_list: List,
            a list containing counts of significant simple and twisted interactions for each n
        sim_dict: Dict,
            a dictionary with ?
        Returns
        ----------------------------
        """
        sim_tab_file_name = self._outprefix + "_sig_interactions_vs_uniform_n.tab"
        signum_list, sim_dict = self.count_simulated_interactions()
        N = len(signum_list)
        with open(sim_tab_file_name, 'wt') as fh:
            print("[INFO] " + "Writing numbers of significant simulated interactions for each n to text file ...")
            for i in range(1, N):
                fh.write(str(i) + "\t" + str(signum_list[i]) + "\t" + str(sim_dict.get(i,0)) + "\n")

    def write_significant_empirical_interactions(self, ei_file:str):
        """
        Get the counts of significant interactions in emprical data for the 
        diachromatic extended interaction file (iefile) and write the to file.
        Uses count_significant_empirical_interactions().
        Parameters
        ----------------------------
        ei_file: str,
            path to a diachromatic extended interaction file .
        n_indef: int,
            todo add definition.
        outprefix: str,
            Prefix for outfile, "{OUTPREFIX}_sig_interactions_vs_empirical_n.tab"
        """
        n_def_list, n_sig_list = self.count_significant_empirical_interactions(eifile=ei_file)
        if len(n_def_list) != len(n_sig_list):
            # should never happen
            raise ValueError("Lengths of n_def_list and n_sig_list do not match")
        N = len(n_def_list)
        emp_tab_file_name = self._out_prefix + "_sig_interactions_vs_empirical_n.tab"
        print("[INFO] " + "Writing numbers of significant empirical interactions for each n to text file ...")
        with open(emp_tab_file_name, 'wt') as fh:
            for i in range(N):
                n_def = n_def_list[i]
                n_sig = n_sig_list[i]
                fh.write(str(i) + "\t" + str(n_def) + "\t" + str(n_sig) + "\n")
