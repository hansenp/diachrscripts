from typing import Dict, List, Tuple
from collections import defaultdict
from scipy.stats import binom
import numpy as np
from .diachr_util import calculate_binomial_p_value, find_indefinable_n, get_n_dict, get_n_dict_definable

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

    If a Diachromatic interaction file is specified, the numbers of significant empirical interactions for each n will also
    be determined and written to file.
    """
    def __init__(
        self, 
        n_max: int=10000, 
        i_num: int=500, 
        p_value_cutoff: float=0.05, 
        out_prefix: str = "OUTPREFIX") -> None:
        """Create new BinomialInteractionModel object.

        Parameters
        ----------------------------
        n_max: int = 10000,
            Causes the class to simulate interactions with 1 to n_max read pairs.
        i_num: int = 500,
            Number of simulated interactions.
        p_value_cutoff: float=0.05, 
            P-value threshold used for the classification of directed and undirected interactions.
        diachromatic_interaction_file: str=None,
            Diachromatic interaction file.
        """
        self._n_max = n_max
        self._i_num = i_num
        self._p_value_cutoff = p_value_cutoff
        self._out_prefix = out_prefix
        # Determine indefinable cutoff for given P-value
        self._n_indef, self._pv_indef = find_indefinable_n(self._p_value_cutoff)
        # Dictionary that keeps track of already calculated P-values
        #    key - a string like 2-7
        #    value - our corresponding binomial p-value
        #    note -- use this as a global variable in this script!
        self._pval_memo = defaultdict(float)
        self._print_params()

    def _print_params(self):
        print("[INFO] " + "Input parameters")
        print("\t[INFO] --out_prefix: " + self._out_prefix)
        print("\t[INFO] --i-num: " + str(self._i_num))
        print("\t[INFO] --n-max: " + str(self._n_max))
        print("\t[INFO] --p-value-cutoff: " + str(self._p_value_cutoff))
        

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

    def count_simulated_interactions(self) -> Tuple[List, Dict]:
        """
        Get the counts of significant interactions in simulated data 
        Parameters
        ----------------------------

        Returns
        ----------------------------
        A tuple of a list with interaction counts and a dictionary with corresponding p values.
        """
        # Dictionary stores the numbers of significant interactions for each n
        N_SIG_DICT_SIM = defaultdict(int)
        # List containing counts of significant simple and twisted interactions for each n
        signum_list = [0] * (self._n_max + 1)
        print("[INFO] " + "Generating random numbers of simple and twisted read pairs ...")
        random_n_vec = np.random.randint(low = self._n_indef, high = self._n_max  + 1, size = self._i_num)
        for n in random_n_vec:
            N_SIG_DICT_SIM[n] += 1
        print("[INFO] " + "Counting significant interactions for each n ...")
        # Iterate dictionary with numbers of interactions for each read pair number n
        for n, i in N_SIG_DICT_SIM.items():
        # Generate random simple read pair counts for current n
            simple_count_list = list(binom.rvs(n, p = 0.5, size = i))
            for simple_count in simple_count_list:
                twisted_count = n - simple_count
                pv = self.binomial_p_value(simple_count=simple_count, twisted_count=twisted_count)
                # Count significant interactions for current n
                if pv <= self._p_value_cutoff:
                    signum_list[n] += 1
        return signum_list, N_SIG_DICT_SIM
        
    

    def write_simulated_interaction_counts(self, outprefix:str="OUTPREFIX"):
        """
        Get the counts of significant interactions in simulated data, and write to file
        Write a text file with significant empirical interactions for each n
        If no value for n can be found, write 0.
        Uses count_simulated_interactions()
        Parameters
        ----------------------------
        signum_list: List,
            a list containing counts of significant simple and twisted interactions for each n
        sim_dict: Dict,
            a dictionary with ?
        outprefix: str,
            Prefix for outfile, "{OUTPREFIX}_sig_interactions_vs_uniform_n.tab"
        Returns
        ----------------------------
        """
        sim_tab_file_name = outprefix + "_sig_interactions_vs_uniform_n.tab"
        signum_list, sim_dict = self.count_simulated_interactions()
        N = len(signum_list)
        with open(sim_tab_file_name, 'wt') as fh:
            print("[INFO] " + "Writing numbers of significant simulated interactions for each n to text file ...")
            for i in range(1, N):
                fh.write(str(i) + "\t" + str(signum_list[i]) + "\t" + str(sim_dict.get(i,0)) + "\n")
                

    def count_significant_empirical_interactions(self, eifile: str, min_dist=20000) -> Tuple[List, List]:
        """
        Get the counts of significant interactions in emprical data for the 
        diachromatic extended interaction file (iefile).
        Parameters
        ----------------------------
        iefile: str,
            path to a diachromatic extended interaction file .
        min_dist: int,
            minimum distance between digests in a digest pair to be counted.

        Returns
        ----------------------------
        A tuple of two lists, n_def and n_sig.
        """
        n_def_list = []
        n_sig_list = []
        # Count significant interactions in empirical data for each n
        N_SIG_DICT_EMP = get_n_dict(ei_file=eifile, status_pair_flag='ALL', min_digest_dist=min_dist, p_value_cutoff=self._p_value_cutoff)
        N_DEF_DICT_EMP = get_n_dict_definable(ei_file=eifile, status_pair_flag='ALL', min_digest_dist=min_dist, n_indef=self._n_indef)
        for n in range(1, 2000):
            n_def = N_DEF_DICT_EMP.get(n,0)
            n_sig = N_SIG_DICT_EMP.get(n,0)
            n_def_list.append(n_def)
            n_sig_list.append(n_sig)
        return n_def_list, n_sig_list


    def write_significant_empirical_interactions(self, ei_file:str, n_indef:int, outprefix: str="OUTPREFIX"):
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
       
    

     



