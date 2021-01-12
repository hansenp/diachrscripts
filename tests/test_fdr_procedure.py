from unittest import TestCase
import os
import sys

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from diachr.diachromatic_interaction_set import DiachromaticInteractionSet
from diachr.randomize_interaction_set import RandomizeInteractionSet

class TestRandomizeFDR(TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Here it is tested whether the P-value thresholds are used correctly during the FDR procedure.
        For this  purpose, we created a test file that contains 10 interactions each from 200 consecutive P-value
        intervals of equal size.
        We read this file into an object of class 'DiachromaticInteractionSet' and pass this object on to an object
        of class 'RandomizeInteractionSet'.
        The class 'RandomizeInteractionSet' contains the function in which the FDR procedure is implemented.
        We call this function with the same arguments for the maximum P-value threshold and step size that we
        used when creating the test file.
        Therefore, we expect the number of significant interactions to increase by exactly 1 as we go from one to the
        next larger P-value threshold. This is what is tested here.
        """

        cls.interaction_set = DiachromaticInteractionSet()
        cls.interaction_set.parse_file("data/test_03/diachromatic_interaction_file_fdr_1.tsv")
        cls.randomize_fdr = RandomizeInteractionSet(interaction_set=cls.interaction_set)

    def test_interactions_numbers_at_different_pval_thresholds(self):

        # The function for the FDR procedure returns a table containing the nnumbers that are being tested here
        fdr_info_dict = self.randomize_fdr.get_pval_tresh_at_chosen_fdr_tresh(
            chosen_fdr_thresh=0.05,
            pval_thresh_max=0.05,
            pval_thresh_step_size=0.00025,
            verbose=False)

        # We test the number of significant interactions for all P-value thresholds
        sig_num_o_expected = 10
        for sig_num_o in fdr_info_dict['RESULTS_TABLE']['SIG_NUM_O']:
            self.assertEqual(sig_num_o_expected,sig_num_o, msg='Number of significant interactions did not increase by 10!')
            sig_num_o_expected += 10
