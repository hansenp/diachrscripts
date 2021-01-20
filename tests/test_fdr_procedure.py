from unittest import TestCase
import os
import sys
import warnings

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from diachr.diachromatic_interaction_set import DiachromaticInteractionSet
from diachr.randomize_interaction_set import RandomizeInteractionSet


class TestRandomizeFDR(TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Here it is tested whether the P-value thresholds are used correctly during the FDR procedure and whether the
        results for two test files with 1,000 and 64,000 interactions. More information about the creation of the
        test files can be found in the Jupyter notebook on the FDR procedure.
        """

        # Mute warnings
        warnings.simplefilter('ignore')

        # Prepare test data to test correct use of P-value thresholds
        cls.interaction_set_1 = DiachromaticInteractionSet()
        cls.interaction_set_1.parse_file("data/test_03/diachromatic_interaction_file_fdr_1.tsv.gz")
        cls.randomize_1 = RandomizeInteractionSet(interaction_set=cls.interaction_set_1)

        # Prepare test data to detect changes in the results for the test file with 1,000 interactions
        cls.interaction_set_1000 = DiachromaticInteractionSet()
        cls.interaction_set_1000.parse_file("data/test_03/diachromatic_interaction_file_fdr_top_1000.tsv.gz")
        cls.randomize_1000 = RandomizeInteractionSet(interaction_set=cls.interaction_set_1000)

        # Prepare test data to detect changes in the results for the test file with 64,000 interactions
        cls.interaction_set_64000 = DiachromaticInteractionSet()
        cls.interaction_set_64000.parse_file("data/test_03/diachromatic_interaction_file_fdr_top_64000.tsv.gz")
        cls.randomize_64000 = RandomizeInteractionSet(interaction_set=cls.interaction_set_64000)

    def test_interactions_numbers_at_different_pval_thresholds(self):
        """
        Here it is tested whether the P-value thresholds are used correctly during the FDR procedure.
        """

        # The function for the FDR procedure returns a dictionary containing the numbers that are being tested here
        fdr_info_dict = self.randomize_1.get_pval_thresh_at_chosen_fdr_thresh(
            chosen_fdr_thresh=0.05,
            pval_thresh_max=0.05,
            pval_thresh_step_size=0.00025,
            verbose=False)

        # We test the number of significant interactions for all P-value thresholds
        sig_num_o_expected = 10
        for sig_num_o in fdr_info_dict['RESULTS_TABLE']['SIG_NUM_O']:
            self.assertEqual(sig_num_o_expected, sig_num_o,
                             msg='Number of significant interactions did not increase by 10!')
            sig_num_o_expected += 10

    def test_top_1000_interactions_results(self):
        """
        Here it is tested whether the results of the FDR procedure for the test file with 1,000 input interactions have
        changed.
        """

        # The function for the FDR procedure returns a dictionary containing the numbers that are being tested here
        fdr_info_dict = self.randomize_1000.get_pval_thresh_at_chosen_fdr_thresh(
            chosen_fdr_thresh=0.05,
            pval_thresh_max=0.05,
            pval_thresh_step_size=0.00025,
            verbose=False)

        # Index of the determined P-value threshold
        result_index = fdr_info_dict['RESULT_INDEX'][0]

        # The determined P-value threshold must be 0.00750
        determined_pval_thresh = fdr_info_dict['RESULTS_TABLE']['PVAL_THRESH'][result_index]
        self.assertAlmostEqual(0.00750, determined_pval_thresh, 5,
                               msg='When this test was created, a different P-value threshold was determined!')

        # The determined -ln(P-value) threshold must be 4.89285
        determined_nnl_pval_thresh = fdr_info_dict['RESULTS_TABLE']['NNL_PVAL_THRESH'][result_index]
        self.assertAlmostEqual(4.89285, determined_nnl_pval_thresh, 5,
                               msg='When this test was created, a different -ln(P-value) threshold was determined!')

        # At the determined P-value threshold the minimum read pair number must be 9
        determined_min_rp_num = fdr_info_dict['RESULTS_TABLE']['MIN_RP_NUM'][result_index]
        self.assertEqual(9, determined_min_rp_num,
                         msg='When this test was created, a different minimum read pair number was determined!')

        # The most extreme P-value with this minimum read pair number must be 0.00391
        determined_min_rp_num_pval = fdr_info_dict['RESULTS_TABLE']['MIN_RP_NUM_PVAL'][result_index]
        self.assertAlmostEqual(0.00391, determined_min_rp_num_pval, 5,
                               msg='When this test was created, a different most extreme P-value for the minimum read '
                                   'pair number was determined!')

        # The number of interactions with 9 or more read pairs must be 686
        determined_min_rp_inter_num = fdr_info_dict['MIN_RP_NUM_INTER_NUM']
        self.assertEqual(686, determined_min_rp_inter_num,
                         msg='When this test was created, a different number of interactions with a minimum read pair '
                             'number was determined!')

        # At the determined P-value threshold the number of significant interactions must be 26
        determined_sig_num_o = fdr_info_dict['RESULTS_TABLE']['SIG_NUM_O'][result_index]
        self.assertEqual(26, determined_sig_num_o,
                         msg='When this test was created, a different number of significant interactions was '
                             'determined!')

        # At the determined P-value threshold the number of randomized significant interactions must be 1
        determined_sig_num_r = fdr_info_dict['RESULTS_TABLE']['SIG_NUM_R'][result_index]
        self.assertEqual(1, determined_sig_num_r,
                         msg='When this test was created, a different number of randomized significant interactions '
                             'was determined!')

        # At the determined P-value threshold the FDR must be 0.03846
        determined_fdr = fdr_info_dict['RESULTS_TABLE']['FDR'][result_index]
        self.assertAlmostEqual(0.03846, determined_fdr, 5, msg='When this test was created, a different FDR was '
                                                               'determined!')

        # For this input file, a warning must be issued because the number of input interactions is smaller than 16,000
        self.assertEqual(1, fdr_info_dict['WARNINGS'][0],
                         msg='When this test was created, a warning was issued because of too few input interactions!')

        # For this input file, a warning must be issued because the number interactions with a minimum number of read
        # pairs is smaller than 1,000
        self.assertEqual(1, fdr_info_dict['WARNINGS'][1],
                         msg='When this test was created, a warning was issued because of too few interactions with a '
                             'minimum number of read pairs!')

        # For this input file, a warning must be issued because the number randomized significant interactions is
        # smaller than 20
        self.assertEqual(1, fdr_info_dict['WARNINGS'][2],
                         msg='When this test was created, a warning was issued because of too few randomized '
                             'significant interactions!')

    def test_top_64000_interactions_results(self):
        """
        Here it is tested whether the results of the FDR procedure for the test file with 64,000 input interactions have
        changed.
        """

        # The function for the FDR procedure returns a dictionary containing the numbers that are being tested here
        fdr_info_dict = self.randomize_64000.get_pval_thresh_at_chosen_fdr_thresh(
            chosen_fdr_thresh=0.05,
            pval_thresh_max=0.05,
            pval_thresh_step_size=0.00025,
            verbose=False)

        # Index of the determined P-value threshold
        result_index = fdr_info_dict['RESULT_INDEX'][0]

        # The determined P-value threshold must be 0.00325
        determined_pval_thresh = fdr_info_dict['RESULTS_TABLE']['PVAL_THRESH'][result_index]
        self.assertAlmostEqual(0.00325, determined_pval_thresh, 5,
                               msg='When this test was created, a different P-value threshold was determined!')

        # The determined -ln(P-value) threshold must be 5.72910
        determined_nnl_pval_thresh = fdr_info_dict['RESULTS_TABLE']['NNL_PVAL_THRESH'][result_index]
        self.assertAlmostEqual(5.72910, determined_nnl_pval_thresh, 5,
                               msg='When this test was created, a different -ln(P-value) threshold was determined!')

        # At the determined P-value threshold the minimum read pair number must be 10
        determined_min_rp_num = fdr_info_dict['RESULTS_TABLE']['MIN_RP_NUM'][result_index]
        self.assertEqual(10, determined_min_rp_num,
                         msg='When this test was created, a different minimum read pair number was determined!')

        # The most extreme P-value with this minimum read pair number must be 0.00195
        determined_min_rp_num_pval = fdr_info_dict['RESULTS_TABLE']['MIN_RP_NUM_PVAL'][result_index]
        self.assertAlmostEqual(0.00195, determined_min_rp_num_pval, 5,
                               msg='When this test was created, a different most extreme P-value for the minimum read '
                                   'pair number was determined!')

        # The number of interactions with 9 or more read pairs must be 38,749
        determined_min_rp_inter_num = fdr_info_dict['MIN_RP_NUM_INTER_NUM']
        self.assertEqual(38749, determined_min_rp_inter_num,
                         msg='When this test was created, a different number of interactions with a minimum read pair '
                             'number was determined!')

        # At the determined P-value threshold the number of significant interactions must be 1,512
        determined_sig_num_o = fdr_info_dict['RESULTS_TABLE']['SIG_NUM_O'][result_index]
        self.assertEqual(1512, determined_sig_num_o,
                         msg='When this test was created, a different number of significant interactions was '
                             'determined!')

        # At the determined P-value threshold the number of randomized significant interactions must be 67
        determined_sig_num_r = fdr_info_dict['RESULTS_TABLE']['SIG_NUM_R'][result_index]
        self.assertEqual(67, determined_sig_num_r,
                         msg='When this test was created, a different number of randomized significant interactions '
                             'was determined!')

        # At the determined P-value threshold the FDR must be 0.04431
        determined_fdr = fdr_info_dict['RESULTS_TABLE']['FDR'][result_index]
        self.assertAlmostEqual(0.04431, determined_fdr, 5, msg='When this test was created, a different FDR was '
                                                               'determined!')

        # For this input file, no warning must be issued because the number of input interactions is smaller than 16,000
        self.assertEqual(0, fdr_info_dict['WARNINGS'][0],
                         msg='When this test was created, a warning was issued because of too few input interactions!')

        # For this input file, no warning must be issued because the number interactions with a minimum number of read
        # pairs is smaller than 1,000
        self.assertEqual(0, fdr_info_dict['WARNINGS'][1],
                         msg='When this test was created, a warning was issued because of too few interactions with a '
                             'minimum number of read pairs!')

        # For this input file, no warning must be issued because the number randomized significant interactions is
        # smaller than 20
        self.assertEqual(0, fdr_info_dict['WARNINGS'][2],
                         msg='When this test was created, a warning was issued because of too few randomized '
                             'significant interactions!')
