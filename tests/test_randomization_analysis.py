from unittest import TestCase
import os
import sys
from numpy import log
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet
from diachr.randomize_interaction_set import RandomizeInteractionSet

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))


class TestRandomizationAnalysis(TestCase):

    @classmethod
    def setUpClass(cls):

        # Prepare interaction set using the test file with 1,000 interactions
        cls.interaction_set_1000 = DiachromaticInteractionSet()
        cls.interaction_set_1000.parse_file("data/test_03/diachromatic_interaction_file_fdr_top_1000.tsv.gz")
        cls.randomize_1000 = RandomizeInteractionSet(random_seed=0)

        # Prepare interaction set using the test file with 64,000 interactions
        cls.interaction_set_64000 = DiachromaticInteractionSet()
        cls.interaction_set_64000.parse_file("data/test_03/diachromatic_interaction_file_fdr_top_64000.tsv.gz")
        cls.randomize_64000 = RandomizeInteractionSet(random_seed=0)

    def test_parallel_processing(self):
        """
        Here it is tested whether the results are independent of the degree of parallel processing.
        """

        # Parameters for randomization
        nominal_alpha = 0.0025
        iter_num = 10

        # Perform randomization analysis without 'multiprocessing' package
        random_analysis_thread_num_0_info_dict = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alphas=[nominal_alpha],
            iter_num=iter_num,
            thread_num=0)

        # Perform randomization analysis with 'multiprocessing' package but only in one process
        random_analysis_thread_num_1_info_dict = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alphas=[nominal_alpha],
            iter_num=iter_num,
            thread_num=1)

        # Perform randomization analysis with 'multiprocessing' package in two parallel processes
        random_analysis_thread_num_2_info_dict = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alphas=[nominal_alpha],
            iter_num=iter_num,
            thread_num=2)

        # Perform randomization analysis with 'multiprocessing' package in three parallel processes
        random_analysis_thread_num_3_info_dict = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alphas=[nominal_alpha],
            iter_num=iter_num,
            thread_num=3)

        # Compare results
        for iter_idx in range(0, iter_num):
            sig_num_r_0 = random_analysis_thread_num_0_info_dict['RESULTS']['SIG_NUM_R_LISTS'][nominal_alpha][iter_idx]
            sig_num_r_1 = random_analysis_thread_num_1_info_dict['RESULTS']['SIG_NUM_R_LISTS'][nominal_alpha][iter_idx]
            sig_num_r_2 = random_analysis_thread_num_2_info_dict['RESULTS']['SIG_NUM_R_LISTS'][nominal_alpha][iter_idx]
            sig_num_r_3 = random_analysis_thread_num_3_info_dict['RESULTS']['SIG_NUM_R_LISTS'][nominal_alpha][iter_idx]
            self.assertEqual(sig_num_r_0, sig_num_r_1,
                             msg='Different numbers of randomized significant interactions depending on '
                                 'parallel processing!')
            self.assertEqual(sig_num_r_0, sig_num_r_2,
                             msg='Different numbers of randomized significant interactions depending on '
                                 'parallel processing!')
            self.assertEqual(sig_num_r_0, sig_num_r_3,
                             msg='Different numbers of randomized significant interactions depending on '
                                 'parallel processing!')

    def test_for_changes_in_results(self):
        """
        Here it is tested whether the results for a certain input have remain unchanged.
        """

        random_analysis_info_dict = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alphas=[0.0025],
            iter_num=100,
            thread_num=0)

        # Find index of the passed nominal alpha
        nominal_alpha_idx = random_analysis_info_dict['INPUT_PARAMETERS']['NOMINAL_ALPHAS'].index(0.0025)

        # The following results were obtained with an earlier version
        expected_sig_num_r_mean = 51.01
        expected_sig_num_r_sd = 6.902891
        expected_z_score = 198.90073

        # The following results are obtained with the current version
        observed_sig_num_r_mean = random_analysis_info_dict['RESULTS']['SIG_NUM_R_MEAN'][nominal_alpha_idx]
        observed_sig_num_r_sd = random_analysis_info_dict['RESULTS']['SIG_NUM_R_SD'][nominal_alpha_idx]
        observed_z_score = random_analysis_info_dict['RESULTS']['Z_SCORE'][nominal_alpha_idx]

        # Compare previous and current results
        self.assertAlmostEqual(expected_sig_num_r_mean, observed_sig_num_r_mean, places=5,
                               msg='A different mean number of randomized significant interactions was obtained for an '
                                   'earlier version!')
        self.assertAlmostEqual(expected_sig_num_r_sd, observed_sig_num_r_sd, places=5,
                               msg='A different standard deviation for the number of randomized significant interactions was '
                                   'obtained for an earlier version!')
        self.assertAlmostEqual(expected_z_score, observed_z_score, places=5,
                               msg='A different Z-score was obtained for an earlier version!')

    def test_results_for_same_and_different_random_seeds(self):
        """
        Here it is tested whether the same results are obtained with the same random seed.and whether different results
        are obtained with different random seeds.
        """

        # Perform analysis for RandomizeInteractionSet object of this class (random seed: 0)
        random_analysis_info_dict = self.randomize_1000.perform_randomization_analysis(
            interaction_set=self.interaction_set_1000,
            nominal_alphas=[0.005],
            iter_num=100,
            thread_num=0)

        # Create a new RandomizeInteractionSet object with the same random seed and perform the same analysis
        randomize_1000_same_seed = RandomizeInteractionSet(random_seed=0)
        random_analysis_info_dict_same_seed = randomize_1000_same_seed.perform_randomization_analysis(
            interaction_set=self.interaction_set_1000,
            nominal_alphas=[0.005],
            iter_num=100,
            thread_num=0)

        # Create a new RandomizeInteractionSet object with the different random seed and perform the same analysis
        randomize_1000_different_seed = RandomizeInteractionSet(random_seed=42)
        random_analysis_info_dict_different_seed = randomize_1000_different_seed.perform_randomization_analysis(
            interaction_set=self.interaction_set_1000,
            nominal_alphas=[0.005],
            iter_num=100,
            thread_num=0)

        # Find indices of the passed nominal alphas
        nominal_alpha_idx = random_analysis_info_dict['INPUT_PARAMETERS']['NOMINAL_ALPHAS'].index(0.005)
        nominal_alpha_idx_same_seed = random_analysis_info_dict_same_seed['INPUT_PARAMETERS']['NOMINAL_ALPHAS'].index(0.005)
        nominal_alpha_idx_different_seed = random_analysis_info_dict_different_seed['INPUT_PARAMETERS']['NOMINAL_ALPHAS'].index(0.005)

        # Compare mean numbers of significant randomized interactions from the three objects
        sig_num_r_mean = random_analysis_info_dict['RESULTS']['SIG_NUM_R_MEAN'][nominal_alpha_idx]
        sig_num_r_mean_same_seed = random_analysis_info_dict_same_seed['RESULTS']['SIG_NUM_R_MEAN'][nominal_alpha_idx_same_seed]
        sig_num_r_mean_different_seed = \
            random_analysis_info_dict_different_seed['RESULTS']['SIG_NUM_R_MEAN'][nominal_alpha_idx_different_seed]
        self.assertEqual(sig_num_r_mean, sig_num_r_mean_same_seed, msg='Get different results with the same random '
                                                                       'seed!')
        self.assertNotEqual(sig_num_r_mean, sig_num_r_mean_different_seed, msg='Get the same results with different '
                                                                               'random seeds!')

    def test_independence_of_randomization_and_parallelization(self):
        """
        Here it is tested whether the same results are obtained for a given random seed, regardless of how many
        processes are used.
        """

        # Perform analysis without multiprocessing (random seed: 0)
        random_analysis_info_dict_0 = self.randomize_1000.perform_randomization_analysis(
            interaction_set=self.interaction_set_1000,
            nominal_alphas=[0.005],
            iter_num=100,
            thread_num=0)

        # Perform analysis with multiprocessing but only one process (random seed: 0)
        random_analysis_info_dict_1 = self.randomize_1000.perform_randomization_analysis(
            interaction_set=self.interaction_set_1000,
            nominal_alphas=[0.005],
            iter_num=100,
            thread_num=1)

        # Perform analysis with multiprocessing but only one process (random seed: 0)
        random_analysis_info_dict_2 = self.randomize_1000.perform_randomization_analysis(
            interaction_set=self.interaction_set_1000,
            nominal_alphas=[0.005],
            iter_num=100,
            thread_num=2)

        # Find indices of the passed nominal alphas
        nominal_alpha_idx_0 = random_analysis_info_dict_0['INPUT_PARAMETERS']['NOMINAL_ALPHAS'].index(0.005)
        nominal_alpha_idx_1 = random_analysis_info_dict_1['INPUT_PARAMETERS']['NOMINAL_ALPHAS'].index(0.005)
        nominal_alpha_idx_2 = random_analysis_info_dict_2['INPUT_PARAMETERS']['NOMINAL_ALPHAS'].index(0.005)

        # Compare mean numbers of significant randomized interactions from the three randomizations
        sig_num_r_mean_0 = random_analysis_info_dict_0['RESULTS']['SIG_NUM_R_MEAN'][nominal_alpha_idx_0]
        sig_num_r_mean_1 = random_analysis_info_dict_1['RESULTS']['SIG_NUM_R_MEAN'][nominal_alpha_idx_1]
        sig_num_r_mean_2 = random_analysis_info_dict_2['RESULTS']['SIG_NUM_R_MEAN'][nominal_alpha_idx_2]
        self.assertEqual(sig_num_r_mean_0, sig_num_r_mean_1, msg='Get different results with the same random depending '
                                                                 'multiprocessing!')
        self.assertEqual(sig_num_r_mean_0, sig_num_r_mean_2, msg='Get different results with the same random depending '
                                                                 'multiprocessing!')

    def test_determine_significant_pvals_at_nominal_alphas_nnl(self):
        """
        Here the function that determines the number of significant P-values for a list of nominal alphas is tested on
        a small example.
        """

        # Create list of nominal alphas
        nominal_alphas = [0.01, 0.02, 0.03, 0.04, 0.05]
        nominal_alphas_nnl = -log(nominal_alphas)

        # Create list of P-values
        p_values = [0.01, 0.011, 0.02, 0.02, 0.03, 0.03, 0.04, 0.041, 0.05, 0.05]
        p_values_nnl = -log(p_values)

        # Test the function implemented inn class RandomizeInteractionSet
        randomize = RandomizeInteractionSet()
        sig_num_list_nnl = randomize._determine_significant_pvals_at_nominal_alphas_nnl(
            nnl_nominal_alphas=nominal_alphas_nnl,
            nnl_p_values=p_values_nnl)

        # Compare the returned list with the expected list
        expected_list = [1, 4, 6, 7, 10]
        for idx in range(0, len(sig_num_list_nnl)):
            self.assertEqual(expected_list[idx], sig_num_list_nnl[idx])

    def test_results_depending_on_how_a_nominal_alpha_was_passed(self):
        """
        Here it is tested whether the results are independent of the lists in which a certain nominal alpha is passed.

        Note: The largest nominal alpha determines the size of 'RP_INTER_DICT' and thus the randomization. In order to
        keep the randomization constant, the largest nominal alpha must be the same in all lists.
        """

        # Create lists of nominal alphas
        nominal_alphas_0 = [0.01, 0.02, 0.05]
        nominal_alphas_1 = [0.01, 0.05]
        nominal_alphas_2 = [0.02, 0.05]
        nominal_alphas_3 = [0.05]

        # Pass all nominal alphas together
        random_analysis_info_dict_0 = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alphas=nominal_alphas_0,
            iter_num=10,
            thread_num=0)

        # Pass only the first and last nominal alpha of the original list
        random_analysis_info_dict_1 = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alphas=nominal_alphas_1,
            iter_num=10,
            thread_num=0)

        # Pass only the second and last nominal alpha of the original list
        random_analysis_info_dict_2 = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alphas=nominal_alphas_2,
            iter_num=10,
            thread_num=0)

        # Pass only the third nominal alpha of the original list
        random_analysis_info_dict_3 = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alphas=nominal_alphas_3,
            iter_num=10,
            thread_num=0)

        # Compare Z-scores for first nominal alpha of the original list
        nominal_alpha_idx = random_analysis_info_dict_0['INPUT_PARAMETERS']['NOMINAL_ALPHAS'].index(nominal_alphas_0[0])
        z_score_0_0 = random_analysis_info_dict_0['RESULTS']['Z_SCORE'][nominal_alpha_idx]
        nominal_alpha_idx = random_analysis_info_dict_1['INPUT_PARAMETERS']['NOMINAL_ALPHAS'].index(nominal_alphas_0[0])
        z_score_1_0 = random_analysis_info_dict_1['RESULTS']['Z_SCORE'][nominal_alpha_idx]
        self.assertEqual(z_score_0_0, z_score_1_0, msg='Results differ depending on how a nominal alpha was passed!')

        # Compare Z-scores for second nominal alpha of the original list
        nominal_alpha_idx = random_analysis_info_dict_0['INPUT_PARAMETERS']['NOMINAL_ALPHAS'].index(nominal_alphas_0[1])
        z_score_0_1 = random_analysis_info_dict_0['RESULTS']['Z_SCORE'][nominal_alpha_idx]
        nominal_alpha_idx = random_analysis_info_dict_2['INPUT_PARAMETERS']['NOMINAL_ALPHAS'].index(nominal_alphas_0[1])
        z_score_2_1 = random_analysis_info_dict_2['RESULTS']['Z_SCORE'][nominal_alpha_idx]
        self.assertEqual(z_score_0_1, z_score_2_1, msg='Results differ depending on how a nominal alpha was passed!')

        # Compare Z-scores for third nominal alpha of the original list
        nominal_alpha_idx = random_analysis_info_dict_0['INPUT_PARAMETERS']['NOMINAL_ALPHAS'].index(nominal_alphas_0[2])
        z_score_0_2 = random_analysis_info_dict_0['RESULTS']['Z_SCORE'][nominal_alpha_idx]
        nominal_alpha_idx = random_analysis_info_dict_1['INPUT_PARAMETERS']['NOMINAL_ALPHAS'].index(nominal_alphas_0[2])
        z_score_1_2 = random_analysis_info_dict_1['RESULTS']['Z_SCORE'][nominal_alpha_idx]
        nominal_alpha_idx = random_analysis_info_dict_2['INPUT_PARAMETERS']['NOMINAL_ALPHAS'].index(nominal_alphas_0[2])
        z_score_2_2 = random_analysis_info_dict_2['RESULTS']['Z_SCORE'][nominal_alpha_idx]
        nominal_alpha_idx = random_analysis_info_dict_3['INPUT_PARAMETERS']['NOMINAL_ALPHAS'].index(nominal_alphas_0[2])
        z_score_3_2 = random_analysis_info_dict_3['RESULTS']['Z_SCORE'][nominal_alpha_idx]
        self.assertEqual(z_score_0_2, z_score_1_2, msg='Results differ depending on how a nominal alpha was passed!')
        self.assertEqual(z_score_0_2, z_score_2_2, msg='Results differ depending on how a nominal alpha was passed!')
        self.assertEqual(z_score_0_2, z_score_3_2, msg='Results differ depending on how a nominal alpha was passed!')
