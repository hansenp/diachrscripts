from unittest import TestCase
import os
import sys
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
        Here it is tested whether the results are of the degree of parallel processing.
        """

        # Parameters for randomization
        nominal_alpha = 0.0025
        iter_num = 10

        # Perform randomization analysis without 'multiprocessing' package
        random_analysis_thread_num_0_info_dict = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alpha=nominal_alpha,
            iter_num=iter_num,
            thread_num=0)

        # Perform randomization analysis with 'multiprocessing' package but only in one process
        random_analysis_thread_num_1_info_dict = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alpha=nominal_alpha,
            iter_num=iter_num,
            thread_num=1)

        # Perform randomization analysis with 'multiprocessing' package in two parallel processes
        random_analysis_thread_num_2_info_dict = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alpha=nominal_alpha,
            iter_num=iter_num,
            thread_num=2)

        # Perform randomization analysis with 'multiprocessing' package in three parallel processes
        random_analysis_thread_num_3_info_dict = self.randomize_64000.perform_randomization_analysis(
            interaction_set=self.interaction_set_64000,
            nominal_alpha=nominal_alpha,
            iter_num=iter_num,
            thread_num=3)

        # Compare results
        for iter_idx in range(0, iter_num):
            sig_num_r_0 = random_analysis_thread_num_0_info_dict['RESULTS']['SIG_NUM_R_LIST'][iter_idx]
            sig_num_r_1 = random_analysis_thread_num_1_info_dict['RESULTS']['SIG_NUM_R_LIST'][iter_idx]
            sig_num_r_2 = random_analysis_thread_num_2_info_dict['RESULTS']['SIG_NUM_R_LIST'][iter_idx]
            sig_num_r_3 = random_analysis_thread_num_3_info_dict['RESULTS']['SIG_NUM_R_LIST'][iter_idx]
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
            nominal_alpha=0.0025,
            iter_num=100,
            thread_num=0)

        # The following results were obtained with an earlier version
        expected_sig_num_r_mean = 51.01
        expected_sig_num_r_sd = 6.902891
        expected_z_score = 198.90073

        # The following results are obtained with the current version
        observed_sig_num_r_mean = random_analysis_info_dict['RESULTS']['SUMMARY']['SIG_NUM_R_MEAN'][0]
        observed_sig_num_r_sd = random_analysis_info_dict['RESULTS']['SUMMARY']['SIG_NUM_R_STD'][0]
        observed_z_score = random_analysis_info_dict['RESULTS']['SUMMARY']['Z_SCORE'][0]

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
            nominal_alpha=0.005,
            iter_num=100,
            thread_num=0)

        # Create a new RandomizeInteractionSet object with the same random seed and perform the same analysis
        randomize_1000_same_seed = RandomizeInteractionSet(random_seed=0)
        random_analysis_info_dict_same_seed = randomize_1000_same_seed.perform_randomization_analysis(
            interaction_set=self.interaction_set_1000,
            nominal_alpha=0.005,
            iter_num=100,
            thread_num=0)

        # Create a new RandomizeInteractionSet object with the different random seed and perform the same analysis
        randomize_1000_different_seed = RandomizeInteractionSet(random_seed=42)
        random_analysis_info_dict_different_seed = randomize_1000_different_seed.perform_randomization_analysis(
            interaction_set=self.interaction_set_1000,
            nominal_alpha=0.005,
            iter_num=100,
            thread_num=0)

        # Compare mean numbers of significant randomized interactions from the three objects
        sig_num_r_mean = random_analysis_info_dict['RESULTS']['SUMMARY']['SIG_NUM_R_MEAN'][0]
        sig_num_r_mean_same_seed = random_analysis_info_dict_same_seed['RESULTS']['SUMMARY']['SIG_NUM_R_MEAN'][0]
        sig_num_r_mean_different_seed = \
            random_analysis_info_dict_different_seed['RESULTS']['SUMMARY']['SIG_NUM_R_MEAN'][0]
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
            nominal_alpha=0.005,
            iter_num=100,
            thread_num=0)

        # Perform analysis with multiprocessing but only one process (random seed: 0)
        random_analysis_info_dict_1 = self.randomize_1000.perform_randomization_analysis(
            interaction_set=self.interaction_set_1000,
            nominal_alpha=0.005,
            iter_num=100,
            thread_num=1)

        # Perform analysis with multiprocessing but only one process (random seed: 0)
        random_analysis_info_dict_2 = self.randomize_1000.perform_randomization_analysis(
            interaction_set=self.interaction_set_1000,
            nominal_alpha=0.005,
            iter_num=100,
            thread_num=2)

        # Compare mean numbers of significant randomized interactions from the three randomizations
        sig_num_r_mean_0 = random_analysis_info_dict_0['RESULTS']['SUMMARY']['SIG_NUM_R_MEAN'][0]
        sig_num_r_mean_1 = random_analysis_info_dict_1['RESULTS']['SUMMARY']['SIG_NUM_R_MEAN'][0]
        sig_num_r_mean_2 = random_analysis_info_dict_2['RESULTS']['SUMMARY']['SIG_NUM_R_MEAN'][0]
        self.assertEqual(sig_num_r_mean_0, sig_num_r_mean_1, msg='Get different results with the same random depending '
                                                                 'multiprocessing!')
        self.assertEqual(sig_num_r_mean_0, sig_num_r_mean_2, msg='Get different results with the same random depending '
                                                                 'multiprocessing!')
