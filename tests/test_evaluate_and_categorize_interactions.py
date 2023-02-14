from unittest import TestCase
import os
import sys
import gzip

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from diachr.diachromatic_interaction_set import DiachromaticInteractionSet


class TestRateAndCategorizeInteractions(TestCase):

    @classmethod
    def setUpClass(cls):
        """
        We use test data to create a DiachromaticInteractionSet object:

               data/test_04/diachromatic_interaction_file.tsv

        For the test file, there are 18 clearly unbalanced and 28 clearly balanced interactions.
        In addition, there are four interactions that do not have enough read pairs to be significant at the given
        P-value threshold. For 16 unbalanced interactions, there are balanced reference interactions with a
        corresponding number of read pairs, but for two interactions there is no matching reference interaction.
        These two interactions are moved to a separate category UX.

        First, the test data is used to create an object of class 'DiachromaticInteractionSet'.
        There are three processing steps that are tested here: (1) classification of unbalanced interactions,
        (2) selection of reference interaction sets, and (3) writing new Diachromatic11 interaction files with
        additional columns for interaction category and P-value.
        """

        cls.interaction_set = DiachromaticInteractionSet(rpc_rule='st')
        test_dir = os.path.dirname(__file__)
        interaction_file = os.path.join(test_dir, "data/test_04/diachromatic_interaction_file.tsv")
        cls.interaction_set.parse_file(interaction_file)
        cls.rate_and_cat_report_dict = cls.interaction_set.evaluate_and_categorize_interactions(0.01)
        cls.select_ref_report_dict = cls.interaction_set.select_reference_interactions()

    def test_rate_and_categorize_report_dict(self):
        """
        In the first processing step, P-values are calculated and, according to a given threshold, divided into
        unbalanced and balanced. Interactions that do not have enough read pairs to be significant are discarded.
        The function for evaluation and categorization returns a dictionary with intermediate results that are tested
        here.
        """

        # Chosen P-value threshold
        self.assertAlmostEqual(0.01, self.rate_and_cat_report_dict['PVAL_THRESH'][0], 7)

        # Minimum number of read pairs required for significance
        self.assertEqual(8, self.rate_and_cat_report_dict['MIN_RP'][0])

        # Largest significant P-value
        self.assertAlmostEqual(0.0078125, self.rate_and_cat_report_dict['MIN_RP_PVAL'][0], 7)

        # Total number of interactions
        self.assertEqual(50, self.rate_and_cat_report_dict['N_PROCESSED'][0])

        # Number of discarded interactions (not enough read pairs)
        self.assertEqual(4, self.rate_and_cat_report_dict['N_DISCARDED'][0])

        # Number of unbalanced interactions
        self.assertEqual(18, self.rate_and_cat_report_dict['N_UNBALANCED'][0])

        # Number of balanced interactions
        self.assertEqual(28, self.rate_and_cat_report_dict['N_BALANCED'][0])

    def test_reference_selection_report_dict(self):
        """
        In a next step, balanced reference interactions (BR) are selected from the balanced interactions (U).
        The function for reference selection returns a dictionary with intermediate results that are tested here.
        """

        # Unbalanced interactions
        self.assertEqual(3, self.select_ref_report_dict['NN']['UR'][0])
        self.assertEqual(3, self.select_ref_report_dict['NE']['UR'][0])
        self.assertEqual(5, self.select_ref_report_dict['EN']['UR'][0])
        self.assertEqual(5, self.select_ref_report_dict['EE']['UR'][0])

        # Unbalanced interactions without reference (104 and 106 read pairs)
        self.assertEqual(1, self.select_ref_report_dict['NE']['UX'][0])
        self.assertEqual(1, self.select_ref_report_dict['EE']['UX'][0])

        # Balanced reference interactions
        self.assertEqual(3, self.select_ref_report_dict['NN']['BR'][0])
        self.assertEqual(3, self.select_ref_report_dict['NE']['BR'][0])
        self.assertEqual(5, self.select_ref_report_dict['EN']['BR'][0])
        self.assertEqual(5, self.select_ref_report_dict['EE']['BR'][0])

        # Balanced interactions not selected as reference
        self.assertEqual(3, self.select_ref_report_dict['NN']['BX'][0])
        self.assertEqual(3, self.select_ref_report_dict['NE']['BX'][0])
        self.assertEqual(3, self.select_ref_report_dict['EN']['BX'][0])
        self.assertEqual(3, self.select_ref_report_dict['EE']['BX'][0])

    def test_reference_selection_created_file(self):
        """
        Finally, the interactions are written to Diachromatic interaction file, which has two additional columns on the
        left for P-value and interaction category ('UX', 'UR', 'BR', 'BX').
        In this test, the file is read in again and the interaction numbers in the various categories are compared
        with the known numbers.
        """

        # Create interaction file with P-values and interaction categories
        self.interaction_set.write_diachromatic_interaction_file(target_file='i_file.tsv')

        # Read created interaction file
        # -----------------------------

        # Nested dictionaries that store the numbers of interactions (value) for different read pair numbers (key)
        rp_inter_dict = {'NN': {'UX': {}, 'UR': {}, 'BR': {}, 'BX': {}},
                         'NE': {'UX': {}, 'UR': {}, 'BR': {}, 'BX': {}},
                         'EN': {'UX': {}, 'UR': {}, 'BR': {}, 'BX': {}},
                         'EE': {'UX': {}, 'UR': {}, 'BR': {}, 'BX': {}}}

        with gzip.open('i_file.tsv', 'rt') as fp:
            for line in fp:
                F = line.rstrip().split('\t')
                interaction_category = F[10]
                enrichment_pair_tag = F[3] + F[7]
                rp_total = int(F[8].split(':')[0]) + int(F[8].split(':')[1]) + int(F[8].split(':')[2]) + int(
                    F[8].split(':')[3])

                if rp_total not in rp_inter_dict[enrichment_pair_tag]:
                    rp_inter_dict[enrichment_pair_tag][interaction_category][rp_total] = 1
                else:
                    rp_inter_dict[enrichment_pair_tag][interaction_category][rp_total] += 1

        # Check values collected from file
        # --------------------------------

        # There must be 3 unbalanced interactions in category NN
        self.assertEqual(3, len(rp_inter_dict['NN']['UR']))

        # There must be 3 unbalanced interactions in category NE
        self.assertEqual(3, len(rp_inter_dict['NE']['UR']))

        # There must be 5 unbalanced interactions in category EN
        self.assertEqual(5, len(rp_inter_dict['EN']['UR']))

        # There must be 5 unbalanced interactions in category EE
        self.assertEqual(5, len(rp_inter_dict['EE']['UR']))

        # For one unbalanced interactions in category NE with 104 read pairs there is no balanced reference interaction
        self.assertEqual(1, len(rp_inter_dict['NE']['UX']))

        # For one unbalanced interactions in category EE with 104 read pairs there is no balanced reference interaction
        self.assertEqual(1, len(rp_inter_dict['EE']['UX']))

        # There must be 3 balanced reference interactions in category NN
        self.assertEqual(3, len(rp_inter_dict['NN']['BR']))

        # There must be 3 balanced reference interactions in category NE (one missing)
        self.assertEqual(3, len(rp_inter_dict['NE']['BR']))

        # There must be 5 balanced reference interactions in category EN
        self.assertEqual(5, len(rp_inter_dict['EN']['BR']))

        # There must be 5 balanced reference interactions in category EE (one missing)
        self.assertEqual(5, len(rp_inter_dict['EE']['BR']))

        # There must be 3 balanced interactions in category NN
        self.assertEqual(3, len(rp_inter_dict['NN']['BX']))

        # There must be 3 balanced interactions in category NE
        self.assertEqual(3, len(rp_inter_dict['NE']['BX']))

        # There must be 3 balanced interactions in category EN
        self.assertEqual(3, len(rp_inter_dict['EN']['BX']))

        # There must be 3 balanced interactions in category EE
        self.assertEqual(3, len(rp_inter_dict['EE']['BX']))

        # Test total read pair numbers in different enrichment categories
        self.assertTrue(101 in rp_inter_dict['NN']['BR'])
        self.assertTrue(102 in rp_inter_dict['NN']['BR'])
        self.assertTrue(103 in rp_inter_dict['NN']['BR'])

        self.assertTrue(101 in rp_inter_dict['NE']['BR'])
        self.assertTrue(102 in rp_inter_dict['NE']['BR'])
        self.assertTrue(103 in rp_inter_dict['NE']['BR'])
        self.assertFalse(104 in rp_inter_dict['NE']['BR'])  # (missing)

        self.assertTrue(101 in rp_inter_dict['EN']['BR'])
        self.assertTrue(102 in rp_inter_dict['EN']['BR'])
        self.assertTrue(103 in rp_inter_dict['EN']['BR'])
        self.assertTrue(104 in rp_inter_dict['EN']['BR'])
        self.assertTrue(105 in rp_inter_dict['EN']['BR'])

        self.assertTrue(101 in rp_inter_dict['EE']['BR'])
        self.assertTrue(102 in rp_inter_dict['EE']['BR'])
        self.assertTrue(103 in rp_inter_dict['EE']['BR'])
        self.assertTrue(104 in rp_inter_dict['EE']['BR'])
        self.assertTrue(105 in rp_inter_dict['EE']['BR'])
        self.assertFalse(106 in rp_inter_dict['EE']['BR'])  # (missing)

        # Remove created interaction file
        os.remove('i_file.tsv')
