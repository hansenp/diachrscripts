from unittest import TestCase
import os
import sys
import gzip
from numpy import log

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

        For the test file, there are 18 clearly directed and 28 clearly undirected interactions.
        In addition, there are four interactions that do not have enough read pairs to be significant at the given
        P-value threshold. For 16 directed interactions, there are undirected reference interactions with a
        corresponding number of read pairs, but for two interactions there is no matching reference interaction.

        First, the test data is used to create an object of class 'DiachromaticInteractionSet'.
        Within this project there are three processing steps that are tested here.
        """

        cls.interaction_set = DiachromaticInteractionSet()
        cls.interaction_set.parse_file("data/test_04/diachromatic_interaction_file.tsv")
        cls.rate_and_cat_report_dict = cls.interaction_set.evaluate_and_categorize_interactions(-log(0.01))
        cls.select_ref_report_dict = cls.interaction_set.select_reference_interactions()

    def test_rate_and_categorize_report_dict(self):
        """
        In the first processing step, P-values are calculated and, according to a given threshold, divided into
        directed and undirected. Interactions that do not have enough read pairs to be significant are discarded.
        The function for evaluation and categorization returns a dictionary with intermediate results that are tested
        here.
        """

        # Chosen P-value threshold
        self.assertAlmostEqual(4.60517018, self.rate_and_cat_report_dict['NLN_PVAL_THRESH'][0], 7)

        # Minimum number of read pairs required for significance
        self.assertEqual(8, self.rate_and_cat_report_dict['MIN_RP'][0])

        # Largest significant P-value
        self.assertAlmostEqual(0.0078125, self.rate_and_cat_report_dict['MIN_RP_PVAL'][0], 7)

        # Total number of interactions
        self.assertEqual(50, self.rate_and_cat_report_dict['N_PROCESSED'][0])

        # Number of discarded interactions (not enough read pairs)
        self.assertEqual(4, self.rate_and_cat_report_dict['N_DISCARDED'][0])

        # Number of directed interactions
        self.assertEqual(18, self.rate_and_cat_report_dict['N_DIRECTED'][0])

        # Number of undirected interactions
        self.assertEqual(28, self.rate_and_cat_report_dict['N_UNDIRECTED'][0])


    def test_reference_selection_report_dict(self):
        """
        In a next step, undirected reference interactions (UIR) are selected from the undirected interactions (UI).
        The function for reference selection returns a dictionary with intermediate results that are tested here.
        """

        # Directed interactions
        self.assertEqual(3, self.select_ref_report_dict['DI_NN'][0])
        self.assertEqual(4, self.select_ref_report_dict['DI_NE'][0])
        self.assertEqual(5, self.select_ref_report_dict['DI_EN'][0])
        self.assertEqual(6, self.select_ref_report_dict['DI_EE'][0])

        # Undirected reference interactions
        self.assertEqual(3, self.select_ref_report_dict['UIR_NN'][0])
        self.assertEqual(3, self.select_ref_report_dict['UIR_NE'][0])
        self.assertEqual(5, self.select_ref_report_dict['UIR_EN'][0])
        self.assertEqual(5, self.select_ref_report_dict['UIR_EE'][0])

        # Missing undirected reference interactions
        self.assertEqual(0, self.select_ref_report_dict['M_UIR_NN'][0])
        self.assertEqual(1, self.select_ref_report_dict['M_UIR_NE'][0])
        self.assertEqual(0, self.select_ref_report_dict['M_UIR_EN'][0])
        self.assertEqual(1, self.select_ref_report_dict['M_UIR_EE'][0])

        # Undirected reference interactions
        self.assertEqual(3, self.select_ref_report_dict['UI_NN'][0])
        self.assertEqual(3, self.select_ref_report_dict['UI_NE'][0])
        self.assertEqual(3, self.select_ref_report_dict['UI_EN'][0])
        self.assertEqual(3, self.select_ref_report_dict['UI_EE'][0])

    def test_reference_selection_created_file(self):
        """
        Finally, the interactions are written to Diachromatic interaction file, which has two additional columns on the
        left for P-value and interaction category ('DI', 'UIR', 'UI').
        Inn this test, the file is read in again and the interaction numbers in the various categories are compared
        with the known numbers.
        """

        # Create interaction file with P-values and interaction categories
        self.interaction_set.write_diachromatic_interaction_file(target_file='i_file.tsv')

        # Read created interaction file
        # -----------------------------

        # Nested dictionaries that store the numbers of interactions (value) for different read pair numbers (key)
        rp_inter_dict = {'NN': {'DI':{},'UIR':{},'UI':{}},
                         'NE': {'DI':{},'UIR':{},'UI':{}},
                         'EN': {'DI':{},'UIR':{},'UI':{}},
                         'EE': {'DI':{},'UIR':{},'UI':{}}}

        with gzip.open('i_file.tsv', 'rt') as fp:
            for line in fp:
                F = line.rstrip().split('\t')
                interaction_category = F[10]
                enrichment_pair_tag = F[3] + F[7]
                rp_total = int(F[8].split(':')[0]) + int(F[8].split(':')[1])

                if rp_total not in rp_inter_dict[enrichment_pair_tag]:
                    rp_inter_dict[enrichment_pair_tag][interaction_category][rp_total] = 1
                else:
                    rp_inter_dict[enrichment_pair_tag][interaction_category][rp_total] += 1

        # Check values collected from file
        # --------------------------------

        # There must be 3 directed interactions in category NN
        self.assertEqual(3, len(rp_inter_dict['NN']['DI']))

        # There must be 4 directed interactions in category NE
        self.assertEqual(4, len(rp_inter_dict['NE']['DI']))

        # There must be 5 directed interactions in category EN
        self.assertEqual(5, len(rp_inter_dict['EN']['DI']))

        # There must be 6 directed interactions in category EE
        self.assertEqual(6, len(rp_inter_dict['EE']['DI']))

        # There must be 3 undirected reference interactions in category NN
        self.assertEqual(3, len(rp_inter_dict['NN']['UIR']))

        # There must be 3 undirected reference interactions in category NE (one missing)
        self.assertEqual(3, len(rp_inter_dict['NE']['UIR']))

        # There must be 5 undirected reference interactions in category EN
        self.assertEqual(5, len(rp_inter_dict['EN']['UIR']))

        # There must be 5 undirected reference interactions in category EE (one missing)
        self.assertEqual(5, len(rp_inter_dict['EE']['UIR']))

        # There must be 3 undirected interactions in category NN
        self.assertEqual(3, len(rp_inter_dict['NN']['UI']))

        # There must be 3 undirected interactions in category NE
        self.assertEqual(3, len(rp_inter_dict['NE']['UI']))

        # There must be 3 undirected interactions in category EN
        self.assertEqual(3, len(rp_inter_dict['EN']['UI']))

        # There must be 3 undirected interactions in category EE
        self.assertEqual(3, len(rp_inter_dict['EE']['UI']))

        # Test total read pair numbers in different enrichment categories
        self.assertTrue(101 in rp_inter_dict['NN']['UIR'])
        self.assertTrue(102 in rp_inter_dict['NN']['UIR'])
        self.assertTrue(103 in rp_inter_dict['NN']['UIR'])

        self.assertTrue(101 in rp_inter_dict['NE']['UIR'])
        self.assertTrue(102 in rp_inter_dict['NE']['UIR'])
        self.assertTrue(103 in rp_inter_dict['NE']['UIR'])
        self.assertFalse(104 in rp_inter_dict['NE']['UIR']) # (missing)

        self.assertTrue(101 in rp_inter_dict['EN']['UIR'])
        self.assertTrue(102 in rp_inter_dict['EN']['UIR'])
        self.assertTrue(103 in rp_inter_dict['EN']['UIR'])
        self.assertTrue(104 in rp_inter_dict['EN']['UIR'])
        self.assertTrue(105 in rp_inter_dict['EN']['UIR'])

        self.assertTrue(101 in rp_inter_dict['EE']['UIR'])
        self.assertTrue(102 in rp_inter_dict['EE']['UIR'])
        self.assertTrue(103 in rp_inter_dict['EE']['UIR'])
        self.assertTrue(104 in rp_inter_dict['EE']['UIR'])
        self.assertTrue(105 in rp_inter_dict['EE']['UIR'])
        self.assertFalse(106 in rp_inter_dict['EE']['UIR']) # (missing)

        # Remove created interaction file
        os.remove('i_file.tsv')
