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
        We use test data to create a DiachromaticInteractionSet object.
        For the test file, there are 18 clearly directed and 28 clearly undirected interactions.
        In addition, there are four interactions that do not have enough read pairs to be significant at the given
        P-value threshold.
        For 16 directed interactions, there are undirected reference interactions with a corresponding number of read
        pairs, but for two interactions there is no matching reference interaction.
        One interaction is in the enrichment category 'NE' and has a total of 104 read pairs.
        The other interaction is in enrichment category 'EE' and has a total of 106 read pairs.
        """

        cls.interaction_set = DiachromaticInteractionSet()
        cls.interaction_set.parse_file("data/test_04/diachromatic_interaction_file.tsv")
        cls.rate_and_cat_report_dict = cls.interaction_set.rate_and_categorize_interactions(-log(0.01))

    def test_rate_and_categorize_report_dict(self):

        # Chosen P-value threshold
        self.assertAlmostEqual(4.60517018, self.rate_and_cat_report_dict['NLN_PVAL_THRESH'], 7)

        # Minimum number of read pairs required for significance
        self.assertEqual(8, self.rate_and_cat_report_dict['MIN_RP'])

        # Largest significant P-value
        self.assertAlmostEqual(0.0078125, self.rate_and_cat_report_dict['MIN_RP_PVAL'], 7)

        # Total number of interactions
        self.assertEqual(50, self.rate_and_cat_report_dict['N_PROCESSED'])

        # Number of discarded interactions (not enough read pairs)
        self.assertEqual(4, self.rate_and_cat_report_dict['N_DISCARDED'])

        # Number of discarded interactions
        self.assertEqual(18, self.rate_and_cat_report_dict['N_DIRECTED'])

        # Number of discarded interactions
        self.assertEqual(28, self.rate_and_cat_report_dict['N_UNDIRECTED'])

    def test_rate_and_categorize_created_file(self):
        pass


    def test_reference_selection_report_dict(self):
        """
        XXX
        """

        # The function for reference selection returns a dictionary with interaction numbers for all categories.
        select_ref_report_dict = self.interaction_set.select_reference_interactions()

        # Directed interactions
        self.assertEqual(3, select_ref_report_dict['DI_NN'])
        self.assertEqual(4, select_ref_report_dict['DI_NE'])
        self.assertEqual(5, select_ref_report_dict['DI_EN'])
        self.assertEqual(6, select_ref_report_dict['DI_EE'])

        # Undirected reference interactions
        self.assertEqual(3, select_ref_report_dict['UIR_NN'])
        self.assertEqual(3, select_ref_report_dict['UIR_NE'])
        self.assertEqual(5, select_ref_report_dict['UIR_EN'])
        self.assertEqual(5, select_ref_report_dict['UIR_EE'])

        # Missing undirected reference interactions
        self.assertEqual(0, select_ref_report_dict['M_UIR_NN'])
        self.assertEqual(1, select_ref_report_dict['M_UIR_NE'])
        self.assertEqual(0, select_ref_report_dict['M_UIR_EN'])
        self.assertEqual(1, select_ref_report_dict['M_UIR_EE'])

        # Undirected reference interactions
        self.assertEqual(3, select_ref_report_dict['UI_NN'])
        self.assertEqual(3, select_ref_report_dict['UI_NE'])
        self.assertEqual(3, select_ref_report_dict['UI_EN'])
        self.assertEqual(3, select_ref_report_dict['UI_EE'])


    def test_reference_selection_created_file(self):
        """
        This test creates an interaction file with P-values and interaction categories and tests the number of
        interactions for all interaction ('DI', 'UIR', 'UI') and enrichment categories ('NN', 'NE', 'EN', 'EE').
        It is also tested whether the interactions in the various categories have the expected number of read pairs.
        """

        # Select reference interactions and write interaction file
        self.interaction_set.select_reference_interactions()
        self.interaction_set.write_diachromatic_interaction_file(target_file_name='i_file.tsv')

        # Read interaction file
        # ---------------------

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

        # There must be 3 undirected interactions in category NE (one missing)
        self.assertEqual(3, len(rp_inter_dict['NE']['UI']))

        # There must be 3 undirected interactions in category EN
        self.assertEqual(3, len(rp_inter_dict['EN']['UI']))

        # There must be 3 undirected interactions in category EE (one missing)
        self.assertEqual(3, len(rp_inter_dict['EE']['UI']))

        # Test total read pair numbers in different enrichment categories
        self.assertTrue(101 in rp_inter_dict['NN']['UIR'])
        self.assertTrue(102 in rp_inter_dict['NN']['UIR'])
        self.assertTrue(103 in rp_inter_dict['NN']['UIR'])

        self.assertTrue(101 in rp_inter_dict['NE']['UIR'])
        self.assertTrue(102 in rp_inter_dict['NE']['UIR'])
        self.assertTrue(103 in rp_inter_dict['NE']['UIR'])

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

        # Remove interaction file
        os.remove('i_file.tsv')
