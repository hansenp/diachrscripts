from unittest import TestCase
import os
import sys
import gzip
from numpy import log

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from diachr.diachromatic_interaction_set import DiachromaticInteractionSet

class TestReferenceInteractionSelection(TestCase):

    @classmethod
    def setUpClass(cls):

        # Use test data to create a DiachromaticInteractionParser object
        cls.interaction_set = DiachromaticInteractionSet()
        cls.interaction_set.parse_file("data/test_04/diachromatic_interaction_file.tsv")
        cls.interaction_set.rate_and_categorize_interactions(-log(0.01))

    def test_missing_ref_info(self):
        """
        For the test file, there are two directed interactions for which no matching undirected interaction
        can be selected. One interaction is in the enrichment category 'NE' and has a total of 104 read pairs.
        The other interaction is in enrichment category 'EE' and has a total of 106 read pairs.
        """

        # The function for reference selection returns a dictionary with information about missing reference interactions
        select_ref_report, select_ref_summary_stat_dict, missing_ref_info = self.interaction_set.select_reference_interactions()

        # Check that there are two missing reference interactions
        self.assertEqual(2, sum(missing_ref_info.values()))

        # Check that these interactions and have correct enrichment categories and read pair numbers
        self.assertTrue("NE:104" in missing_ref_info)
        self.assertTrue("EE:106" in missing_ref_info)

    def test_interaction_file_after_reference_selection(self):
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
