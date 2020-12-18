from unittest import TestCase
import os
import sys
from numpy import log

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from diachr.diachromatic_interaction_parser import DiachromaticInteractionParser

class TestReferenceInteractionSelection(TestCase):

    @classmethod
    def setUpClass(cls):

        # Use test data to create a DiachromaticInteractionParser object
        cls.interaction_set = DiachromaticInteractionParser()
        cls.interaction_set.parse_file("data/test_04/diachromatic_interaction_file.tsv")
        cls.interaction_set.rate_and_categorize_interactions(-log(0.01))

    def test_num_of_missing_reference_interactions(self):
        """
        For the test file, there are two directed interactions for which no matching undirected interaction
        can be selected. One interaction is in the enrichment category 'NE' and has a total of 104 read pairs.
        The other interaction is in enrichment category 'EE' and has a total of 106 read pairs.
        """

        # The function for reference selection returns a dictionary with information about missing reference interactions
        report, missing_ref_info  = self.interaction_set.select_reference_interactions()

        # Check that there are two missing reference interactions
        self.assertEqual(2,sum(missing_ref_info.values()))

        # Check that these interactions and have correct enrichment categories and read pair numbers
        self.assertTrue("NE:104" in missing_ref_info)
        self.assertTrue("EE:106" in missing_ref_info)

        # Write interaction file
        self.interaction_set.write_diachromatic_interaction_file(target_file_name='i_file.tsv')

        # Remove interaction file
        os.remove('i_file.tsv')




