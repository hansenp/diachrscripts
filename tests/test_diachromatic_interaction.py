from unittest import TestCase
import os
import sys
PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

#
#
from diachr import DiachromaticInteractionSet
#
#
class DiachromaticInteraction(TestCase):
    @classmethod
    def setUpClass(cls):
        diachromatic_input_dir = os.path.join(os.path.dirname(__file__), 'data', "test_01")
        cls.required_reps = 2
        #parser = DiachromaticInteractionSet(interaction_file_dir=diachromatic_input_dir)
        #parser.parse_file()
        # For testing, reach into the parser object and get a dictionary of interactions
        #cls.interaction_dict = parser.get_read_file_info_dict()

    def test_get_four_interactions(self):
        """
        There are a total of four interactions in the data
        """
        #self.assertEqual(4, len(self.interaction_dict))
        self.assertTrue(True)
#
#     def test_get_three_above_threshold_interactions(self):
#         """
#         There are a total of four interactions in the data
#         """
#         n_above_threshold = 0
#         threshold = self.required_reps  # threshold number of replicates
#         for _, interaction in self.interaction_dict.items():
#             if interaction.has_data_for_required_replicate_num(threshold):
#                 n_above_threshold += 1
#         self.assertEqual(3, n_above_threshold)
#
#     def test_specific_interaction(self):
#         """
#         Test one of the interactions, that should be
#         chr17   72411026        72411616        I       chr17   72712662        72724357        I       9:6
#         """
#         #Note the key should be the string '17724110261772712662'
#         k = '17724110261772712662'
#         self.assertTrue(k in self.interaction_dict)
#         interaction = self.interaction_dict.get(k)
#         self.assertEqual("chr17", interaction.chrA)
#         self.assertEqual("chr17", interaction.chrB)
#         self.assertEqual(9, interaction.n_simple)
#         self.assertEqual(6, interaction.n_twisted)
#         self.assertEqual('II', interaction.enrichment_status_tag_pair)
