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
        parser = DiachromaticInteractionSet()
        for file in os.listdir(diachromatic_input_dir):
            ifile_path = os.path.join(diachromatic_input_dir, file)
            if file.endswith(".tsv.gz"):
                parser.parse_file(i_file=ifile_path)
        # For testing, reach into the parser object and get a dictionary of interactions
        cls.interaction_dict = parser.interaction_dict

    def test_get_four_interactions(self):
        """
        There are a total of four interactions in the data
        """
        n_interactions = len(self.interaction_dict)
        self.assertEqual(4, n_interactions)


    def test_get_three_above_threshold_interactions(self):
        """
        There are a total of four interactions in the data
        """
        n_above_threshold = 0
        threshold = self.required_reps  # threshold number of replicates
        for _, interaction in self.interaction_dict.items():
            if interaction.has_data_for_required_replicate_num(threshold):
                n_above_threshold += 1
        self.assertEqual(3, n_above_threshold)

    def test_specific_interaction(self):
        """
        Test one of the interactions, that should be
                 chr17   72411026        72411616        I       chr17   72712662        72724357        I       9:6
        """
        k = '17:72411026:17:72712662'
        self.assertTrue(k in self.interaction_dict)
        interaction = self.interaction_dict.get(k)
        self.assertEqual("chr17", interaction.chrA)
        self.assertEqual("chr17", interaction.chrB)
        self.assertEqual(9, interaction.n_simple)
        self.assertEqual(6, interaction.n_twisted)
        self.assertEqual('NN', interaction.enrichment_status_tag_pair)
