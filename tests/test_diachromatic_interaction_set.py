from unittest import TestCase
import os
import sys

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from diachr.diachromatic_interaction_set import DiachromaticInteractionSet


class TestDiachromaticInteractionSet(TestCase):

    @classmethod
    def setUpClass(cls):
        """
        We use test data to create a DiachromaticInteractionSet object:

               data/test_04/diachromatic_interaction_file.tsv
        """

        cls.interaction_set_1 = DiachromaticInteractionSet()
        cls.interaction_set_2 = DiachromaticInteractionSet()
        test_dir = os.path.dirname(__file__)
        interaction_file = os.path.join(test_dir, "data/test_04/diachromatic_interaction_file.tsv")
        cls.interaction_set_1.parse_file(interaction_file)
        cls.interaction_set_2.parse_file(interaction_file)

    def test_shuffle_inter_dict(self):
        """
        Test whether executing the function 'DiachromaticInteractionSet.shuffle_inter_dict' has the effect that
        interactions are reproducibly processed in a different order.
        """

        # Before shuffling, the first digests of the first three interactions are on chromosomes 'chr14', 'chr8' and
        # 'chr15' (as in the file).
        i = 0
        for d_inter in self.interaction_set_1.interaction_list:
            if i == 0:
                self.assertEqual('chr14', d_inter.chrA)
            if i == 1:
                self.assertEqual('chr8', d_inter.chrA)
            if i == 2:
                self.assertEqual('chr15', d_inter.chrA)
            i += 1

        # After shuffling with random seed 1, the first digests of the first three interactions are on chromosomes
        # 'chr15', 'chr10' and 'chr7'.
        self.interaction_set_1.shuffle_inter_dict(random_seed=1)
        i = 0
        for d_inter in self.interaction_set_1.interaction_list:
            if i == 0:
                self.assertEqual('chr15', d_inter.chrA)
            if i == 1:
                self.assertEqual('chr10', d_inter.chrA)
            if i == 2:
                self.assertEqual('chr7', d_inter.chrA)
            i += 1

        # After shuffling another DiachromaticInteractionSet with random seed 1, the first digests of the first three
        # interactions are on chromosomes 'chr15', 'chr10' and 'chr7'.
        self.interaction_set_2.shuffle_inter_dict(random_seed=1)
        i = 0
        for d_inter in self.interaction_set_2.interaction_list:
            if i == 0:
                self.assertEqual('chr15', d_inter.chrA)
            if i == 1:
                self.assertEqual('chr10', d_inter.chrA)
            if i == 2:
                self.assertEqual('chr7', d_inter.chrA)
            i += 1
