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
        cls.interaction_set_3 = DiachromaticInteractionSet()
        cls.interaction_set_st = DiachromaticInteractionSet(rpc_rule='st')
        cls.interaction_set_ht = DiachromaticInteractionSet(rpc_rule='ht')
        cls.interaction_set_4 = DiachromaticInteractionSet()
        test_dir = os.path.dirname(__file__)
        interaction_file = os.path.join(test_dir, "data/test_04/diachromatic_interaction_file.tsv")
        cls.interaction_set_1.parse_file(interaction_file)
        cls.interaction_set_2.parse_file(interaction_file)
        cls.interaction_set_3.parse_file(interaction_file)
        cls.interaction_set_st.parse_file(interaction_file)
        cls.interaction_set_ht.parse_file(interaction_file)
        cls.interaction_set_4.parse_file(interaction_file)

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

    def test_read_pair_count_rules_two_rule(self):
        """
        How the four counts are used to caalculate the P-value depends on the RPC rule, which is an attribute of
        'DiachromaticInteractionSet'. Here it is tested whether the aassignment of counts is independent of the RPC
        rule and whether the binomial P-values are calculated correctly.
        """

        # Calculate P-values for the two interaction sets
        self.interaction_set_st.evaluate_and_categorize_interactions(0.01)
        self.interaction_set_ht.evaluate_and_categorize_interactions(0.01)

        # Get lists of interactions to test
        d_inter_list_st = list(self.interaction_set_st.interaction_list)
        d_inter_list_ht = list(self.interaction_set_ht.interaction_list)

        # The four counts must be independent of the RPC rule
        for i in range(0, len(d_inter_list_st)):
            self.assertEqual(d_inter_list_st[i]._simple_1,d_inter_list_ht[i]._simple_1)
            self.assertEqual(d_inter_list_st[i]._simple_2, d_inter_list_ht[i]._simple_2)
            self.assertEqual(d_inter_list_st[i]._twisted_1, d_inter_list_ht[i]._twisted_1)
            self.assertEqual(d_inter_list_st[i]._twisted_2, d_inter_list_ht[i]._twisted_2)

        # But the P-values have to be calculated differently
        self.assertAlmostEqual(28.09, d_inter_list_st[0]._log10_pval, places=2)
        self.assertAlmostEqual(30.40, d_inter_list_ht[0]._log10_pval, places=2)

        self.assertAlmostEqual(26.68, d_inter_list_st[1]._log10_pval, places=2)
        self.assertAlmostEqual(30.71, d_inter_list_ht[1]._log10_pval, places=2)

        self.assertAlmostEqual(25.44, d_inter_list_st[2]._log10_pval, places=2)
        self.assertAlmostEqual(31.01, d_inter_list_ht[2]._log10_pval, places=2)
