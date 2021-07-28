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

    def test_remove_digest_length_outliers(self):
        """
        Test whether interactions with extreme digest pair length are removed properly.
        """

        # First five interactions in the test file:
        # chr14	43059116	43059494	N	chr14	43101212	43101810	N	100:0:1:0	378	598 <- removed
        # chr8	129042054	129044258	N	chr8	129121269	129121986	N	100:0:2:0	2204	717 <- removed in second pass
        # chr15	73467156	73468652	N	chr15	73526903	73528438	N	100:0:3:0	1496	1535
        # chr17	72411026	72411616	N	chr17	72712662	72724357	E	100:0:1:0	590	11695 <- removed
        # chr18	38724804	38726198	N	chr18	76794986	76803172	E	100:0:2:0	1394	8186 <- removed
        # chr11	114362648	114362686	N	chr11	114396073	114404234	E	100:0:3:0	38	8161 <- removed
        # chr15	56158017	56158267	N	chr15	56462978	56465983	E	100:0:4:0	250	3005 <- removed
        # chr14	34714080	34716362	E	chr14	50135355	50139051	N	100:0:1:0	2282	3696

        # Remove outliers
        self.interaction_set_4.remove_digest_length_outliers(dg_min_len=500, dg_max_len=10000, dg_min_len_q=0.25)

        # Test first three interactions
        d_inter_list = list(self.interaction_set_4.interaction_list)
        self.assertEqual('8:129042054:8:129121269', d_inter_list[0].key)
        self.assertEqual('15:73467156:15:73526903', d_inter_list[1].key)
        self.assertEqual('14:34714080:14:50135355', d_inter_list[2].key)

        # Remove one additional outlier
        self.interaction_set_4.remove_digest_length_outliers(dg_min_len=718, dg_max_len=10000, dg_min_len_q=0.25)
        d_inter_list = list(self.interaction_set_4.interaction_list)

        # Test first three interactions
        self.assertEqual('15:73467156:15:73526903', d_inter_list[0].key)
        self.assertEqual('14:34714080:14:50135355', d_inter_list[1].key)
