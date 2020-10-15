from unittest import TestCase
import os
import sys
PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))


from diachr2 import EnhancedInteractionParser, EnhancedInteraction


class TestEnhancedInteractionParser(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.enhanced_interaction_file = os.path.join(os.path.dirname(__file__), 'data', "small_enhanced_interaction_file.tsv")
        cls.parser = EnhancedInteractionParser(cls.enhanced_interaction_file)
        # Also just get individual lines for testing
        cls.lines = []
        with open(cls.enhanced_interaction_file) as f:
            for line in f:
                cls.lines.append(line)

    def test_zeroth_line(self):
        """
        Test we get the right information from the zeroth line
        chr6:89716777-89716882;chr6:89916736-89920828	199854	DIII	;	2:18	II	8.51	-1/-1	;
        """
        parser = EnhancedInteractionParser(self.enhanced_interaction_file)
        ei = parser.parse_line(self.lines[0])
        self.assertEqual('chr6', ei.chr_a)
        self.assertEqual(89716777, ei.sta_a)
        self.assertEqual(89716882, ei.end_a)
        self.assertEqual('chr6', ei.chr_b)
        self.assertEqual(89916736, ei.sta_b)
        self.assertEqual(89920828, ei.end_b)
        self.assertEqual("",ei.syms_a)
        self.assertEqual("",ei.syms_b)
        self.assertEqual("",ei.tsss_a)
        self.assertEqual("",ei.tsss_b)
        self.assertEqual("II", ei.enrichment_pair_tag)
        self.assertEqual("-1/-1", ei.strand_pair_tag)
        self.assertEqual("DIII", ei.interaction_category)
        self.assertAlmostEqual(8.51, ei.neg_log_p_value)
        self.assertEqual(20, ei.rp_total)  # 2:18 means 20 read pairs in total
        self.assertEqual(199854, ei.i_dist)
        #  interaction_category, neg_log_p_value, rp_total, i_dist 