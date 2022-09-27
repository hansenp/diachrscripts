from unittest import TestCase
import os
import sys
from diachr.TIMViz import CHCTadViz

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))


class TestCHCTadViz(TestCase):

    @classmethod
    def setUpClass(cls):

        test_dir = os.path.dirname(__file__)

        # Prepare test data
        interaction_file = os.path.join(test_dir, "data/CHCTadViz_i_test_file.tsv")
        tad_file = os.path.join(test_dir, "data/CHCTadViz_t_test_file.bed")
        cls.CHCTadViz = CHCTadViz(i_file=interaction_file, t_file=tad_file)

    def test_conversion(self):
        """
        Test whether the conversion from genomic coordinates to plot coordinates and vice versa works correctly.
        """

        # Set parameters required for conversion
        figure_size = 10
        begin = 12400000
        end = 16300000
        self.CHCTadViz._figure_size = figure_size
        self.CHCTadViz._min = begin
        self.CHCTadViz._max = end
        self.CHCTadViz._span = end - begin
        self.CHCTadViz._factor = self.CHCTadViz._figure_size / self.CHCTadViz._span

        # Somewhere at the beginning of the region
        pos_1 = begin + 20000
        coord = self.CHCTadViz.pos_to_coordinate(pos_1)
        pos_2 = self.CHCTadViz.coordinate_to_pos(coord)
        self.assertEqual(pos_1, pos_2)

        # Somewhere at the end of the region
        pos_1 = end - 20000
        coord = self.CHCTadViz.pos_to_coordinate(pos_1)
        pos_2 = self.CHCTadViz.coordinate_to_pos(coord)
        self.assertEqual(pos_1, pos_2)

        # In the middle
        pos_1 = begin + ((end - begin)/2)
        coord = self.CHCTadViz.pos_to_coordinate(pos_1)
        pos_2 = self.CHCTadViz.coordinate_to_pos(coord)
        self.assertEqual(pos_1, pos_2)
