from unittest import TestCase
import os
import sys

#PACKAGE_PARENT = '..'
#SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.#path.expanduser(__file__))))
#sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from diachr import TadBoundarySet


class TestTadToundary(TestCase):
    @classmethod
    def setUpClass(cls):
        test_dir = os.path.dirname(__file__)
        print("TESTDIR", test_dir)
        tad_file = os.path.join(test_dir, "data/tad_boundary/tads_test.bed")
        cls.tadBS = TadBoundarySet(tad_region_bed_file=tad_file)
    
    def test_file_not_null(self):
        self.assertIsNotNone(self.tadBS)

    def test_get_dist1(self):
        """
        Our test file has these TADs
        chr1	8346815	8969941
        chr1	9112441	9494941
        find the distance to some nearby points
        """
        chrom = 'chr1'
        pos = 8969941 + 1
        dist = self.tadBS.get_distance_to_nearest_tad_boundary(chr_key=chrom, coord=pos)
        self.assertEqual(1, dist)
        pos = 8969941 + 2
        dist = self.tadBS.get_distance_to_nearest_tad_boundary(chr_key=chrom, coord=pos)
        self.assertEqual(2, dist)
        pos = 9112441 - 1
        dist = self.tadBS.get_distance_to_nearest_tad_boundary(chr_key=chrom, coord=pos)
        self.assertEqual(1, dist)
        ## What happens if we are within a TAD?
        pos = 8969941 - 11
        dist = self.tadBS.get_distance_to_nearest_tad_boundary(chr_key=chrom, coord=pos)
        self.assertEqual(11, dist)


    def test_chromosome_with_no_tads(self):
        """
        -1 is the signal that a chromosome has no TADs (in this case,
        the chromosome will not be in the map.) In analysis code, 
        such interactions are simply omitted.
        """
        invalid_chromosome = 'chr42'
        pos = 42
        dist = self.tadBS.get_distance_to_nearest_tad_boundary(invalid_chromosome, pos)
        self.assertEqual(-1, dist)

   



