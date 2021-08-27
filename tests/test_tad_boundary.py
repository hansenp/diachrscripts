from unittest import TestCase
import os

# PACKAGE_PARENT = '..'
# SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.#path.expanduser(__file__))))
# sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from diachr import TadBoundarySet


class TestTadBoundary(TestCase):
    @classmethod
    def setUpClass(cls):
        test_dir = os.path.dirname(__file__)
        tad_boundary_file = os.path.join(test_dir, "data/tad_boundary/tad_boundary_test.bed")
        cls.tadBS = TadBoundarySet(tad_boundary_bed_file=tad_boundary_file)

    def test_file_not_null(self):
        self.assertIsNotNone(self.tadBS)

    def test_get_dist1(self):
        """
        Our test file has these TAD boundary positions
        chr1	8346815	8346816
        chr1	9112441	9112442
        find the distance to some nearby points
        """
        chrom = 'chr1'
        pos = 8346815 + 1
        dist = self.tadBS.get_distance_to_nearest_tad_boundary(chr_key=chrom, coord=pos)
        self.assertEqual(1, dist)
        pos = 8346815 + 2
        dist = self.tadBS.get_distance_to_nearest_tad_boundary(chr_key=chrom, coord=pos)
        self.assertEqual(2, dist)
        pos = 8346815 - 1
        dist = self.tadBS.get_distance_to_nearest_tad_boundary(chr_key=chrom, coord=pos)
        self.assertEqual(1, dist)
        pos = 8346815 - 11
        dist = self.tadBS.get_distance_to_nearest_tad_boundary(chr_key=chrom, coord=pos)
        self.assertEqual(11, dist)

    def test_chromosome_with_no_tad_boundaries(self):
        """
        -1 is the signal that a chromosome has no TADs (in this case,
        the chromosome will not be in the map.) In analysis code, 
        such interactions are simply omitted.
        """
        invalid_chromosome = 'chr42'
        pos = 42
        dist = self.tadBS.get_distance_to_nearest_tad_boundary(invalid_chromosome, pos)
        self.assertEqual(-1, dist)

    def test_get_number_of_boundaries_spanned_by_region(self):
        """
        Test a region on chromosome 'chr1' that spans three tad boundaries. The first position of the
        region is one position behind a TAD boundary and the last position is one position before a
        TAD boundary.
        chr1	8346815	8346816
        chr1	9112441	9112442
        chr1	9724942	9724943
        chr1	10187442	10187443
        chr1	11021193	11021194
        """

        tad_boundary_pos_1 = 8346815 + 1  # One position behind a TAD boundary
        # 3 tad boundaries in between
        tad_boundary_pos_2 = 11021193  # One position before a TAD boundary

        # Expect three TAD boundaries in between
        num_of_spanned_boundaries = self.tadBS.get_number_of_boundaries_spanned_by_region(
            'chr1', tad_boundary_pos_1, tad_boundary_pos_2)
        self.assertEqual(3, num_of_spanned_boundaries)

        # Moving start position of region one position to the left
        tad_boundary_pos_1 -= 1

        # Expect four TAD boundaries in between
        num_of_spanned_boundaries = self.tadBS.get_number_of_boundaries_spanned_by_region(
            'chr1', tad_boundary_pos_1, tad_boundary_pos_2)
        self.assertEqual(4, num_of_spanned_boundaries)

        # Moving end position of region one position to the right
        tad_boundary_pos_2 += 1

        # Expect five TAD boundaries in between
        num_of_spanned_boundaries = self.tadBS.get_number_of_boundaries_spanned_by_region(
            'chr1', tad_boundary_pos_1, tad_boundary_pos_2)
        self.assertEqual(5, num_of_spanned_boundaries)

        # Move start position between second and third boundary
        tad_boundary_pos_1 = 9112441 + 10

        # Expect three TAD boundaries in between
        num_of_spanned_boundaries = self.tadBS.get_number_of_boundaries_spanned_by_region(
            'chr1', tad_boundary_pos_1, tad_boundary_pos_2)
        self.assertEqual(3, num_of_spanned_boundaries)

    def test_chromosome_with_no_tad_boundaries_spanned_boundaries(self):
        """
        -1 is the signal that a chromosome has no TAD boundary (in this case,
        the chromosome will not be in the map.) In analysis code,
        such interactions are simply omitted.
        """
        invalid_chromosome = 'chr42'
        pos = 42
        num_of_spanned_boundaries = self.tadBS.get_number_of_boundaries_spanned_by_region(invalid_chromosome, pos, pos + 100)
        self.assertEqual(-1, num_of_spanned_boundaries)



