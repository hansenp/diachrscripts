from unittest import TestCase
import os
import sys

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from diachr.diachromatic_interaction_set import DiachromaticInteractionSet

class TestCombineInteractions(TestCase):

    @classmethod
    def setUpClass(cls):

        # Use test toy example data to create a DiachromaticInteractionParser object
        cls.interaction_set = DiachromaticInteractionSet()
        current_dir =  os.path.dirname(__file__)
        cls.interaction_set.parse_file(os.path.join(current_dir,"data/test_01/diachromatic_interaction_file_r1.tsv.gz")) # one interaction
        cls.interaction_set.parse_file(os.path.join(current_dir,"data/test_01/diachromatic_interaction_file_r2.tsv.gz")) # two interactions
        cls.interaction_set.parse_file(os.path.join(current_dir,"data/test_01/diachromatic_interaction_file_r3.tsv.gz")) # three interactions
        cls.interaction_set.parse_file(os.path.join(current_dir,"data/test_01/diachromatic_interaction_file_r4.tsv.gz")) # four interactions

        # Use real data sample files (top 64,000 of Megakaryocyte data) to create a DiachromaticInteractionParser object
        cls.interaction_set_real_top_64000 = DiachromaticInteractionSet()
        cls.interaction_set_real_top_64000.parse_file(os.path.join(current_dir,"data/test_01/real_data/diachromatic_interaction_file_top_64000_r1.tsv.gz")) # 64,000 interaction
        cls.interaction_set_real_top_64000.parse_file(os.path.join(current_dir,"data/test_01/real_data/diachromatic_interaction_file_top_64000_r2.tsv.gz")) # 64,000 interactions
        cls.interaction_set_real_top_64000.parse_file(os.path.join(current_dir,"data/test_01/real_data/diachromatic_interaction_file_top_64000_r3.tsv.gz")) # 64,000 interactions
        cls.interaction_set_real_top_64000.parse_file(os.path.join(current_dir,"data/test_01/real_data/diachromatic_interaction_file_top_64000_r4.tsv.gz")) # 64,000 interactions

        # Read real data files in reverse order
        cls.interaction_set_real_top_64000_rev = DiachromaticInteractionSet()
        cls.interaction_set_real_top_64000_rev.parse_file(os.path.join(current_dir,"data/test_01/real_data/diachromatic_interaction_file_top_64000_r4.tsv.gz")) # 64,000 interaction
        cls.interaction_set_real_top_64000_rev.parse_file(os.path.join(current_dir,"data/test_01/real_data/diachromatic_interaction_file_top_64000_r3.tsv.gz")) # 64,000 interactions
        cls.interaction_set_real_top_64000_rev.parse_file(os.path.join(current_dir,"data/test_01/real_data/diachromatic_interaction_file_top_64000_r2.tsv.gz")) # 64,000 interactions
        cls.interaction_set_real_top_64000_rev.parse_file(os.path.join(current_dir,"data/test_01/real_data/diachromatic_interaction_file_top_64000_r1.tsv.gz")) # 64,000 interactions

    def test_num_of_combined_inter(self):
        """
        We have parsed 1 + 2 + 3 + 4 = 10 interactions from four interactions files.
        We combine interaction that have the same coordinates.
        For the test data, this must result in 4 combined interactions.
        """

        i_num = len(self.interaction_set.interaction_list)
        self.assertEqual(4, i_num) # equal up to 5 significant digits


    def test_num_of_combined_inter_with_required_num_of_rep(self):
        """
        We filter combined interactions based on the number of replicates in which they occur.
        This function tests the number of combined interactions for different numbers of required replicates.
        """

        required_replicates = 1 # There must be 4 interactions that occur in at least one replicate
        n_inter_with_required_replicate_num = 0
        for i in self.interaction_set.interaction_list:
            if i.has_data_for_required_replicate_num(required_replicates):
                n_inter_with_required_replicate_num += 1
        self.assertEqual(4, n_inter_with_required_replicate_num)

        required_replicates = 2 # There must be 3 interactions that occur in at least two replicates
        n_inter_with_required_replicate_num = 0
        for i in self.interaction_set.interaction_list:
            if i.has_data_for_required_replicate_num(required_replicates):
                n_inter_with_required_replicate_num += 1
        self.assertEqual(3, n_inter_with_required_replicate_num)

        required_replicates = 3 # There must be 2 interactions that occur in at least three replicates
        n_inter_with_required_replicate_num = 0
        for i in self.interaction_set.interaction_list:
            if i.has_data_for_required_replicate_num(required_replicates):
                n_inter_with_required_replicate_num += 1
        self.assertEqual(2, n_inter_with_required_replicate_num)

        required_replicates = 4  # There is only one interactions that occurs in all four replicates
        n_inter_with_required_replicate_num = 0
        for i in self.interaction_set.interaction_list:
            if i.has_data_for_required_replicate_num(required_replicates):
                n_inter_with_required_replicate_num += 1
        self.assertEqual(1, n_inter_with_required_replicate_num)


    def test_summation_of_simple_and_twisted_rp_counts(self):
        """
        For combined interactions, we add up the read pair counts from different replicates separately
        for simple and twisted.
        """

        for i in self.interaction_set.interaction_list:

            # The interaction on chromosome 'chr1' occurs in all four replicates
            # Simple:Twisted counts: 2:1 + 2:1 + 2:1 + 2:1 -> 8:4
            if i.chrA == "chr1":
                self.assertEqual(8, i.n_simple)
                self.assertEqual(4, i.n_twisted)

            # The interaction on chromosome 'chr17' occurs in three replicates
            # Simple:Twisted counts: 3:2 + 3:2 + 3:2 -> 9:6
            if i.chrA == "chr17":
                self.assertEqual(9, i.n_simple)
                self.assertEqual(6, i.n_twisted)

            # The interaction on chromosome 'chr7' occurs in two replicates
            # Simple:Twisted counts: 4:3 + 4:3 -> 8:6
            if i.chrA == "chr7":
                self.assertEqual(8, i.n_simple)
                self.assertEqual(6, i.n_twisted)

            # The interaction on chromosome 'chr11' occurs in one replicate
            # Simple:Twisted counts: 5:4 -> 5:4
            if i.chrA == "chr11":
                self.assertEqual(5, i.n_simple)
                self.assertEqual(4, i.n_twisted)

    def test_for_changes_in_results(self):
        """
        Here it is tested whether the results for a certain input have remain unchanged.
        """

        # Test by how many interactions the interaction set grew with each read file
        self.assertEqual(64000, self.interaction_set_real_top_64000._read_file_info_dict['I_SET_SIZE'][0])
        self.assertEqual(113804, self.interaction_set_real_top_64000._read_file_info_dict['I_SET_SIZE'][1])
        self.assertEqual(157600, self.interaction_set_real_top_64000._read_file_info_dict['I_SET_SIZE'][2])
        self.assertEqual(207364, self.interaction_set_real_top_64000._read_file_info_dict['I_SET_SIZE'][3])

        # Test number of interactions that occur in at least one replicate
        self.interaction_set_real_top_64000.write_diachromatic_interaction_file(required_replicates=1, target_file='foo')
        os.remove('foo')
        self.assertEqual(0, self.interaction_set_real_top_64000._write_file_info_dict['N_INCOMPLETE_DATA'][0])
        self.assertEqual(207364, self.interaction_set_real_top_64000._write_file_info_dict['N_COMPLETE_DATA'][0])

        # Test number of interactions that occur in at least two replicates
        self.interaction_set_real_top_64000.write_diachromatic_interaction_file(required_replicates=2, target_file='foo')
        os.remove('foo')
        self.assertEqual(174777, self.interaction_set_real_top_64000._write_file_info_dict['N_INCOMPLETE_DATA'][0])
        self.assertEqual(32587, self.interaction_set_real_top_64000._write_file_info_dict['N_COMPLETE_DATA'][0])

        # Test number of interactions that occur in at least three replicates
        self.interaction_set_real_top_64000.write_diachromatic_interaction_file(required_replicates=3, target_file='foo')
        os.remove('foo')
        self.assertEqual(195738, self.interaction_set_real_top_64000._write_file_info_dict['N_INCOMPLETE_DATA'][0])
        self.assertEqual(11626, self.interaction_set_real_top_64000._write_file_info_dict['N_COMPLETE_DATA'][0])

        # Test number of interactions that occur in all four replicates
        self.interaction_set_real_top_64000.write_diachromatic_interaction_file(required_replicates=4, target_file='foo')
        os.remove('foo')
        self.assertEqual(202941, self.interaction_set_real_top_64000._write_file_info_dict['N_INCOMPLETE_DATA'][0])
        self.assertEqual(4423, self.interaction_set_real_top_64000._write_file_info_dict['N_COMPLETE_DATA'][0])

    def test_influence_of_reading_order_on_final_result(self):

        # Regardless of the reading order, the interaction set must end up with the same number of interactions
        self.assertEqual(self.interaction_set_real_top_64000._read_file_info_dict['I_SET_SIZE'][3],
                         self.interaction_set_real_top_64000_rev._read_file_info_dict['I_SET_SIZE'][3])

        # Test number of interactions that occur in at least one, two, three and four replicates
        for required_replicates in [1, 2, 3, 4]:
            self.interaction_set_real_top_64000.write_diachromatic_interaction_file(required_replicates=required_replicates, target_file='foo')
            os.remove('foo')
            self.interaction_set_real_top_64000_rev.write_diachromatic_interaction_file(required_replicates=required_replicates, target_file='foo')
            os.remove('foo')
            self.assertEqual(self.interaction_set_real_top_64000._write_file_info_dict['N_INCOMPLETE_DATA'][0],
                             self.interaction_set_real_top_64000_rev._write_file_info_dict['N_INCOMPLETE_DATA'][0])
            self.assertEqual(self.interaction_set_real_top_64000._write_file_info_dict['N_COMPLETE_DATA'][0],
                             self.interaction_set_real_top_64000_rev._write_file_info_dict['N_COMPLETE_DATA'][0])
