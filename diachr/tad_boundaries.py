from bisect import bisect_left
import random


class TadBoundarySet:
    # Class to store TAD boundaries with a function that returns the next boundary for a given coordinate

    # Attributes

    n_tads = 0  # Total number of TAD boundaries

    chr_tad_boundary_dict = {}  # Dictionary that has a list of boundary coordinates for each chromosome

    # Initializer

    def __init__(self, tad_boundary_bed_file: str, chr_size_file: str = None):

        self.chr_tad_boundary_dict = {}
        with open(tad_boundary_bed_file, 'rt') as fp:
            next(fp)
            for line in fp:

                # Parse line
                chr_key, boundary_1, boundary_plus_one = line.rstrip().split('\t')

                # Count TAD boundaries
                self.n_tads = self.n_tads + 2

                # Add boundary coordinates to dictionary
                if chr_key in self.chr_tad_boundary_dict:
                    self.chr_tad_boundary_dict[chr_key].append(int(boundary_1))
                else:
                    self.chr_tad_boundary_dict[chr_key] = [int(boundary_1)]

        # Sort lists with coordinates for each chromosome
        for chr_key in self.chr_tad_boundary_dict:
            self.chr_tad_boundary_dict[chr_key] = sorted(self.chr_tad_boundary_dict[chr_key])

        # Create dictionary for chromosome sizes required for randomization
        if chr_size_file is not None:
            self.chr_size_dict = {}
            with open(chr_size_file, 'rt') as fp:
                for line in fp:
                    chr_key, size = line.rstrip().rsplit('\t')
                    self.chr_size_dict[chr_key] = int(size)

    def get_nearest_tad_boundary(self, chr_key, coord):

        # Get the index of x (if in the list) or the index of next larger coordinate
        boundary_to_the_right_idx = bisect_left(self.chr_tad_boundary_dict[chr_key], int(coord))

        if len(self.chr_tad_boundary_dict[chr_key]) <= boundary_to_the_right_idx:
            # There is no boundary to the right
            return self.chr_tad_boundary_dict[chr_key][boundary_to_the_right_idx - 1]

        if boundary_to_the_right_idx == 0:
            # There is no boundary to the left
            return self.chr_tad_boundary_dict[chr_key][boundary_to_the_right_idx]

        if self.chr_tad_boundary_dict[chr_key][boundary_to_the_right_idx] == coord:
            # coord is a TAD boundary and the distance is zero
            return coord
        else:
            # coord is between two TAD boundaries
            boundary_to_the_left_idx = boundary_to_the_right_idx - 1

            if (self.chr_tad_boundary_dict[chr_key][boundary_to_the_right_idx] - coord) < (
                    coord - self.chr_tad_boundary_dict[chr_key][boundary_to_the_left_idx]):
                # The boundary to the right is closer
                return self.chr_tad_boundary_dict[chr_key][boundary_to_the_right_idx]
            else:
                # The boundary to the left is closer
                return self.chr_tad_boundary_dict[chr_key][boundary_to_the_left_idx]

    def get_distance_to_nearest_tad_boundary(self, chr_key, coord):

        # Check whether there are TAD boundaries on this chromosome
        if chr_key not in self.chr_tad_boundary_dict:
            return -1

        # Determine distance to next TAD boundary
        nearest_tad_boundary_coord = self.get_nearest_tad_boundary(chr_key, coord)
        if coord <= nearest_tad_boundary_coord:
            return nearest_tad_boundary_coord - coord
        else:
            return coord - nearest_tad_boundary_coord

    def get_distance_to_nearest_tad_boundary(self, chr_key, coord):

        # Check whether there are TAD boundaries on this chromosome
        if chr_key not in self.chr_tad_boundary_dict:
            return -1

        # Determine distance to next TAD boundary
        nearest_tad_boundary_coord = self.get_nearest_tad_boundary(chr_key, coord)
        if coord <= nearest_tad_boundary_coord:
            return nearest_tad_boundary_coord - coord
        else:
            return coord - nearest_tad_boundary_coord

    def get_number_of_boundaries_spanned_by_region(self, chr_key, sta, end):

        # Check whether there are TAD boundaries on this chromosome
        if chr_key not in self.chr_tad_boundary_dict:
            return -1

        num_of_spanned_boundaries = 0

        # Get the index of sta (if in the list) or the index of next larger coordinate
        boundary_to_the_right_idx = bisect_left(self.chr_tad_boundary_dict[chr_key], int(sta))

        if len(self.chr_tad_boundary_dict[chr_key]) <= boundary_to_the_right_idx:
            # There is no boundary to the right
            return num_of_spanned_boundaries

        if self.chr_tad_boundary_dict[chr_key][boundary_to_the_right_idx] == sta:
            # sta is a TAD boundary
            num_of_spanned_boundaries += 1
            boundary_to_the_right_idx += 1

        while boundary_to_the_right_idx < len(self.chr_tad_boundary_dict[chr_key]) and \
                self.chr_tad_boundary_dict[chr_key][boundary_to_the_right_idx] < end:
            num_of_spanned_boundaries += 1
            boundary_to_the_right_idx += 1

        return num_of_spanned_boundaries

    def randomize_tad_boundaries(self, random_seed: int):

        random.seed(random_seed)

        for chr_key in self.chr_tad_boundary_dict.keys():
            # Get number of TAD boundaries on this chromosome
            tb_num = len(self.chr_tad_boundary_dict[chr_key])

            # Draw the a number of random positions that corresponds to the number of TAD boundaries on this chromosome
            random_tb_list = random.sample(range(0, self.chr_size_dict[chr_key]), tb_num)

            # Sort random list and replace original list
            self.chr_tad_boundary_dict[chr_key] = sorted(random_tb_list)

    def randomize_tad_boundaries_2(self, random_seed: int, random_range: int = 100000):

        random.seed(random_seed)

        for chr_key in self.chr_tad_boundary_dict.keys():
            random_tb_list = []
            for tb in self.chr_tad_boundary_dict[chr_key]:
                if self.chr_size_dict[chr_key] < tb + random_range:
                    tb_random = random.randint(tb - random_range, self.chr_size_dict[chr_key])
                elif tb - random_range < 0:
                    tb_random = random.randint(0, tb + random_range)
                else:
                    tb_random = random.randint(tb - random_range, tb + random_range)
                random_tb_list.append(tb_random)

            # Sort random list and replace original list
            self.chr_tad_boundary_dict[chr_key] = sorted(random_tb_list)
