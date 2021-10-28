from bisect import bisect_left
import random
import copy


class TadBoundarySet:
    # Class to store TAD boundaries with a function that returns the next boundary for a given coordinate

    # Attributes

    n_tad_boundaries = 0  # Total number of TAD boundaries

    chr_tad_boundary_dict = {}  # Dictionary that has a list of boundary coordinates for each chromosome

    chr_size_dict = {}

    # Initializer
    def __init__(self, tad_boundary_bed_file: str, chr_size_file: str = None):

        self.chr_tad_boundary_dict = {}
        with open(tad_boundary_bed_file, 'rt') as fp:
            next(fp)
            for line in fp:

                # Parse line
                chr_key, boundary_1, boundary_2 = line.rstrip().split('\t')
                boundary_1 = int(boundary_1)
                boundary_2 = int(boundary_2)

                # Add boundary coordinates to dictionary
                if chr_key in self.chr_tad_boundary_dict:
                    self.chr_tad_boundary_dict[chr_key].append(int(boundary_1))
                    self.n_tad_boundaries += 1
                    if 1 < boundary_2 - boundary_1:
                        self.chr_tad_boundary_dict[chr_key].append(int(boundary_2))
                        self.n_tad_boundaries += 1
                else:
                    self.chr_tad_boundary_dict[chr_key] = [int(boundary_1)]
                    self.n_tad_boundaries += 1
                    if 1 < boundary_2 - boundary_1:
                        self.chr_tad_boundary_dict[chr_key].append(int(boundary_2))
                        self.n_tad_boundaries += 1

        # Sort lists with coordinates for each chromosome
        for chr_key in self.chr_tad_boundary_dict:
            self.chr_tad_boundary_dict[chr_key] = sorted(self.chr_tad_boundary_dict[chr_key])

        # Create dictionary for chromosome sizes required for randomization
        if chr_size_file is not None:
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

    def get_randomized_tad_boundary_set(self, random_seed: int, random_range: int = 0):
        """
        This function implements two ways of randomizing TAD boundaries.
        If the 'random_range' is zero, then, for each chromosome,
        a corresponding number of TAD boundaries is randomly drawn from all positions of the chromosome.
        If the 'random_range' is greater than zero, then for each TAD boundary,
        a position is randomly drawn from the surrounding sequence.
        If the 'random_range' is greater than zero,
        then an error is reported and the TadBoundarySet remains unchanged.
        """

        random.seed(random_seed)

        random_tbs = copy.deepcopy(self)

        if random_range == 0:

            # Draw position from all positions of a chromosome
            for chr_key in random_tbs.chr_tad_boundary_dict.keys():
                # Get number of TAD boundaries on this chromosome
                tb_num = len(random_tbs.chr_tad_boundary_dict[chr_key])

                # Draw the a number of random positions
                # that corresponds to the number of TAD boundaries on this chromosome
                random_tb_list = random.sample(range(0, random_tbs.chr_size_dict[chr_key]), tb_num)

                # Sort random list and replace original list
                random_tbs.chr_tad_boundary_dict[chr_key] = sorted(random_tb_list)

            return random_tbs

        elif 0 < random_range:

            # Draw position from the surrounding region of each TAD boundary
            for chr_key in random_tbs.chr_tad_boundary_dict.keys():
                random_tb_list = []
                for tb in random_tbs.chr_tad_boundary_dict[chr_key]:
                    if random_tbs.chr_size_dict[chr_key] < tb + random_range:
                        # Not enough space for surrounding region to the right
                        tb_random = random.randint(tb - random_range, random_tbs.chr_size_dict[chr_key])
                    elif tb - random_range < 0:
                        # Not enough space for surrounding region to the left
                        tb_random = random.randint(0, tb + random_range)
                    else:
                        tb_random = random.randint(tb - random_range, tb + random_range)
                    random_tb_list.append(tb_random)

                # Sort random list and replace original list
                random_tbs.chr_tad_boundary_dict[chr_key] = sorted(random_tb_list)

            return random_tbs

        else:
            print("[ERROR] \'random_range\' must be greater 0!")
            return None
