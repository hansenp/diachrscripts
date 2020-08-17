"""
This script takes as input an enhanced interaction file and determines two digest sets:

   1. Digests that are involved in at least one directed interaction.
   2. Digest that are involved in undirected reference interactions (and in no directed interaction).

The enhanced interaction file is then iterated a second time read pair numbers and interaction distances are collected.
For each digest of the two sets, we keep track of the following:

    # Strength

      1. Total number of read pairs that end in this digest.

      2. Total number of simple read pairs that end in this digest.

      3. Total number of twisted read pairs that in this digest.

    # Dispersion

       1. Total number of interactions that end in this digest.

       2. Total number of simple interactions that end in this digest.

       3. Total number of simple interactions that end in this digest.

   # Range

      1. Median interaction distance at the level of read pairs.

         a. Median distances of simple read pairs.

         b. Median distances of twisted read pairs.

         c. Maximum distance of simple read pairs.

         d. Maximum distance of simple read pairs.

         e. q1, q2, ...

For each feature, one BedGraph file will be created.

"""

import argparse
import gzip
import diachrscripts_toolkit as dclass
import numpy
import copy


### Parse command line
######################

parser = argparse.ArgumentParser(description='Get distances to nearest TAD boundary for diegstes involved in directed interactions and undirected reference interactions not involved in directed interactions.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--enhanced-interaction-file', help='Enhanced interaction file supplemented with digest associated gene symbols and TSS as well as directionality P-values.', required=True)
parser.add_argument('--gopher-digest-file', help='File that contains the coordinates of all digests for a given genome build.', required=True)
parser.add_argument('--data-set-track-tag', help='A tag that is used to create track names and descriptions.')
parser.add_argument('--digest-black-list', help='BED file with digests to be excluded.')
parser.add_argument('--invert-black-list', help='Use digests in black list only.', action='store_true', default=False)

args = parser.parse_args()
out_prefix = args.out_prefix
enhanced_interaction_file = args.enhanced_interaction_file
gopher_digest_file = args.gopher_digest_file
data_set_track_tag = args.data_set_track_tag
digest_black_list = args.digest_black_list
invert_black_list = args.invert_black_list

print("[INFO] " + "Input parameters")
print("\t[INFO] --out-prefix: " + out_prefix)
print("\t[INFO] --enhanced-interaction-file: " + enhanced_interaction_file)
print("\t[INFO] --gopher-digest-file: " + gopher_digest_file)
print("\t[INFO] --data-set-track-tag: " + data_set_track_tag)
if digest_black_list != None:
    print("\t[INFO] --digest-black-list: " + digest_black_list)
    print("\t[INFO] --invert-black-list: " + str(invert_black_list))
else:
    if invert_black_list:
        print("\t[ERROR] No black list given to be inverted!")
        exit(0)

# Define digest class
class DigestRegion:

    # Basic attributes
    chromosome = None
    start = None
    end = None
    length = None

    # Interaction numbers
    n_interaction_left = 0
    n_interaction_right = 0

    # Read numbers
    n_simple_left_read_pairs = 0
    n_twisted_left_read_pairs = 0
    n_simple_right_read_pairs = 0
    n_twisted_right_read_pairs = 0

    # Interaction distances
    interaction_distances = None

    # Read pair distances
    simple_left_read_pair_distances = None
    simple_right_read_pair_distances = None
    twisted_left_read_pair_distances = None
    twisted_right_read_pair_distances = None

    # Initializer
    def __init__(self, coordinates):

        # Digest coordinates
        self.chromosome = coordinates.split('\t')[0]
        self.start = int(coordinates.split('\t')[1])
        self.end = int(coordinates.split('\t')[2])

        # Median interaction distances
        self.simple_left_read_pair_distances = []
        self.simple_right_read_pair_distances = []
        self.twisted_left_read_pair_distances = []
        self.twisted_right_read_pair_distances = []

        self.interaction_distances = []

    def is_left(self, digest_coordinates):
        other_digest_chromosome = digest_coordinates.split('\t')[0]
        other_digest_start = int(digest_coordinates.split('\t')[1])
        if other_digest_chromosome != self.chromosome:
            print("[Warning] Left and right is not defined for trans-chromosomal interactions! Will return 'False'!")
            return False
        if other_digest_start < self.start:
            return True
        else:
            return False

    # Interaction numbers

    def get_total_number_of_interactions(self):
        return self.n_interaction_left + self.n_interaction_right

    def get_total_number_of_left_interactions(self):
        return self.n_interaction_left

    def get_total_number_of_right_interactions(self):
        return self.n_interaction_right

    # Read pair numbers

    def get_total_number_of_read_pairs(self):
        return self.n_simple_left_read_pairs + self.n_twisted_left_read_pairs + self.n_simple_right_read_pairs + self.n_twisted_right_read_pairs

    def get_total_number_of_left_read_pairs(self):
        return self.n_simple_left_read_pairs + self.n_twisted_left_read_pairs

    def get_total_number_of_right_read_pairs(self):
        return self.n_simple_right_read_pairs + self.n_twisted_right_read_pairs

    # Dispersion

    def get_total_dispersion(self):
        if 0 < self.n_simple_left_read_pairs + self.n_twisted_left_read_pairs + self.n_simple_right_read_pairs + self.n_twisted_right_read_pairs:
            return (self.n_interaction_left + self.n_interaction_right) / (self.n_simple_left_read_pairs + self.n_twisted_left_read_pairs + self.n_simple_right_read_pairs + self.n_twisted_right_read_pairs)
        else:
            return 0

    def get_left_dispersion(self):
        if 0 < self.n_simple_left_read_pairs + self.n_twisted_left_read_pairs:
            return self.n_interaction_left / (self.n_simple_left_read_pairs + self.n_twisted_left_read_pairs)
        else:
            return 0

    def get_right_dispersion(self):
        if 0 < self.n_simple_right_read_pairs + self.n_twisted_right_read_pairs:
            return self.n_interaction_right / (self.n_simple_right_read_pairs + self.n_twisted_right_read_pairs)
        else:
            return 0

    # Median interaction distances

    def get_median_interaction_distance(self):
        return numpy.median(self.interaction_distances)

    # Median read pair distances

    def get_median_read_pair_distance(self):
        return numpy.median(self.simple_left_read_pair_distances + self.twisted_left_read_pair_distances + self.simple_right_read_pair_distances + self.twisted_right_read_pair_distances)

    def get_median_simple_read_pair_distance(self):
        if 0 < len(self.simple_left_read_pair_distances + self.simple_right_read_pair_distances):
            return numpy.median(self.simple_left_read_pair_distances + self.simple_right_read_pair_distances)
        else:
            return 0

    def get_median_twisted_read_pair_distance(self):
        if 0 < len(self.twisted_left_read_pair_distances + self.twisted_right_read_pair_distances):
            return numpy.median(self.twisted_left_read_pair_distances + self.twisted_right_read_pair_distances)
        else:
            return 0

    def get_median_left_read_pair_distance(self):
        if 0 < len(self.simple_left_read_pair_distances + self.twisted_left_read_pair_distances):
            return numpy.median(self.simple_left_read_pair_distances + self.twisted_left_read_pair_distances)
        else:
            return 0

    def get_median_right_read_pair_distance(self):
        if 0 < len(self.simple_right_read_pair_distances + self.twisted_right_read_pair_distances):
            return numpy.median(self.simple_right_read_pair_distances + self.twisted_right_read_pair_distances)
        else:
            return 0

    def get_median_simple_left_read_pair_distance(self):
        return numpy.median(self.simple_left_read_pair_distances)

    def get_median_simple_right_read_pair_distance(self):
        return numpy.median(self.simple_right_read_pair_distances)

    def get_median_twisted_left_read_pair_distance(self):
        return numpy.median(self.twisted_left_read_pair_distances)

    def get_median_twisted_right_read_pair_distance(self):
        return numpy.median(self.twisted_right_read_pair_distances)

    def get_fraction_of_simple_read_pairs(self):
        if 0 < self.n_simple_left_read_pairs + self.n_twisted_left_read_pairs + self.n_simple_right_read_pairs + self.n_twisted_right_read_pairs:
            return float(1.0 * (self.n_simple_left_read_pairs + self.n_simple_right_read_pairs) / (self.n_simple_left_read_pairs + self.n_twisted_left_read_pairs + self.n_simple_right_read_pairs + self.n_twisted_right_read_pairs))
        else:
            return 0.0

    def get_fraction_of_twisted_read_pairs(self):
        if 0 < self.n_simple_left_read_pairs + self.n_twisted_left_read_pairs + self.n_simple_right_read_pairs + self.n_twisted_right_read_pairs:
            return float(1.0 * (self.n_twisted_left_read_pairs + self.n_twisted_right_read_pairs) / (self.n_simple_left_read_pairs + self.n_twisted_left_read_pairs + self.n_simple_right_read_pairs + self.n_twisted_right_read_pairs))
        else:
            return 0.0

# Init digest sets
digests_from_dir_inter = dict()
digests_from_uie_inter = dict()
digests_from_uii_inter = dict()
digest_to_be_excluded = set()

# Read digest regions to be excluded from file
if digest_black_list != None:
    with open(digest_black_list, 'rt') as fp:
        line = fp.readline().rstrip()
        while line:

            # Parse line
            chr, sta, end = line.split('\t')

            # Add to set
            digest_to_be_excluded.add(chr + '\t' + str(sta) + '\t' + str(end))

            # Go to next line
            line = fp.readline().rstrip()

# Determine sets of directed and undirected interacting digests
print("[INFO] Iterating enhanced interaction file and init digest objects with coordinates ...")
with gzip.open(enhanced_interaction_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()
    while line:

        # Parse enhanced interactions line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total, i_dist = \
            dclass.parse_enhanced_interaction_line_with_gene_symbols(line)

        coord_key_da = chr_a + '\t' + str(sta_a) + '\t' + str(end_a)
        coord_key_db = chr_b + '\t' + str(sta_b) + '\t' + str(end_b)

        if not invert_black_list:

            # Add digest to set for directed interactions
            if interaction_category == "DIII" or interaction_category == "DI":
                if coord_key_da not in digest_to_be_excluded:
                    digests_from_dir_inter[coord_key_da] = DigestRegion(coord_key_da)
                if coord_key_db not in digest_to_be_excluded:
                    digests_from_dir_inter[coord_key_db] = DigestRegion(coord_key_db)

            # Add digest to set for undirected interactions
            if interaction_category == "UIRII" or interaction_category == "UIE":
                if coord_key_da not in digest_to_be_excluded:
                    digests_from_uie_inter[coord_key_da] = DigestRegion(coord_key_da)
                if coord_key_db not in digest_to_be_excluded:
                    digests_from_uie_inter[coord_key_db] = DigestRegion(coord_key_db)

            # Add digest to set for undirected interactions
            if interaction_category == "UIRII" or interaction_category == "UII":
                if coord_key_da not in digest_to_be_excluded:
                    digests_from_uii_inter[coord_key_da] = DigestRegion(coord_key_da)
                if coord_key_db not in digest_to_be_excluded:
                    digests_from_uii_inter[coord_key_db] = DigestRegion(coord_key_db)

        else:

            # Add digest to set for directed interactions
            if interaction_category == "DIII" or interaction_category == "DI":
                if coord_key_da in digest_to_be_excluded:
                    digests_from_dir_inter[coord_key_da] = DigestRegion(coord_key_da)
                if coord_key_db in digest_to_be_excluded:
                    digests_from_dir_inter[coord_key_db] = DigestRegion(coord_key_db)

            # Add digest to set for undirected interactions
            if interaction_category == "UIRII" or interaction_category == "UIE":
                if coord_key_da in digest_to_be_excluded:
                    digests_from_uie_inter[coord_key_da] = DigestRegion(coord_key_da)
                if coord_key_db in digest_to_be_excluded:
                    digests_from_uie_inter[coord_key_db] = DigestRegion(coord_key_db)

            # Add digest to set for undirected interactions
            if interaction_category == "UIRII" or interaction_category == "UII":
                if coord_key_da in digest_to_be_excluded:
                    digests_from_uii_inter[coord_key_da] = DigestRegion(coord_key_da)
                if coord_key_db in digest_to_be_excluded:
                    digests_from_uii_inter[coord_key_db] = DigestRegion(coord_key_db)

        line = fp.readline()

    fp.close()

print("\tNumber of digests from DI interactions: " + str(len(digests_from_dir_inter)))
print("\tNumber of digests from UIE interactions: " + str(len(digests_from_uie_inter)))
print("\tNumber of digests from UII interactions: " + str(len(digests_from_uii_inter)))
print("\t\tIntersect of DI and UIE: " + str(len(digests_from_dir_inter.keys() & digests_from_uie_inter.keys())))
print("\t\tIntersect of DI and UII: " + str(len(digests_from_dir_inter.keys() & digests_from_uii_inter.keys())))
print("\t\tIntersect of UIE and UII: " + str(len(digests_from_uie_inter.keys() & digests_from_uii_inter.keys())))

print("[INFO] ... done.")

# Prepare three digest sets
digests_set_dir = dict()
digests_set_undir = dict()
digests_set_amb = dict()
digests_set_union = dict()

# Digests involved in directed interactions only
for dig_key in digests_from_dir_inter.keys() - (digests_from_uie_inter.keys() | digests_from_uii_inter.keys()):
    digests_set_dir[dig_key] = copy.deepcopy(digests_from_dir_inter[dig_key])

# Digests involved in undirected interactions only
for dig_key in (digests_from_uie_inter.keys() | digests_from_uii_inter.keys()) - digests_from_dir_inter.keys():
    if dig_key in digests_from_uie_inter:
        digests_set_undir[dig_key] = copy.deepcopy(digests_from_uie_inter[dig_key])
    if dig_key in digests_from_uii_inter:
        digests_set_undir[dig_key] = copy.deepcopy(digests_from_uii_inter[dig_key])

# Digests involved in directed and undirected interactions
for dig_key in digests_from_dir_inter.keys() & (digests_from_uie_inter.keys() | digests_from_uii_inter.keys()):
    digests_set_amb[dig_key] = copy.deepcopy(digests_from_dir_inter[dig_key])

# Digest that are involved in directed or undirected interactions
for dig_key in digests_set_dir.keys() | digests_set_undir.keys() | digests_set_amb.keys():
    if dig_key in digests_from_dir_inter:
        digests_set_union[dig_key] = copy.deepcopy(digests_from_dir_inter[dig_key])
    if dig_key in digests_set_undir:
        digests_set_union[dig_key] = copy.deepcopy(digests_set_undir[dig_key])
    if dig_key in digests_set_amb:
        digests_set_union[dig_key] = copy.deepcopy(digests_set_amb[dig_key])


print("Number of digests involved in directed interactions only: " + str(len(digests_set_dir)))
print("Number of digests involved in undirected interactions only: " + str(len(digests_set_undir)))
print("Number of digests involved in directed and undirected interactions: " + str(len(digests_set_amb)))
print("Number of digests involved in interactions (union): " + str(len(digests_set_union)))

digests_from_dir_inter = digests_set_dir
digests_from_uie_inter = digests_set_undir
digests_from_uii_inter = digests_set_amb
digests_from_inter = digests_set_union

# Collect information about digests
print("[INFO] Iterating enhanced interaction file to collect information about digests ...")
with gzip.open(enhanced_interaction_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()
    while line:

        # Parse enhanced interactions line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total, i_dist = \
            dclass.parse_enhanced_interaction_line_with_gene_symbols(line)

        coord_key_da = chr_a + '\t' + str(sta_a) + '\t' + str(end_a)
        coord_key_db = chr_b + '\t' + str(sta_b) + '\t' + str(end_b)

        x = line.split('\t')
        n_simple = int(x[4].split(':')[0])
        n_twisted = int(x[4].split(':')[1])

        # Union of digests


        # Directed interaction
        if interaction_category == "DI":

            # First digest of directed interaction
            if coord_key_da in digests_from_dir_inter:

                digests_from_dir_inter[coord_key_da].interaction_distances.append(i_dist)

                # Left or right
                if digests_from_dir_inter[coord_key_da].is_left(coord_key_db):

                    # Interaction numbers
                    digests_from_dir_inter[coord_key_da].n_interaction_left += 1

                    # Read pair numbers
                    digests_from_dir_inter[coord_key_da].n_simple_left_read_pairs += n_simple
                    digests_from_dir_inter[coord_key_da].n_twisted_left_read_pairs += n_twisted

                    # Read pair distances
                    for i in range(1,n_simple):
                        digests_from_dir_inter[coord_key_da].simple_left_read_pair_distances.append(int(i_dist))
                    for i in range(1, n_twisted):
                        digests_from_dir_inter[coord_key_da].twisted_left_read_pair_distances.append(int(i_dist))

                else:

                    # Interaction numbers
                    digests_from_dir_inter[coord_key_da].n_interaction_right += 1

                    # Read pair numbers
                    digests_from_dir_inter[coord_key_da].n_simple_right_read_pairs += n_simple
                    digests_from_dir_inter[coord_key_da].n_twisted_right_read_pairs += n_twisted

                    # Read pair distances
                    for i in range(1,n_simple):
                        digests_from_dir_inter[coord_key_da].simple_right_read_pair_distances.append(int(i_dist))
                    for i in range(1, n_twisted):
                        digests_from_dir_inter[coord_key_da].twisted_right_read_pair_distances.append(int(i_dist))

            # Second digest of directed interaction
            if coord_key_db in digests_from_dir_inter:

                digests_from_dir_inter[coord_key_db].interaction_distances.append(int(i_dist))

                # Left or right
                if digests_from_dir_inter[coord_key_db].is_left(coord_key_da):

                    # Interaction numbers
                    digests_from_dir_inter[coord_key_db].n_interaction_left += 1

                    # Read pair numbers
                    digests_from_dir_inter[coord_key_db].n_simple_left_read_pairs += n_simple
                    digests_from_dir_inter[coord_key_db].n_twisted_left_read_pairs += n_twisted

                    # Read pair distances
                    for i in range(1,n_simple):
                        digests_from_dir_inter[coord_key_db].simple_left_read_pair_distances.append(int(i_dist))
                    for i in range(1, n_twisted):
                        digests_from_dir_inter[coord_key_db].twisted_left_read_pair_distances.append(int(i_dist))
                else:

                    # Interaction numbers
                    digests_from_dir_inter[coord_key_db].n_interaction_right += 1

                    # Read pair numbers
                    digests_from_dir_inter[coord_key_db].n_simple_right_read_pairs += n_simple
                    digests_from_dir_inter[coord_key_db].n_twisted_right_read_pairs += n_twisted

                    # Read pair distances
                    for i in range(1,n_simple):
                        digests_from_dir_inter[coord_key_db].simple_right_read_pair_distances.append(int(i_dist))
                    for i in range(1, n_twisted):
                        digests_from_dir_inter[coord_key_db].twisted_right_read_pair_distances.append(int(i_dist))

        # Exclusive undirected interaction
        if interaction_category == "UIE":

            # First digest of interaction
            if coord_key_da in digests_from_uie_inter:

                digests_from_uie_inter[coord_key_da].interaction_distances.append(int(i_dist))

                # Left or right
                if digests_from_uie_inter[coord_key_da].is_left(coord_key_db):

                    # Interaction numbers
                    digests_from_uie_inter[coord_key_da].n_interaction_left += 1

                    # Read pair numbers
                    digests_from_uie_inter[coord_key_da].n_simple_left_read_pairs += n_simple
                    digests_from_uie_inter[coord_key_da].n_twisted_left_read_pairs += n_twisted

                    # Read pair distances
                    for i in range(1,n_simple):
                        digests_from_uie_inter[coord_key_da].simple_left_read_pair_distances.append(int(i_dist))
                    for i in range(1, n_twisted):
                        digests_from_uie_inter[coord_key_da].twisted_left_read_pair_distances.append(int(i_dist))
                else:

                    # Interaction numbers
                    digests_from_uie_inter[coord_key_da].n_interaction_right += 1

                    # Read pair numbers
                    digests_from_uie_inter[coord_key_da].n_simple_right_read_pairs += n_simple
                    digests_from_uie_inter[coord_key_da].n_twisted_right_read_pairs += n_twisted

                    # Read pair distances
                    for i in range(1,n_simple):
                        digests_from_uie_inter[coord_key_da].simple_right_read_pair_distances.append(int(i_dist))
                    for i in range(1, n_twisted):
                        digests_from_uie_inter[coord_key_da].twisted_right_read_pair_distances.append(int(i_dist))

            # Second digest of interaction
            if coord_key_db in digests_from_uie_inter:

                digests_from_uie_inter[coord_key_db].interaction_distances.append(int(i_dist))

                # Left or right
                if digests_from_uie_inter[coord_key_db].is_left(coord_key_da):

                    # Interaction numbers
                    digests_from_uie_inter[coord_key_db].n_interaction_left += 1

                    # Read pair numbers
                    digests_from_uie_inter[coord_key_db].n_simple_left_read_pairs += n_simple
                    digests_from_uie_inter[coord_key_db].n_twisted_left_read_pairs += n_twisted

                    # Read pair distances
                    for i in range(1,n_simple):
                        digests_from_uie_inter[coord_key_db].simple_left_read_pair_distances.append(int(i_dist))
                    for i in range(1, n_twisted):
                        digests_from_uie_inter[coord_key_db].twisted_left_read_pair_distances.append(int(i_dist))
                else:

                    # Interaction numbers
                    digests_from_uie_inter[coord_key_db].n_interaction_right += 1

                    # Read pair numbers
                    digests_from_uie_inter[coord_key_db].n_simple_right_read_pairs += n_simple
                    digests_from_uie_inter[coord_key_db].n_twisted_right_read_pairs += n_twisted

                    # Read pair distances
                    for i in range(1,n_simple):
                        digests_from_uie_inter[coord_key_db].simple_right_read_pair_distances.append(int(i_dist))
                    for i in range(1, n_twisted):
                        digests_from_uie_inter[coord_key_db].twisted_right_read_pair_distances.append(int(i_dist))

        # Inclusive undirected interaction
        if interaction_category == "UII":

            # First digest of interaction
            if coord_key_da in digests_from_uii_inter:

                digests_from_uii_inter[coord_key_da].interaction_distances.append(int(i_dist))

                # Left or right
                if digests_from_uii_inter[coord_key_da].is_left(coord_key_db):

                    # Interaction numbers
                    digests_from_uii_inter[coord_key_da].n_interaction_left += 1

                    # Read pair numbers
                    digests_from_uii_inter[coord_key_da].n_simple_left_read_pairs += n_simple
                    digests_from_uii_inter[coord_key_da].n_twisted_left_read_pairs += n_twisted

                    # Read pair distances
                    for i in range(1,n_simple):
                        digests_from_uii_inter[coord_key_da].simple_left_read_pair_distances.append(int(i_dist))
                    for i in range(1, n_twisted):
                        digests_from_uii_inter[coord_key_da].twisted_left_read_pair_distances.append(int(i_dist))
                else:

                    # Interaction numbers
                    digests_from_uii_inter[coord_key_da].n_interaction_right += 1

                    # Read pair numbers
                    digests_from_uii_inter[coord_key_da].n_simple_right_read_pairs += n_simple
                    digests_from_uii_inter[coord_key_da].n_twisted_right_read_pairs += n_twisted

                    # Read pair distances
                    for i in range(1,n_simple):
                        digests_from_uii_inter[coord_key_da].simple_right_read_pair_distances.append(int(i_dist))
                    for i in range(1, n_twisted):
                        digests_from_uii_inter[coord_key_da].twisted_right_read_pair_distances.append(int(i_dist))

            # Second digest of interaction
            if coord_key_db in digests_from_uii_inter:

                digests_from_uii_inter[coord_key_db].interaction_distances.append(int(i_dist))

                # Left or right
                if digests_from_uii_inter[coord_key_db].is_left(coord_key_da):

                    # Interaction numbers
                    digests_from_uii_inter[coord_key_db].n_interaction_left += 1

                    # Read pair numbers
                    digests_from_uii_inter[coord_key_db].n_simple_left_read_pairs += n_simple
                    digests_from_uii_inter[coord_key_db].n_twisted_left_read_pairs += n_twisted

                    # Read pair distances
                    for i in range(1,n_simple):
                        digests_from_uii_inter[coord_key_db].simple_left_read_pair_distances.append(int(i_dist))
                    for i in range(1, n_twisted):
                        digests_from_uii_inter[coord_key_db].twisted_left_read_pair_distances.append(int(i_dist))
                else:

                    # Interaction numbers
                    digests_from_uii_inter[coord_key_db].n_interaction_right += 1

                    # Read pair numbers
                    digests_from_uii_inter[coord_key_db].n_simple_right_read_pairs += n_simple
                    digests_from_uii_inter[coord_key_db].n_twisted_right_read_pairs += n_twisted

                    # Read pair distances
                    for i in range(1,n_simple):
                        digests_from_uii_inter[coord_key_db].simple_right_read_pair_distances.append(int(i_dist))
                    for i in range(1, n_twisted):
                        digests_from_uii_inter[coord_key_db].twisted_right_read_pair_distances.append(int(i_dist))

        line = fp.readline()

print("[INFO] ... done.")


# Create BedGraph tracks
########################                

# Interaction distance
# --------------------

# Median interaction distance for directed and undirected interactions
track_name = data_set_track_tag + " - MRPD - DI - S"
track_description = track_name
track_color = "255,161,0"
track_altColor ="0,100,200"
track_maxHeightPixels = 60
bedgraph_stream_simple_rp_dist_median_dir_digests_output = open(out_prefix + "_simple_rp_median_dist_di.bedgraph", 'wt')
bedgraph_stream_simple_rp_dist_median_dir_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_name = data_set_track_tag + " - MRPD - DI - T"
track_description = track_name
bedgraph_stream_twisted_rp_dist_median_dir_digests_output = open(out_prefix + "_twisted_rp_median_dist_di.bedgraph", 'wt')
bedgraph_stream_twisted_rp_dist_median_dir_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_name = data_set_track_tag + " - MRPD - UIE - S"
track_description = track_name
track_color = "0,100,200"
bedgraph_stream_simple_rp_dist_median_undir_digests_output = open(out_prefix + "_simple_rp_median_dist_uie.bedgraph", 'wt')
bedgraph_stream_simple_rp_dist_median_undir_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_name = data_set_track_tag + " - MRPD - UIE - T"
track_description = track_name
bedgraph_stream_twisted_rp_dist_median_undir_digests_output = open(out_prefix + "_twisted_rp_median_dist_uie.bedgraph", 'wt')
bedgraph_stream_twisted_rp_dist_median_undir_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')


track_name = data_set_track_tag + " - MRPD - UII - S"
track_description = track_name
track_color = "0,100,200"
bedgraph_stream_simple_rp_dist_median_uii_digests_output = open(out_prefix + "_simple_rp_median_dist_uii.bedgraph", 'wt')
bedgraph_stream_simple_rp_dist_median_uii_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_name = data_set_track_tag + " - MRPD - UII - T"
track_description = track_name
bedgraph_stream_twisted_rp_dist_median_uii_digests_output = open(out_prefix + "_twisted_rp_median_dist_uii.bedgraph", 'wt')
bedgraph_stream_twisted_rp_dist_median_uii_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

# DI
track_name = data_set_track_tag + " - MRPD - DI - L"
track_description = track_name
track_color = "255,161,0"
bedgraph_stream_left_rp_dist_median_di_digests_output = open(out_prefix + "_left_rp_median_dist_di.bedgraph", 'wt')
bedgraph_stream_left_rp_dist_median_di_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_name = data_set_track_tag + " - MRPD - DI - R"
track_description = track_name
bedgraph_stream_right_rp_dist_median_di_digests_output = open(out_prefix + "_right_rp_median_dist_di.bedgraph", 'wt')
bedgraph_stream_right_rp_dist_median_di_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

# UII
track_name = data_set_track_tag + " - MRPD - UII - L"
track_description = track_name
track_color = "0,100,200"
bedgraph_stream_left_rp_dist_median_uii_digests_output = open(out_prefix + "_left_rp_median_dist_uii.bedgraph", 'wt')
bedgraph_stream_left_rp_dist_median_uii_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_name = data_set_track_tag + " - MRPD - UII - R"
track_description = track_name
bedgraph_stream_right_rp_dist_median_uii_digests_output = open(out_prefix + "_right_rp_median_dist_uii.bedgraph", 'wt')
bedgraph_stream_right_rp_dist_median_uii_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_name = data_set_track_tag + " - MRPD - UII - LR"
track_description = track_name
bedgraph_stream_left_right_rp_dist_median_uii_digests_output = open(out_prefix + "_left_right_rp_median_dist_uii.bedgraph", 'wt')
bedgraph_stream_left_right_rp_dist_median_uii_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')


# UIE
track_name = data_set_track_tag + " - MRPD - UIE - L"
track_description = track_name
track_color = "0,100,200"
bedgraph_stream_left_rp_dist_median_uie_digests_output = open(out_prefix + "_left_rp_median_dist_uie.bedgraph", 'wt')
bedgraph_stream_left_rp_dist_median_uie_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_name = data_set_track_tag + " - MRPD - UIE - R"
track_description = track_name
bedgraph_stream_right_rp_dist_median_uie_digests_output = open(out_prefix + "_right_rp_median_dist_uie.bedgraph", 'wt')
bedgraph_stream_right_rp_dist_median_uie_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

# Read pair numbers

# DI
track_name = data_set_track_tag + " - RPN - DI - L"
track_description = track_name
track_color = "255,161,0"
bedgraph_stream_left_rp_number_di_digests_output = open(out_prefix + "_left_rp_number_di.bedgraph", 'wt')
bedgraph_stream_left_rp_number_di_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_name = data_set_track_tag + " - RPN - DI - R"
track_description = track_name
bedgraph_stream_right_rp_number_di_digests_output = open(out_prefix + "_right_rp_number_di.bedgraph", 'wt')
bedgraph_stream_right_rp_number_di_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

# UII
track_name = data_set_track_tag + " - RPN - UII - L"
track_description = track_name
track_color = "0,100,200"
bedgraph_stream_left_rp_number_uii_digests_output = open(out_prefix + "_left_rp_number_uii.bedgraph", 'wt')
bedgraph_stream_left_rp_number_uii_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_name = data_set_track_tag + " - RPN - UII - R"
track_description = track_name
bedgraph_stream_right_rp_number_uii_digests_output = open(out_prefix + "_right_rp_number_uii.bedgraph", 'wt')
bedgraph_stream_right_rp_number_uii_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

# UIE
track_name = data_set_track_tag + " - RPN - UIE - L"
track_description = track_name
track_color = "0,100,200"
bedgraph_stream_left_rp_number_uie_digests_output = open(out_prefix + "_left_rp_number_uie.bedgraph", 'wt')
bedgraph_stream_left_rp_number_uie_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_name = data_set_track_tag + " - RPN - UIE - R"
track_description = track_name
bedgraph_stream_right_rp_number_uie_digests_output = open(out_prefix + "_right_rp_number_uie.bedgraph", 'wt')
bedgraph_stream_right_rp_number_uie_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

# Dispersion

# DI
track_name = data_set_track_tag + " - DISP - DI - L"
track_description = track_name
track_color = "255,161,0"
bedgraph_stream_left_disp_di_digests_output = open(out_prefix + "_left_dispersion_di.bedgraph", 'wt')
bedgraph_stream_left_disp_di_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_name = data_set_track_tag + " - DISP - DI - R"
track_description = track_name
bedgraph_stream_right_disp_di_digests_output = open(out_prefix + "_right_dispersion_di.bedgraph", 'wt')
bedgraph_stream_right_disp_di_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

# UII
track_name = data_set_track_tag + " - DISP - UII - L"
track_description = track_name
track_color = "0,100,200"
bedgraph_stream_left_disp_uii_digests_output = open(out_prefix + "_left_dispersion_uii.bedgraph", 'wt')
bedgraph_stream_left_disp_uii_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_name = data_set_track_tag + " - DISP - UII - R"
track_description = track_name
bedgraph_stream_right_disp_uii_digests_output = open(out_prefix + "_right_dispersion_uii.bedgraph", 'wt')
bedgraph_stream_right_disp_uii_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

# UIE
track_name = data_set_track_tag + " - DISP - UIE - L"
track_description = track_name
track_color = "0,100,200"
bedgraph_stream_left_disp_uie_digests_output = open(out_prefix + "_left_dispersion_uie.bedgraph", 'wt')
bedgraph_stream_left_disp_uie_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_name = data_set_track_tag + " - DISP - UIE - R"
track_description = track_name
bedgraph_stream_right_disp_uie_digests_output = open(out_prefix + "_right_dispersion_uie.bedgraph", 'wt')
bedgraph_stream_right_disp_uie_digests_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')


print("[INFO] Iterating digest file to create BedGraph tracks for interaction distances ...")
with open(gopher_digest_file, 'rt') as fp:

    line = fp.readline()
    line = fp.readline()

    n_progress = 0
    while line:

        n_progress += 1
        if n_progress % 100000 == 0:
            print("\tProcessed " + str(n_progress) + " digests ...")

        # Get key for current digest region
        fields = line.split('\t')
        coord_key = fields[0] + '\t' + fields[1] + '\t' + fields[2]

        # Get info for directed digests
        if not coord_key in digests_from_dir_inter:
            simple_rp_dist_median_dir = 0
            twisted_rp_dist_median_dir = 0
            left_rp_dist_median_dir = 0
            right_rp_dist_median_dir = 0
            left_rp_number_dir = 0
            right_rp_number_dir = 0
            left_disp_dir = 0
            right_disp_dir = 0
        else:
            simple_rp_dist_median_dir = digests_from_dir_inter[coord_key].get_median_simple_read_pair_distance()
            twisted_rp_dist_median_dir = digests_from_dir_inter[coord_key].get_median_twisted_read_pair_distance()

            left_rp_dist_median_dir = digests_from_dir_inter[coord_key].get_median_left_read_pair_distance()
            right_rp_dist_median_dir = digests_from_dir_inter[coord_key].get_median_right_read_pair_distance()

            left_rp_number_dir = digests_from_dir_inter[coord_key].get_total_number_of_left_read_pairs()
            right_rp_number_dir = digests_from_dir_inter[coord_key].get_total_number_of_right_read_pairs()

            left_disp_dir = digests_from_dir_inter[coord_key].get_left_dispersion()
            right_disp_dir = digests_from_dir_inter[coord_key].get_right_dispersion()

        bedgraph_stream_simple_rp_dist_median_dir_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(simple_rp_dist_median_dir) + '\n')
        bedgraph_stream_twisted_rp_dist_median_dir_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(twisted_rp_dist_median_dir) + '\n')

        bedgraph_stream_left_rp_dist_median_di_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(left_rp_dist_median_dir) + '\n')
        bedgraph_stream_right_rp_dist_median_di_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(right_rp_dist_median_dir) + '\n')

        bedgraph_stream_left_rp_number_di_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(left_rp_number_dir) + '\n')
        bedgraph_stream_right_rp_number_di_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(right_rp_number_dir) + '\n')

        bedgraph_stream_left_disp_di_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(left_disp_dir) + '\n')
        bedgraph_stream_right_disp_di_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(right_disp_dir) + '\n')

        # Get info for UIE digests
        if not coord_key in digests_from_uie_inter:
            simple_rp_dist_median_undir = 0
            twisted_rp_dist_median_undir = 0
            left_rp_dist_median_undir = 0
            right_rp_dist_median_undir = 0
            left_rp_number_undir = 0
            right_rp_number_undir = 0
            left_disp_undir = 0
            right_disp_undir = 0
        else:
            simple_rp_dist_median_undir = digests_from_uie_inter[coord_key].get_median_simple_read_pair_distance()
            twisted_rp_dist_median_undir = digests_from_uie_inter[coord_key].get_median_twisted_read_pair_distance()
            left_rp_dist_median_undir = digests_from_uie_inter[coord_key].get_median_left_read_pair_distance()
            right_rp_dist_median_undir = digests_from_uie_inter[coord_key].get_median_right_read_pair_distance()
            left_rp_number_undir = digests_from_uie_inter[coord_key].get_total_number_of_left_read_pairs()
            right_rp_number_undir = digests_from_uie_inter[coord_key].get_total_number_of_right_read_pairs()
            left_disp_undir = digests_from_uie_inter[coord_key].get_left_dispersion()
            right_disp_undir = digests_from_uie_inter[coord_key].get_right_dispersion()

        bedgraph_stream_simple_rp_dist_median_undir_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(-simple_rp_dist_median_undir) + '\n')
        bedgraph_stream_twisted_rp_dist_median_undir_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(-twisted_rp_dist_median_undir) + '\n')

        bedgraph_stream_left_rp_dist_median_uie_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(-left_rp_dist_median_undir) + '\n')
        bedgraph_stream_right_rp_dist_median_uie_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(-right_rp_dist_median_undir) + '\n')

        bedgraph_stream_left_rp_number_uie_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(left_rp_number_undir) + '\n')
        bedgraph_stream_right_rp_number_uie_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(right_rp_number_undir) + '\n')

        bedgraph_stream_left_disp_uie_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(left_disp_undir) + '\n')
        bedgraph_stream_right_disp_uie_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(right_disp_undir) + '\n')

        # Get info for UII digests
        if not coord_key in digests_from_uii_inter:
            simple_rp_dist_median_uii = 0
            twisted_rp_dist_median_uii = 0
            left_rp_dist_median_uii = 0
            right_rp_dist_median_uii = 0
            left_rp_number_uii = 0
            right_rp_number_uii = 0
            left_disp_uii = 0
            right_disp_uii = 0
        else:
            simple_rp_dist_median_uii = digests_from_uii_inter[coord_key].get_median_simple_read_pair_distance()
            twisted_rp_dist_median_uii = digests_from_uii_inter[coord_key].get_median_twisted_read_pair_distance()
            left_rp_dist_median_uii = digests_from_uii_inter[coord_key].get_median_left_read_pair_distance()
            right_rp_dist_median_uii = digests_from_uii_inter[coord_key].get_median_right_read_pair_distance()
            left_rp_number_uii = digests_from_uii_inter[coord_key].get_total_number_of_left_read_pairs()
            right_rp_number_uii = digests_from_uii_inter[coord_key].get_total_number_of_right_read_pairs()
            left_disp_uii = digests_from_uii_inter[coord_key].get_left_dispersion()
            right_disp_uii = digests_from_uii_inter[coord_key].get_right_dispersion()

        bedgraph_stream_simple_rp_dist_median_uii_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(-simple_rp_dist_median_uii) + '\n')
        bedgraph_stream_twisted_rp_dist_median_uii_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(-twisted_rp_dist_median_uii) + '\n')

        bedgraph_stream_left_rp_dist_median_uii_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(-left_rp_dist_median_uii) + '\n')
        bedgraph_stream_right_rp_dist_median_uii_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(-right_rp_dist_median_uii) + '\n')
        bedgraph_stream_left_right_rp_dist_median_uii_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(min(left_rp_dist_median_uii,right_rp_dist_median_uii)) + '\n')

        bedgraph_stream_left_rp_number_uii_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(left_rp_number_uii) + '\n')
        bedgraph_stream_right_rp_number_uii_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(right_rp_number_uii) + '\n')

        bedgraph_stream_left_disp_uii_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(left_disp_uii) + '\n')
        bedgraph_stream_right_disp_uii_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(right_disp_uii) + '\n')

        line = fp.readline()

bedgraph_stream_simple_rp_dist_median_dir_digests_output.close()
bedgraph_stream_twisted_rp_dist_median_dir_digests_output.close()

bedgraph_stream_simple_rp_dist_median_undir_digests_output.close()
bedgraph_stream_twisted_rp_dist_median_undir_digests_output.close()

bedgraph_stream_simple_rp_dist_median_uii_digests_output.close()
bedgraph_stream_twisted_rp_dist_median_uii_digests_output.close()

bedgraph_stream_left_rp_dist_median_di_digests_output.close()
bedgraph_stream_right_rp_dist_median_di_digests_output.close()

bedgraph_stream_left_rp_dist_median_uii_digests_output.close()
bedgraph_stream_right_rp_dist_median_uii_digests_output.close()
bedgraph_stream_left_right_rp_dist_median_uii_digests_output.close()

bedgraph_stream_left_rp_dist_median_uie_digests_output.close()
bedgraph_stream_right_rp_dist_median_uie_digests_output.close()

bedgraph_stream_left_rp_number_di_digests_output.close()
bedgraph_stream_right_rp_number_di_digests_output.close()

bedgraph_stream_left_rp_number_uie_digests_output.close()
bedgraph_stream_right_rp_number_uie_digests_output.close()

bedgraph_stream_left_rp_number_uii_digests_output.close()
bedgraph_stream_right_rp_number_uii_digests_output.close()

bedgraph_stream_left_disp_di_digests_output.close()
bedgraph_stream_right_disp_di_digests_output.close()

bedgraph_stream_left_disp_uii_digests_output.close()
bedgraph_stream_right_disp_uii_digests_output.close()

bedgraph_stream_left_disp_uie_digests_output.close()
bedgraph_stream_right_disp_uie_digests_output.close()


print("[INFO] ... done.")

exit(0)

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################

# Interaction strength as measured by read pair numbers
# -----------------------------------------------------

# Difference of the total numbers of read pairs for directed and undirected for each digest
bedgraph_stream_total_rp_dir_minus_undir_inter_digests_output = open(out_prefix + "_total_rp_dir_minus_undir_inter_digests.bedgraph", 'wt')
bedgraph_stream_total_rp_dir_minus_undir_inter_digests_output.write("track type=bedGraph name=\"Read pairs - Directed minus undirected (DI and UIE - NCD4)\" description=\"Read pairs - Directed minus undirected (DI and UIE - NCD4)\" visibility=full color=255,161,0 altColor=0,100,200 priority=20 maxHeightPixels=60" + '\n')

# Difference of the fractions of read pairs for directed and undirected for each digest
bedgraph_stream_frac_rp_dir_minus_undir_inter_digests_output = open(out_prefix + "_frac_rp_dir_minus_undir_inter_digests.bedgraph", 'wt')
bedgraph_stream_frac_rp_dir_minus_undir_inter_digests_output.write("track type=bedGraph name=\"Read pairs - Directed minus undirected fraction (DI and UIE - NCD4)\" description=\"Read pairs - Directed minus undirected fraction (DI and UIE - NCD4)\" visibility=full color=255,161,0 altColor=0,100,200 priority=20 maxHeightPixels=30 viewLimits=-1.0:1.0" + '\n')

# Difference of the simple and twisted fractions for each digest
bedgraph_stream_frac_rp_dir_simple_minus_twisted_inter_digests_output = open(out_prefix + "_frac_rp_simple_minus_twisted_dir_inter_digests.bedgraph", 'wt')
bedgraph_stream_frac_rp_dir_simple_minus_twisted_inter_digests_output.write("track type=bedGraph name=\"Directed - Simple minus twisted fraction (DI and UIE - NCD4)\" description=\"Directed - Simple minus twisted fraction (DI and UIE - NCD4)\" visibility=full color=232,174,174 altColor=0,138,138 priority=20 maxHeightPixels=30 viewLimits=-1.0:1.0" + '\n')

# Difference of the simple and twisted fractions for each digest
bedgraph_stream_frac_rp_undir_simple_minus_twisted_inter_digests_output = open(out_prefix + "_frac_rp_simple_minus_twisted_undir_inter_digests.bedgraph", 'wt')
bedgraph_stream_frac_rp_undir_simple_minus_twisted_inter_digests_output.write("track type=bedGraph name=\"Undirected - Simple minus twisted fraction (DI and UIE - NCD4)\" description=\"Undirected - Simple minus twisted fraction (DI and UIE - NCD4)\" visibility=full color=232,174,174 altColor=0,138,138 priority=20 maxHeightPixels=30 viewLimits=-1.0:1.0" + '\n')


print("[INFO] Iterating digest file to create BedGraph tracks for read pair numbers ...")
with open(gopher_digest_file, 'rt') as fp:

    line = fp.readline()
    line = fp.readline()

    while line:

        # Get key for current digest region
        fields = line.split('\t')
        coord_key = fields[0] + '\t' + fields[1] + '\t' + fields[2]

        # Get info for directed digests
        if not coord_key in digests_from_dir_inter:
            total_rp_dir = 0
            frac_rp_simple_dir = 0.0
            frac_rp_twisted_dir = 0.0
        else:
            total_rp_dir = digests_from_dir_inter[coord_key].get_total_number_of_read_pairs()
            frac_rp_simple_dir = digests_from_dir_inter[coord_key].get_fraction_of_simple_read_pairs()
            frac_rp_twisted_dir = digests_from_dir_inter[coord_key].get_fraction_of_twisted_read_pairs()

        # Get info for undirected digests
        if not coord_key in digests_from_uie_inter:
            total_rp_undir = 0
            frac_rp_simple_undir = 0.0
            frac_rp_twisted_undir = 0.0
        else:
            total_rp_undir = digests_from_uie_inter[coord_key].get_total_number_of_read_pairs()
            frac_rp_simple_undir = digests_from_uie_inter[coord_key].get_fraction_of_simple_read_pairs()
            frac_rp_twisted_undir = digests_from_uie_inter[coord_key].get_fraction_of_twisted_read_pairs()

        # Difference between total numbers of read pairs from directed and undirected interactions
        total_rp_dir_minus_undir = total_rp_dir - total_rp_undir
        bedgraph_stream_total_rp_dir_minus_undir_inter_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(fields[2]) + '\t' + str(total_rp_dir_minus_undir) + '\n')

        # Difference between fractions of read pairs from directed and undirected interactions
        if 0 < total_rp_dir + total_rp_undir:
            frac_rp_dir = 1.0 * total_rp_dir / (total_rp_dir + total_rp_undir)
            frac_rp_undir = 1.0 * total_rp_undir / (total_rp_dir + total_rp_undir)
            frac_rp_dir_minus_undir = frac_rp_dir - frac_rp_undir
        else:
            frac_rp_dir_minus_undir = 0.0

        bedgraph_stream_frac_rp_dir_minus_undir_inter_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(fields[2]) + '\t' + str(frac_rp_dir_minus_undir) + '\n')

        # Difference between fractions of simple and twisted read pairs from directed interactions
        frac_rp_simple_minus_twisted_dir = frac_rp_simple_dir - frac_rp_twisted_dir
        bedgraph_stream_frac_rp_dir_simple_minus_twisted_inter_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(fields[2]) + '\t' + str(frac_rp_simple_minus_twisted_dir) + '\n')

        # Difference between fractions of simple and twisted read pairs from undirected interactions
        frac_rp_simple_minus_twisted_undir = frac_rp_simple_undir - frac_rp_twisted_undir
        bedgraph_stream_frac_rp_undir_simple_minus_twisted_inter_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(fields[2]) + '\t' + str(frac_rp_simple_minus_twisted_undir) + '\n')

        line = fp.readline()

fp.close()

bedgraph_stream_total_rp_dir_minus_undir_inter_digests_output.close()
bedgraph_stream_frac_rp_dir_minus_undir_inter_digests_output.close()
bedgraph_stream_frac_rp_dir_simple_minus_twisted_inter_digests_output.close()
bedgraph_stream_frac_rp_undir_simple_minus_twisted_inter_digests_output.close()


print("[INFO] ... done.")


# Interaction dispersion as measured by interaction numbers
# ---------------------------------------------------------

# Difference between total numbers of directed and undirected interactions
bedgraph_stream_total_ia_dir_minus_undir_inter_digests_output = open(out_prefix + "_total_ia_dir_minus_undir_inter_digests.bedgraph", 'wt')
bedgraph_stream_total_ia_dir_minus_undir_inter_digests_output.write("track type=bedGraph name=\"Interactions - Directed minus undirected (DI and UIE - NCD4)\" description=\"Interactions - Directed minus undirected (DI and UIE - NCD4)\" visibility=full color=255,161,0 altColor=0,100,200 priority=20 maxHeightPixels=60" + '\n')

# Difference between read pair normalized numbers of directed and undirected interactions
bedgraph_stream_total_ia_rp_norm_dir_minus_undir_inter_digests_output = open(out_prefix + "_total_ia_rp_norm_dir_minus_undir_inter_digests.bedgraph", 'wt')
bedgraph_stream_total_ia_rp_norm_dir_minus_undir_inter_digests_output.write("track type=bedGraph name=\"Interactions rp normalized - Directed minus undirected (DI and UIE - NCD4)\" description=\"Interactions rp normalized - Directed minus undirected (DI and UIE - NCD4)\" visibility=full color=255,161,0 altColor=0,100,200 priority=20 maxHeightPixels=60" + '\n')

print("[INFO] Iterating digest file to create BedGraph tracks for interaction numbers ...")
with open(gopher_digest_file, 'rt') as fp:

    line = fp.readline()
    line = fp.readline()

    while line:

        # Get key for current digest region
        fields = line.split('\t')
        coord_key = fields[0] + '\t' + fields[1] + '\t' + fields[2]

        # Get info for directed digests
        if not coord_key in digests_from_dir_inter:
            total_ia_dir = 0
            total_rp_dir = 0
            total_ia_rp_norm_dir = 0.0
        else:
            total_ia_dir = digests_from_dir_inter[coord_key].n_interaction
            total_rp_dir = digests_from_dir_inter[coord_key].get_total_number_of_read_pairs()
            total_ia_rp_norm_dir = float(1.0 * total_ia_dir / total_rp_dir)

        # Get info for undirected digests
        if not coord_key in digests_from_uie_inter:
            total_ia_undir = 0
            total_rp_undir = 0
            total_ia_rp_norm_undir = 0.0
        else:
            total_ia_undir = digests_from_uie_inter[coord_key].n_interaction
            total_rp_undir = digests_from_uie_inter[coord_key].get_total_number_of_read_pairs()
            total_ia_rp_norm_undir = float(1.0 * total_ia_undir / total_rp_undir)

        # Difference between total numbers of directed and undirected interactions
        total_ia_dir_minus_undir = total_ia_dir - total_ia_undir
        bedgraph_stream_total_ia_dir_minus_undir_inter_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(fields[2]) + '\t' + str(total_ia_dir_minus_undir) + '\n')

        # Difference between read pair normalized numbers of directed and undirected interactions
        total_ia_rp_norm_dir_minus_undir = total_ia_rp_norm_dir - total_ia_rp_norm_undir
        bedgraph_stream_total_ia_rp_norm_dir_minus_undir_inter_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(fields[2]) + '\t' + str(total_ia_rp_norm_dir_minus_undir) + '\n')

        line = fp.readline()

bedgraph_stream_total_ia_dir_minus_undir_inter_digests_output.close()
bedgraph_stream_total_ia_rp_norm_dir_minus_undir_inter_digests_output.close()

print("[INFO] ... done.")


