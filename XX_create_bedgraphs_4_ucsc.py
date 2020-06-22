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


### Parse command line
######################

parser = argparse.ArgumentParser(description='Get distances to nearest TAD boundary for diegstes involved in directed interactions and undirected reference interactions not involved in directed interactions.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--enhanced-interaction-file', help='Enhanced interaction file supplemented with digest associated gene symbols and TSS as well as directionality P-values.', required=True)
parser.add_argument('--gopher-digest-file', help='File that contains the coordinates of all digests for a given genome build.', required=True)

args = parser.parse_args()
out_prefix = args.out_prefix
enhanced_interaction_file = args.enhanced_interaction_file
gopher_digest_file = args.gopher_digest_file

print("[INFO] " + "Input parameters")
print("\t[INFO] --out-prefix: " + out_prefix)
print("\t[INFO] --enhanced-interaction-file: " + enhanced_interaction_file)
print("\t[INFO] --gopher-digest-file: " + gopher_digest_file)

# Define digest class
class DigestRegion:

    # Attributes
    coordinates = None

    n_simple_read_pairs = 0
    n_twisted_read_pairs = 0

    interaction_distances = None
    read_pair_distances = None
    n_interaction = 0


    # Initializer
    def __init__(self, coordinates):
        self.coordinates = coordinates
        self.read_pair_distances = []
        self.interaction_distances = []

    def get_total_number_of_read_pairs(self):
        return self.n_simple_read_pairs + self.n_twisted_read_pairs

    def get_fraction_of_simple_read_pairs(self):
        if 0 < self.n_simple_read_pairs + self.n_twisted_read_pairs:
            return float(1.0 * self.n_simple_read_pairs / (self.n_simple_read_pairs + self.n_twisted_read_pairs))
        else:
            return 0.0

    def get_fraction_of_twisted_read_pairs(self):
        if 0 < self.n_simple_read_pairs + self.n_twisted_read_pairs:
            return float(1.0 * self.n_twisted_read_pairs / (self.n_simple_read_pairs + self.n_twisted_read_pairs))
        else:
            return 0.0

    def get_median_interaction_distance(self):
        return numpy.median(self.interaction_distances)

    def get_q1_read_pair_distance(self):
        return numpy.quantile(self.read_pair_distances, .25)

    def get_median_read_pair_distance(self):
        return numpy.median(self.read_pair_distances)

    def get_q3_read_pair_distance(self):
        return numpy.quantile(self.read_pair_distances, .75)

# Init digest sets
digests_from_dir_inter = dict()
digests_from_undir_inter = dict()

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

        # Add digest to set for directed interactions
        if interaction_category == "DIII" or interaction_category == "DI":
            digests_from_dir_inter[coord_key_da] = DigestRegion(coord_key_da)
            digests_from_dir_inter[coord_key_db] = DigestRegion(coord_key_db)

        # Add digest to set for undirected interactions
        if interaction_category == "UIRII" or interaction_category == "UIE":
            digests_from_undir_inter[coord_key_da] = DigestRegion(coord_key_da)
            digests_from_undir_inter[coord_key_db] = DigestRegion(coord_key_db)

        line = fp.readline()

    fp.close()

print("\tNumber of digests from directed interactions: " + str(len(digests_from_dir_inter)))
print("\tNumber of digests from undirected interactions: " + str(len(digests_from_undir_inter)))
print("\t\tIntersect: " + str(len(digests_from_dir_inter.keys() & digests_from_undir_inter.keys())))

print("[INFO] ... done.")


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

        if interaction_category == "DIII" or interaction_category == "DI":

            if coord_key_da in digests_from_dir_inter:
                digests_from_dir_inter[coord_key_da].n_interaction += 1
                digests_from_dir_inter[coord_key_da].n_simple_read_pairs += n_simple
                digests_from_dir_inter[coord_key_da].n_twisted_read_pairs += n_twisted
                digests_from_dir_inter[coord_key_da].interaction_distances.append(i_dist)
                for i in range(1,(n_simple + n_twisted)):
                    digests_from_dir_inter[coord_key_da].read_pair_distances.append(int(i_dist))

            if coord_key_db in digests_from_dir_inter:
                digests_from_dir_inter[coord_key_db].n_interaction += 1
                digests_from_dir_inter[coord_key_db].n_simple_read_pairs += n_simple
                digests_from_dir_inter[coord_key_db].n_twisted_read_pairs += n_twisted
                digests_from_dir_inter[coord_key_db].interaction_distances.append(int(i_dist))
                for i in range(1,(n_simple + n_twisted)):
                    digests_from_dir_inter[coord_key_db].read_pair_distances.append(int(i_dist))

        if interaction_category == "UIRII" or interaction_category == "UIE":

            if coord_key_da in digests_from_undir_inter:
                digests_from_undir_inter[coord_key_da].n_interaction += 1
                digests_from_undir_inter[coord_key_da].n_simple_read_pairs += n_simple
                digests_from_undir_inter[coord_key_da].n_twisted_read_pairs += n_twisted
                digests_from_undir_inter[coord_key_da].interaction_distances.append(int(i_dist))
                for i in range(1,(n_simple + n_twisted)):
                    digests_from_undir_inter[coord_key_da].read_pair_distances.append(int(i_dist))

            if coord_key_db in digests_from_undir_inter:
                digests_from_undir_inter[coord_key_db].n_interaction += 1
                digests_from_undir_inter[coord_key_db].n_simple_read_pairs += n_simple
                digests_from_undir_inter[coord_key_db].n_twisted_read_pairs += n_twisted
                digests_from_undir_inter[coord_key_db].interaction_distances.append(int(i_dist))
                for i in range(1,(n_simple + n_twisted)):
                    digests_from_undir_inter[coord_key_db].read_pair_distances.append(int(i_dist))

        line = fp.readline()

print("[INFO] ... done.")


# Create BedGraph tracks
########################

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
        if not coord_key in digests_from_undir_inter:
            total_rp_undir = 0
            frac_rp_simple_undir = 0.0
            frac_rp_twisted_undir = 0.0
        else:
            total_rp_undir = digests_from_undir_inter[coord_key].get_total_number_of_read_pairs()
            frac_rp_simple_undir = digests_from_undir_inter[coord_key].get_fraction_of_simple_read_pairs()
            frac_rp_twisted_undir = digests_from_undir_inter[coord_key].get_fraction_of_twisted_read_pairs()

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
        if not coord_key in digests_from_undir_inter:
            total_ia_undir = 0
            total_rp_undir = 0
            total_ia_rp_norm_undir = 0.0
        else:
            total_ia_undir = digests_from_undir_inter[coord_key].n_interaction
            total_rp_undir = digests_from_undir_inter[coord_key].get_total_number_of_read_pairs()
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

# Interaction distance
# --------------------

# Q1 interaction distance for directed and undirected interactions
bedgraph_stream_rp_dist_q1_dir_digests_output = open(out_prefix + "_rp_q1_dist_dir_inter_digests.bedgraph", 'wt')
bedgraph_stream_rp_dist_q1_dir_digests_output.write("track type=bedGraph name=\"Q1 read pair distance - Directed (DI and UIE - MK)\" description=\"Q1 read pair distance - Directed (DI and UIE - MK)\" visibility=full color=255,161,0 altColor=0,100,200 priority=20 maxHeightPixels=60" + '\n')

bedgraph_stream_rp_dist_q1_undir_digests_output = open(out_prefix + "_rp_q1_dist_undir_inter_digests.bedgraph", 'wt')
bedgraph_stream_rp_dist_q1_undir_digests_output.write("track type=bedGraph name=\"Q1 read pair distance - Undirected (DI and UIE - MK)\" description=\"Q1 read pair distance - Undirected (DI and UIE - MK)\" visibility=full color=0,100,200 altColor=0,100,200 priority=20 maxHeightPixels=60" + '\n')

# Median interaction distance for directed and undirected interactions
bedgraph_stream_rp_dist_median_dir_digests_output = open(out_prefix + "_rp_median_dist_dir_inter_digests.bedgraph", 'wt')
bedgraph_stream_rp_dist_median_dir_digests_output.write("track type=bedGraph name=\"Median read pair distance - Directed (DI and UIE - NCD4)\" description=\"Median read pair distance - Directed (DI and UIE - NCD4)\" visibility=full color=255,161,0 altColor=0,100,200 priority=20 maxHeightPixels=60" + '\n')

bedgraph_stream_rp_dist_median_undir_digests_output = open(out_prefix + "_rp_median_dist_undir_inter_digests.bedgraph", 'wt')
bedgraph_stream_rp_dist_median_undir_digests_output.write("track type=bedGraph name=\"Median read pair distance - Undirected (DI and UIE - NCD4)\" description=\"Median read pair distance - Undirected (DI and UIE - NCD4)\" visibility=full color=0,100,200 altColor=0,100,200 priority=20 maxHeightPixels=60" + '\n')

# Q2 interaction distance for directed and undirected interactions
bedgraph_stream_rp_dist_q3_dir_digests_output = open(out_prefix + "_rp_q3_dist_dir_inter_digests.bedgraph", 'wt')
bedgraph_stream_rp_dist_q3_dir_digests_output.write("track type=bedGraph name=\"Q3 read pair distance - Directed (DI and UIE - MK)\" description=\"Q3 read pair distance - Directed (DI and UIE - MK)\" visibility=full color=255,161,0 altColor=0,100,200 priority=20 maxHeightPixels=60" + '\n')

bedgraph_stream_rp_dist_q3_undir_digests_output = open(out_prefix + "_rp_q3_dist_undir_inter_digests.bedgraph", 'wt')
bedgraph_stream_rp_dist_q3_undir_digests_output.write("track type=bedGraph name=\"Q3 read pair distance - Undirected (DI and UIE - MK)\" description=\"Q3 read pair distance - Undirected (DI and UIE - MK)\" visibility=full color=0,100,200 altColor=0,100,200 priority=20 maxHeightPixels=60" + '\n')


print("[INFO] Iterating digest file to create BedGraph tracks for interaction distances ...")
with open(gopher_digest_file, 'rt') as fp:

    line = fp.readline()
    line = fp.readline()

    n_progress = 0
    while line:

        n_progress += 1
        if n_progress % 100000 == 0:
            print(n_progress)

        # Get key for current digest region
        fields = line.split('\t')
        coord_key = fields[0] + '\t' + fields[1] + '\t' + fields[2]

        # Get info for directed digests
        if not coord_key in digests_from_dir_inter:
            rp_dist_q1_dir = 0
            rp_dist_median_dir = 0
            rp_dist_q3_dir = 0
        else:
            rp_dist_q1_dir = digests_from_dir_inter[coord_key].get_q1_read_pair_distance()
            rp_dist_median_dir = digests_from_dir_inter[coord_key].get_median_read_pair_distance()
            rp_dist_q3_dir = digests_from_dir_inter[coord_key].get_q3_read_pair_distance()

        # Get info for undirected digests
        if not coord_key in digests_from_undir_inter:
            rp_dist_q1_undir = 0
            rp_dist_median_undir = 0
            rp_dist_q3_undir = 0
        else:
            rp_dist_q1_undir = digests_from_undir_inter[coord_key].get_q1_read_pair_distance()
            rp_dist_median_undir = digests_from_undir_inter[coord_key].get_median_read_pair_distance()
            rp_dist_q3_undir = digests_from_undir_inter[coord_key].get_q3_read_pair_distance()

        bedgraph_stream_rp_dist_q1_dir_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(rp_dist_q1_dir) + '\n')
        bedgraph_stream_rp_dist_q1_undir_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(-rp_dist_q1_undir) + '\n')

        bedgraph_stream_rp_dist_median_dir_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(rp_dist_median_dir) + '\n')
        bedgraph_stream_rp_dist_median_undir_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(-rp_dist_median_undir) + '\n')

        bedgraph_stream_rp_dist_q3_dir_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(rp_dist_q3_dir) + '\n')
        bedgraph_stream_rp_dist_q3_undir_digests_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(-rp_dist_q3_undir) + '\n')

        line = fp.readline()

bedgraph_stream_rp_dist_q1_dir_digests_output.close()
bedgraph_stream_rp_dist_q1_undir_digests_output.close()

bedgraph_stream_rp_dist_median_dir_digests_output.close()
bedgraph_stream_rp_dist_median_undir_digests_output.close()

bedgraph_stream_rp_dist_q3_dir_digests_output.close()
bedgraph_stream_rp_dist_q3_undir_digests_output.close()

print("[INFO] ... done.")