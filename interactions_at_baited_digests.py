"""
This script takes as input an enhanced interaction file and a BED file with baited digests and generates BED files
with regions that are spanned by interactions and end in baited digests.

The goal is to get a representation in the UCSC genome browser in which interactions that end in neighbouring baited
digests can be distinguished easily.

"""

import argparse
import gzip
import diachrscripts_toolkit as dclass
import random

parser = argparse.ArgumentParser(description='Get BED file with regions spanned by interactions that end in baited digests.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--enhanced-interaction-file', help='Enhanced interaction file supplemented with digest associated gene symbols and TSS as well as directionality P-values.', required=True)
parser.add_argument('--data-set-track-tag', help='A tag that is used to create track names and descriptions.')
parser.add_argument('--digest-black-list', help='BED file with digests to be excluded.')

args = parser.parse_args()
out_prefix = args.out_prefix
enhanced_interaction_file = args.enhanced_interaction_file
data_set_track_tag = args.data_set_track_tag
digest_black_list = args.digest_black_list

print("[INFO] " + "Input parameters")
print("\t[INFO] --out-prefix: " + out_prefix)
print("\t[INFO] --enhanced-interaction-file: " + enhanced_interaction_file)
print("\t[INFO] --data-set-track-tag: " + data_set_track_tag)
print("\t[INFO] --digest-black-list: " + digest_black_list)


# Define baited digest objects
class BaitedDigest:

    # Consecutive number for consecutive baited digests
    bait_id = None

    # Digest coordinates
    chr = None
    sta = None
    end = None

    # Lists of interactions ending in baited digest

    di_left_interactions = None
    di_right_interactions = None

    ui_left_interactions = None
    ui_right_interactions = None

    # Initializer
    def __init__(self, bait_id, chr, sta, end):

        self.bait_id = bait_id

        self.chr = chr
        self.sta = sta
        self.end = end

        self.di_left_interactions = []
        self.di_right_interactions = []

        self.ui_left_interactions = []
        self.ui_right_interactions = []

    def get_all_pairwise_differences_of_interaction_distances(self, ori_dir):

        if ori_dir  == "ui_left":
            interactions = self.ui_left_interactions
        if ori_dir  == "ui_right":
            interactions = self.ui_right_interactions
        if ori_dir  == "di_left":
            interactions = self.di_left_interactions
        if ori_dir  == "di_right":
            interactions = self.di_right_interactions

        # Get array of all interaction distances
        distances = []
        for interaction in interactions:
            chr, sta, end, x = interaction.split('\t')
            dist = int(end) - int(sta)
            distances.append(dist)

        # Sort distances by size
        distances_sorted = sorted(distances, reverse=True)

        # Calculate all pairwise distances
        dist_diff = []
        for i in range(0,len(distances_sorted)):
            for j in range(i+1, len(distances_sorted)):
                diff = distances_sorted[i] - distances_sorted[j]
                dist_diff.append(diff)
        return dist_diff


# Intit dictionary with baited digest objects
bait_id = 0
baited_digests = dict()
with open(digest_black_list, 'rt') as fp:
    line = fp.readline().rstrip()
    while line:

        # Parse line
        chr, sta, end = line.split('\t')

        # Add to dictionary
        baited_digests[chr + '\t' + str(sta) + '\t' + str(end)] = BaitedDigest(bait_id, chr, sta, end)
        bait_id += 1

        # Go to next line
        line = fp.readline().rstrip()

print("Initialized " + str(bait_id) + " baited digest objects.")

# Iterate interaction file and add interactions to baited digests
bed_stream_not_baited_interactions_output = open(out_prefix + "_not_baited_interactions.bed", 'wt')
bed_stream_not_baited_interactions_output.write("track type=bed name=\"" + "Not baited interactions" + "\" description=\"" + "Not baited interactions" + "\" visibility=2 itemRgb=\"On\"" + '\n')

print("[INFO] Iterating enhanced interaction file to collect information about digests ...")
with gzip.open(enhanced_interaction_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()
    while line:

        n_interaction_total += 1

        # Parse enhanced interactions line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total, i_dist = \
            dclass.parse_enhanced_interaction_line_with_gene_symbols(line)

        st_counts = line.split('\t')[4]
        simple_cnt = int(st_counts.split(":")[0])
        twisted_cnt = int(st_counts.split(":")[1])
        if simple_cnt < twisted_cnt:
            st_counts = "T|" + st_counts
        else:
            st_counts = "S|" + st_counts

        # Get digest coordinates from enhanced interaction file
        coord_key_da = chr_a + '\t' + str(sta_a) + '\t' + str(end_a)
        coord_key_db = chr_b + '\t' + str(sta_b) + '\t' + str(end_b)

        if interaction_category == 'DI':

            #if coord_key_da in baited_digests: # right
            if enrichment_pair_tag == 'EN':  # right
                baited_digests[coord_key_da].di_right_interactions.append(chr_a + '\t' + str(sta_a) + '\t' + str(end_b) + '\t' + st_counts)

            #if coord_key_db in baited_digests: # left
            if enrichment_pair_tag == 'NE':  # left
                baited_digests[coord_key_db].di_left_interactions.append(chr_a + '\t' + str(sta_a) + '\t' + str(end_b) + '\t' + st_counts)

        else:

            #if coord_key_da in baited_digests: # right
            if enrichment_pair_tag == 'EN':  # right
                baited_digests[coord_key_da].ui_right_interactions.append(chr_a + '\t' + str(sta_a) + '\t' + str(end_b) + '\t' + st_counts)

            #if coord_key_db in baited_digests: # left
            if enrichment_pair_tag == 'NE':  # left
                baited_digests[coord_key_db].ui_left_interactions.append(chr_a + '\t' + str(sta_a) + '\t' + str(end_b) + '\t' + st_counts)

        if coord_key_da not in baited_digests and coord_key_db not in baited_digests:
            bed_stream_not_baited_interactions_output.write(chr_a + '\t' + str(sta_a) + '\t' + str(end_b) + '\n')

        line = fp.readline()

bed_stream_not_baited_interactions_output.close()

# Sort dictionary by bait ID
baited_digests_bid = dict()
for key in baited_digests:
    baited_digests_bid[baited_digests[key].bait_id] = baited_digests[key]


# Write regions spanned by interactions to BED files
# chr2:26,547,643-31,761,880
bed_stream_baited_interactions_1_di_output = open(out_prefix + "_baited_interactions_1_di.bed", 'wt')
bed_stream_baited_interactions_1_di_output.write("track type=bed name=\"" + "Baited interactions 1 - DI" + "\" description=\"" + "Baited interactions 1 - DI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

bed_stream_baited_interactions_2_di_output = open(out_prefix + "_baited_interactions_2_di.bed", 'wt')
bed_stream_baited_interactions_2_di_output.write("track type=bed name=\"" + "Baited interactions 2 - DI" + "\" description=\"" + "Baited interactions 2 - DI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

bed_stream_baited_interactions_3_di_output = open(out_prefix + "_baited_interactions_3_di.bed", 'wt')
bed_stream_baited_interactions_3_di_output.write("track type=bed name=\"" + "Baited interactions 3 - DI" + "\" description=\"" + "Baited interactions 3 - DI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

bed_stream_baited_interactions_4_di_output = open(out_prefix + "_baited_interactions_4_di.bed", 'wt')
bed_stream_baited_interactions_4_di_output.write("track type=bed name=\"" + "Baited interactions 4 - DI" + "\" description=\"" + "Baited interactions 4 - DI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

bed_stream_baited_interactions_5_di_output = open(out_prefix + "_baited_interactions_5_di.bed", 'wt')
bed_stream_baited_interactions_5_di_output.write("track type=bed name=\"" + "Baited interactions 5 - DI" + "\" description=\"" + "Baited interactions 5 - DI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

bed_stream_baited_interactions_6_di_output = open(out_prefix + "_baited_interactions_6_di.bed", 'wt')
bed_stream_baited_interactions_6_di_output.write("track type=bed name=\"" + "Baited interactions 6 - DI" + "\" description=\"" + "Baited interactions 6 - DI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

bed_stream_baited_interactions_7_di_output = open(out_prefix + "_baited_interactions_7_di.bed", 'wt')
bed_stream_baited_interactions_7_di_output.write("track type=bed name=\"" + "Baited interactions 7 - DI" + "\" description=\"" + "Baited interactions 7 - DI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

bed_stream_baited_interactions_8_di_output = open(out_prefix + "_baited_interactions_8_di.bed", 'wt')
bed_stream_baited_interactions_8_di_output.write("track type=bed name=\"" + "Baited interactions 8 - DI" + "\" description=\"" + "Baited interactions 8 - DI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

bed_stream_baited_interactions_9_di_output = open(out_prefix + "_baited_interactions_9_di.bed", 'wt')
bed_stream_baited_interactions_9_di_output.write("track type=bed name=\"" + "Baited interactions 9 - DI" + "\" description=\"" + "Baited interactions 9 - DI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

bed_stream_baited_interactions_10_di_output = open(out_prefix + "_baited_interactions_10_di.bed", 'wt')
bed_stream_baited_interactions_10_di_output.write("track type=bed name=\"" + "Baited interactions 10 - DI" + "\" description=\"" + "Baited interactions 10 - DI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

# XXXX
bed_stream_array_di = []

bed_stream_baited_interactions_1_ui_output = open(out_prefix + "_baited_interactions_1_ui.bed", 'wt')
bed_stream_baited_interactions_1_ui_output.write("track type=bed name=\"" + "Baited interactions 1 - UI" + "\" description=\"" + "Baited interactions 1 - UI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

bed_stream_baited_interactions_2_ui_output = open(out_prefix + "_baited_interactions_2_ui.bed", 'wt')
bed_stream_baited_interactions_2_ui_output.write("track type=bed name=\"" + "Baited interactions 2 - UI" + "\" description=\"" + "Baited interactions 2 - UI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

bed_stream_baited_interactions_3_ui_output = open(out_prefix + "_baited_interactions_3_ui.bed", 'wt')
bed_stream_baited_interactions_3_ui_output.write("track type=bed name=\"" + "Baited interactions 3 - UI" + "\" description=\"" + "Baited interactions 3 - UI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

bed_stream_baited_interactions_4_ui_output = open(out_prefix + "_baited_interactions_4_ui.bed", 'wt')
bed_stream_baited_interactions_4_ui_output.write("track type=bed name=\"" + "Baited interactions 4 - UI" + "\" description=\"" + "Baited interactions 4 - UI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

bed_stream_baited_interactions_5_ui_output = open(out_prefix + "_baited_interactions_5_ui.bed", 'wt')
bed_stream_baited_interactions_5_ui_output.write("track type=bed name=\"" + "Baited interactions 5 - UI" + "\" description=\"" + "Baited interactions 5 - UI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

bed_stream_baited_interactions_6_ui_output = open(out_prefix + "_baited_interactions_6_ui.bed", 'wt')
bed_stream_baited_interactions_6_ui_output.write("track type=bed name=\"" + "Baited interactions 6 - UI" + "\" description=\"" + "Baited interactions 6 - UI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

bed_stream_baited_interactions_7_ui_output = open(out_prefix + "_baited_interactions_7_ui.bed", 'wt')
bed_stream_baited_interactions_7_ui_output.write("track type=bed name=\"" + "Baited interactions 7 - UI" + "\" description=\"" + "Baited interactions 7 - UI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

bed_stream_baited_interactions_8_ui_output = open(out_prefix + "_baited_interactions_8_ui.bed", 'wt')
bed_stream_baited_interactions_8_ui_output.write("track type=bed name=\"" + "Baited interactions 8 - UI" + "\" description=\"" + "Baited interactions 8 - UI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

bed_stream_baited_interactions_9_ui_output = open(out_prefix + "_baited_interactions_9_ui.bed", 'wt')
bed_stream_baited_interactions_9_ui_output.write("track type=bed name=\"" + "Baited interactions 9 - UI" + "\" description=\"" + "Baited interactions 9 - UI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

bed_stream_baited_interactions_10_ui_output = open(out_prefix + "_baited_interactions_10_ui.bed", 'wt')
bed_stream_baited_interactions_10_ui_output.write("track type=bed name=\"" + "Baited interactions 10 - UI" + "\" description=\"" + "Baited interactions 10 - UI" + "\" visibility=2 itemRgb=\"On\"" + '\n')

dist_diff_ui_left = []
dist_diff_ui_right = []
dist_diff_di_left = []
dist_diff_di_right = []
bed_file_number = 1
for key in sorted(baited_digests_bid):

    # Generate random numbers
    r_number_r = random.randint(155, 255) # 100 colors
    r_number_b = random.randint(0, 100)

    # Directed interactions right
    red_color = r_number_r
    green_color = 355 - r_number_r - r_number_b
    blue_color = r_number_b

    dist_diff_ui_left.extend(baited_digests_bid[key].get_all_pairwise_differences_of_interaction_distances("ui_left"))
    dist_diff_ui_right.extend(baited_digests_bid[key].get_all_pairwise_differences_of_interaction_distances("ui_right"))
    dist_diff_di_left.extend(baited_digests_bid[key].get_all_pairwise_differences_of_interaction_distances("di_left"))
    dist_diff_di_right.extend(baited_digests_bid[key].get_all_pairwise_differences_of_interaction_distances("di_right"))

    strand = "+"
    for i_coords in baited_digests_bid[key].di_right_interactions:

        i_chr, i_sta, i_end, i_st = i_coords.split('\t')

        if bed_file_number == 1:
            bed_stream_baited_interactions_1_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 2:
            bed_stream_baited_interactions_2_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 3:
            bed_stream_baited_interactions_3_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 4:
            bed_stream_baited_interactions_4_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 5:
            bed_stream_baited_interactions_5_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 6:
            bed_stream_baited_interactions_6_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 7:
            bed_stream_baited_interactions_7_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 8:
            bed_stream_baited_interactions_8_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 9:
            bed_stream_baited_interactions_9_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 10:
            bed_stream_baited_interactions_10_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

    # Directed interactions left
    red_color = r_number_b
    green_color = 355 - r_number_r - r_number_b
    blue_color = r_number_r
    strand = "-"
    for i_coords in baited_digests_bid[key].di_left_interactions:

        i_chr, i_sta, i_end, i_st = i_coords.split('\t')

        if bed_file_number == 1:
            bed_stream_baited_interactions_1_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 2:
            bed_stream_baited_interactions_2_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 3:
            bed_stream_baited_interactions_3_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 4:
            bed_stream_baited_interactions_4_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 5:
            bed_stream_baited_interactions_5_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 6:
            bed_stream_baited_interactions_6_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 7:
            bed_stream_baited_interactions_7_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 8:
            bed_stream_baited_interactions_8_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 9:
            bed_stream_baited_interactions_9_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 10:
            bed_stream_baited_interactions_10_di_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

    # Undirected interactions right
    red_color = r_number_r
    green_color = 355 - r_number_r - r_number_b
    blue_color = r_number_b
    strand = "+"
    for i_coords in baited_digests_bid[key].ui_right_interactions:

        i_chr, i_sta, i_end, i_st = i_coords.split('\t')

        if bed_file_number == 1:
            bed_stream_baited_interactions_1_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 2:
            bed_stream_baited_interactions_2_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 3:
            bed_stream_baited_interactions_3_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 4:
            bed_stream_baited_interactions_4_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 5:
            bed_stream_baited_interactions_5_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 6:
            bed_stream_baited_interactions_6_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 7:
            bed_stream_baited_interactions_7_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 8:
            bed_stream_baited_interactions_8_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 9:
            bed_stream_baited_interactions_9_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 10:
            bed_stream_baited_interactions_10_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

    # Undirected interactions left
    red_color = r_number_b
    green_color = 355 - r_number_r - r_number_b
    blue_color = r_number_r
    strand = "-"
    for i_coords in baited_digests_bid[key].ui_left_interactions:

        i_chr, i_sta, i_end, i_st = i_coords.split('\t')

        if bed_file_number == 1:
            bed_stream_baited_interactions_1_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 2:
            bed_stream_baited_interactions_2_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 3:
            bed_stream_baited_interactions_3_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 4:
            bed_stream_baited_interactions_4_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 5:
            bed_stream_baited_interactions_5_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 6:
            bed_stream_baited_interactions_6_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 7:
            bed_stream_baited_interactions_7_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 8:
            bed_stream_baited_interactions_8_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 9:
            bed_stream_baited_interactions_9_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

        if bed_file_number == 10:
            bed_stream_baited_interactions_10_ui_output.write(i_coords + "|" + str(baited_digests_bid[key].bait_id) + "\t0\t" + strand + "\t" + str(i_sta) + '\t' + str(i_end) + "\t" + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

    bed_file_number +=1
    if bed_file_number == 11:
        bed_file_number = 1

bed_stream_baited_interactions_1_di_output.close()
bed_stream_baited_interactions_2_di_output.close()
bed_stream_baited_interactions_3_di_output.close()
bed_stream_baited_interactions_4_di_output.close()
bed_stream_baited_interactions_5_di_output.close()
bed_stream_baited_interactions_6_di_output.close()
bed_stream_baited_interactions_7_di_output.close()
bed_stream_baited_interactions_8_di_output.close()
bed_stream_baited_interactions_9_di_output.close()
bed_stream_baited_interactions_10_di_output.close()

bed_stream_baited_interactions_1_ui_output.close()
bed_stream_baited_interactions_2_ui_output.close()
bed_stream_baited_interactions_3_ui_output.close()
bed_stream_baited_interactions_4_ui_output.close()
bed_stream_baited_interactions_5_ui_output.close()
bed_stream_baited_interactions_6_ui_output.close()
bed_stream_baited_interactions_7_ui_output.close()
bed_stream_baited_interactions_8_ui_output.close()
bed_stream_baited_interactions_9_ui_output.close()
bed_stream_baited_interactions_10_ui_output.close()


n_di = 0
n_ui = 0
n_di_total = 0
n_ui_total = 0
for key in baited_digests:
    n_di = len(baited_digests[key].di_left_interactions) + len(baited_digests[key].di_right_interactions)
    n_ui = len(baited_digests[key].ui_left_interactions) + len(baited_digests[key].ui_right_interactions)
    n_di_total = n_di_total + n_di
    n_ui_total = n_ui_total + n_ui
    baited_digests_bid[baited_digests[key].bait_id] = baited_digests[key]
    if 10 < (n_ui+1)/(n_di+1) and 20<n_ui:
        print("------")
        print(baited_digests[key].chr + "\t" + str(baited_digests[key].sta) + "\t" + str(baited_digests[key].end))
        print(str(len(baited_digests[key].di_left_interactions)) + "\t" + str(len(baited_digests[key].di_right_interactions)))
        print(str(len(baited_digests[key].ui_left_interactions)) + "\t" + str(len(baited_digests[key].ui_right_interactions)))

print("Number of Baits: " + str(bait_id))
print("n_di_total: " + str(n_di_total))
print("n_ui_total: " + str(n_ui_total))
print("n_di_total+n_ui_total: " + str(n_di_total+n_ui_total))
print("n_interaction_total: " + str(n_interaction_total))
# chr5:156,927,010-160,559,341
#print(n_di_total/(n_di_total+n_ui_total))

tab_stream_pairwise_dist_diff_ui_left_output = open(out_prefix + "_pairwise_dist_diff_ui_left.tab", 'wt')
for dist_diff in sorted(dist_diff_ui_left):
    tab_stream_pairwise_dist_diff_ui_left_output.write(str(dist_diff) + '\n')
tab_stream_pairwise_dist_diff_ui_left_output.close()

tab_stream_pairwise_dist_diff_ui_right_output = open(out_prefix + "_pairwise_dist_diff_ui_right.tab", 'wt')
for dist_diff in sorted(dist_diff_ui_right):
    tab_stream_pairwise_dist_diff_ui_right_output.write(str(dist_diff) + '\n')
tab_stream_pairwise_dist_diff_ui_right_output.close()

tab_stream_pairwise_dist_diff_di_left_output = open(out_prefix + "_pairwise_dist_diff_di_left.tab", 'wt')
for dist_diff in sorted(dist_diff_di_left):
    tab_stream_pairwise_dist_diff_di_left_output.write(str(dist_diff) + '\n')
tab_stream_pairwise_dist_diff_di_left_output.close()

tab_stream_pairwise_dist_diff_di_right_output = open(out_prefix + "_pairwise_dist_diff_di_right.tab", 'wt')
for dist_diff in sorted(dist_diff_di_right):
    tab_stream_pairwise_dist_diff_di_right_output.write(str(dist_diff) + '\n')
tab_stream_pairwise_dist_diff_di_right_output.close()
