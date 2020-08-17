"""
This script takes as input an enhanced interaction file and determines the densities of different digest types along
the genome. First the interaction file is iterated and two digest sets are determined:

   1. Digests that are involved in at least one directed interaction.
   2. Digests that are involved in at least one undirected interaction.

Because a digest can be involved in more than one interaction, these two digest sets overlap. The following subsets and
union sets can be derived from this:

   1. Digests that are involved in directed interaction only.
   2. Digests that are involved in directed and undirected interactions (intersect).
   3. Digests that are involved in undirected interaction only.
   4. Digests that are involved in at least one directed interaction (same set as above).
   5. Digests that are involved in at least one undirected interaction (same set as above).
   6. Digests that are involved in directed or undirected interactions (union).
   7. Digests that are not involved in any interaction.

Optionally, a blacklist with coordinates of digests that should be excluded from the analysis can be passed to the
script (8). It is also determined how the various digest types are affected by the blacklisting (for 1, 2 and 3).

A digest file must be passed to the script, which contains the sorted coordinates of all digests in the genome
(created with GOPHER). This file is iterated and the density of 5 different digest types is determined within a sliding
window of x digests (by default x = 25):

   1. DDAD: Digests that are involved in at least one directed interaction (4).
   2. UD: Digests that are involved in undirected interactions only (3).
   3. CD: Digests that are involved in directed or undirected interactions (6).
   4. ND: Digests that are not involved in any interaction (7).
   5. BL: Blacklisted digests (8).

The density for a specific digest type in a window at a specific position (digest) is defined as the number of digests
of this type divided by the total number of digests in the window, which is constant (default is 25).

Since the digests are of different sizes, the width of the window also varies. Therefore, the density is normalized to
the current width of the window. The normalized density is assigned to the digest in the middle of the window (e.g. the
13th digest with a window width of 25 digests) and written to a BedGraph file.

The BedGraph file can be displayed in UCSC's Gennom browser and the values from the fourth column can be used to
calculate correlations between the various digest types.

"""


import argparse
import gzip
import diachrscripts_toolkit as dclass


### Parse command line
######################

parser = argparse.ArgumentParser(description='Determine densities along the genome.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--enhanced-interaction-file', help='Enhanced interaction file with interactions already divided into directional (DI, DIAA, DIAI or DIII) and undirected interactions (all others) based on directionality P-value threshold.', required=True)
parser.add_argument('--gopher-digest-file', help='File produced by GOPHER that contains the coordinates of all digests for a given genome build.', required=True)
parser.add_argument('--data-set-track-tag', help='A tag that is used to create track names and descriptions.')
parser.add_argument('--digest-black-list', help='BED file with digests to be excluded.')
parser.add_argument('--invert-black-list', help='Use digests in blacklist only.', action='store_true', default=False)
parser.add_argument('--digest-density-range', help='Number of consecutive digests in the sliding window.', default=25)

args = parser.parse_args()
out_prefix = args.out_prefix
enhanced_interaction_file = args.enhanced_interaction_file
gopher_digest_file = args.gopher_digest_file
data_set_track_tag = args.data_set_track_tag
digest_black_list = args.digest_black_list
invert_black_list = args.invert_black_list
digest_density_range = int(args.digest_density_range)

print("[INFO] " + "Input parameters")
print("\t[INFO] --out-prefix: " + out_prefix)
print("\t[INFO] --enhanced-interaction-file: " + enhanced_interaction_file)
print("\t[INFO] --gopher-digest-file: " + gopher_digest_file)
print("\t[INFO] --digest-density-range: " + str(digest_density_range))
print("\t[INFO] --data-set-track-tag: " + data_set_track_tag)
if digest_black_list != None:
    print("\t[INFO] --digest-black-list: " + digest_black_list)
    print("\t[INFO] --invert-black-list: " + str(invert_black_list))
else:
    if invert_black_list:
        print("\t[ERROR] No blacklist given to be inverted!")
        exit(0)

# Tags for directed interaction categories
directed_interaction_category_tags = ["DI","DIII","DIAI","DIAA"]

# Read digest regions to be excluded from file
digest_to_be_excluded = set()
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


### Determine digest types
##########################

print("[INFO] Iterating enhanced interaction to determine digest types ...")

# Count directed interactions and directed interactions with at least one baited digest
di_cnt = 0
di_bl_cnt = 0

# Count undirected interactions and undirected interactions with at least one baited digest
ui_cnt = 0
ui_bl_cnt = 0

# Collect digest coordinates before blacklisting of baited digests
ddad_key_set = set()
udad_key_set = set()

with gzip.open(enhanced_interaction_file, 'rt') as fp:

    line = fp.readline()
    n_progress = 0
    while line:

        # Report progress
        n_progress += 1
        if n_progress % 1000000 == 0:
            print("\tProcessed " + str(n_progress) + " interactions ...")

        # Parse enhanced interaction line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total, i_dist = \
            dclass.parse_enhanced_interaction_line_with_gene_symbols(line)

        # Get coordinates of associated digests
        coord_key_da = chr_a + '\t' + str(sta_a) + '\t' + str(end_a)
        coord_key_db = chr_b + '\t' + str(sta_b) + '\t' + str(end_b)

        # Add associated digests to sets
        if interaction_category in directed_interaction_category_tags:
            ddad_key_set.add(coord_key_da)
            ddad_key_set.add(coord_key_db)
            di_cnt += 1
            if coord_key_da in digest_to_be_excluded or coord_key_db in digest_to_be_excluded:
                di_bl_cnt += 1
        else:
            udad_key_set.add(coord_key_da)
            udad_key_set.add(coord_key_db)
            ui_cnt += 1
            if coord_key_da in digest_to_be_excluded or coord_key_db in digest_to_be_excluded:
                ui_bl_cnt += 1

        line = fp.readline()

    fp.close()

print("[INFO] ... done.")

# Derive subsets into three disjoint subsets
dd_key_set = ddad_key_set - udad_key_set
ad_key_set = ddad_key_set & udad_key_set
ud_key_set = udad_key_set - ddad_key_set
cd_key_set = ddad_key_set | udad_key_set

# Print out some statistics
print("[INFO] Statistics about interactions and digest types ...")
print("\t[INFO] Number of all DI interactions: " + str(di_cnt))
print("\t[INFO] Number of all UI interactions: " + str(ui_cnt))
print("\t[INFO] ------------------------------------------------------- " )
print("\t[INFO] Number of DI interactions with baited digest: " + str(di_bl_cnt) + " (" + "{:.2%}".format(di_bl_cnt/di_cnt) + ")")
print("\t[INFO] Number of UI interactions with baited digest: " + str(ui_bl_cnt) + " (" + "{:.2%}".format(ui_bl_cnt/ui_cnt) + ")")
print("\t[INFO] ------------------------------------------------------- " )
print("\t[INFO] Number of DD digests before blacklisting: " + str(len(dd_key_set)))
print("\t[INFO] Number of AD digests before blacklisting: " + str(len(ad_key_set)))
print("\t[INFO] Number of UD digests before blacklisting: " + str(len(ud_key_set)))
print("\t[INFO] Number of CD digests before blacklisting: " + str(len(cd_key_set)))

# Subtract blacklisted digests
dd_key_set_bl = dd_key_set - digest_to_be_excluded
ddad_key_set_bl = ddad_key_set - digest_to_be_excluded
ad_key_set_bl = ad_key_set - digest_to_be_excluded
ud_key_set_bl = ud_key_set - digest_to_be_excluded
cd_key_set_bl = cd_key_set - digest_to_be_excluded

print("\t[INFO] ------------------------------------------------------- " )
print("\t[INFO] Number of DD digests after blacklisting: " + str(len(dd_key_set_bl)) + " (" + "{:.2%}".format(len(dd_key_set_bl)/len(dd_key_set)) + ")")
print("\t[INFO] Number of UD digests after blacklisting: " + str(len(ud_key_set_bl)) + " (" + "{:.2%}".format(len(ud_key_set_bl)/len(ud_key_set)) + ")")
print("\t[INFO] Number of AD digests after blacklisting: " + str(len(ad_key_set_bl)) + " (" + "{:.2%}".format(len(ad_key_set_bl)/len(ad_key_set)) + ")")
print("\t[INFO] Number of CD digests after blacklisting: " + str(len(cd_key_set_bl)) + " (" + "{:.2%}".format(len(cd_key_set_bl)/(len(dd_key_set) + len(ud_key_set) + len(ad_key_set))) + ")")


### Determine digest densities and write to BedGraph files
##########################################################

# Define UCSC track options
track_color = "255,161,0"
track_altColor ="0,100,200"
track_maxHeightPixels = 60

track_name = data_set_track_tag + " - DENSITY - DDAD"
track_description = track_name
bedgraph_stream_ddad_digests_density_output = open(out_prefix + "_ddad_digests_density.bedgraph", 'wt')
bedgraph_stream_ddad_digests_density_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_name = data_set_track_tag + " - DENSITY - UD"
track_description = track_name
bedgraph_stream_ud_digests_density_output = open(out_prefix + "_ud_digests_density.bedgraph", 'wt')
bedgraph_stream_ud_digests_density_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_name = data_set_track_tag + " - DENSITY - CD"
track_description = track_name
bedgraph_stream_cd_digests_density_output = open(out_prefix + "_cd_digests_density.bedgraph", 'wt')
bedgraph_stream_cd_digests_density_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_name = data_set_track_tag + " - DENSITY - ND"
track_description = track_name
bedgraph_stream_nd_digests_density_output = open(out_prefix + "_nd_digests_density.bedgraph", 'wt')
bedgraph_stream_nd_digests_density_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_name = data_set_track_tag + " - DENSITY - BL"
track_description = track_name
bedgraph_stream_bl_digests_density_output = open(out_prefix + "_bl_digests_density.bedgraph", 'wt')
bedgraph_stream_bl_digests_density_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

print("[INFO] Iterating digest file to create BedGraph tracks for interaction distances ...")

# Init variables for sliding digest window

digest_queue = []
last_chr = "chr1"
current_window_size = 0

ddad_cnt = 0
ddad_density = 0.0
ud_cnt = 0
ud_density = 0.0
cd_cnt = 0
cd_density = 0.0
nd_cnt = 0
nd_density = 0.0
bl_cnt = 0
bl_density = 0.0
chr_d_cnt = 0

with open(gopher_digest_file, 'rt') as fp:

    line = fp.readline()
    line = fp.readline() # Skip first line

    n_progress = 0
    while line:

        # Report progress
        n_progress += 1
        if n_progress % 100000 == 0:
            print("\tProcessed " + str(n_progress) + " digests ...")

        # Get key for current digest region in GOPHER's digest file
        fields = line.split('\t')
        coord_key = fields[0] + '\t' + fields[1] + '\t' + fields[2]

        # Set everything to zero when encountering a new chromosome
        if last_chr != fields[0]:

            # Empty queue
            while 0<len(digest_queue):
                digest_queue.pop(0)

            # Reset variables for sliding digest window
            ddad_cnt = 0
            ddad_density = 0.0
            ud_cnt = 0
            ud_density = 0.0
            cd_cnt = 0
            cd_density = 0.0
            bl_cnt = 0
            bl_density = 0.0
            nd_cnt = 0
            nd_density = 0.0
            current_window_size = 0
            last_chr = fields[0]
            chr_d_cnt = 0

        # Fill queue
        digest_queue.append(coord_key)
        chr_d_cnt +=1

        # Update current window size
        current_window_size = current_window_size + (int(fields[2])-int(fields[1]))

        # Directed
        if coord_key in ddad_key_set_bl:
            ddad_cnt += 1

        # Undirected
        if coord_key in ud_key_set_bl:
            ud_cnt += 1

        # Directed or undirected
        if coord_key in cd_key_set_bl:
            cd_cnt +=1

        # Not interacting
        if coord_key not in cd_key_set:
            nd_cnt = nd_cnt + 1

        # Blacklisted
        if coord_key in digest_to_be_excluded:
            bl_cnt += 1

        # Update digest densities
        if len(digest_queue) == digest_density_range + 1:

            # Get old digest from queue
            old_digest_from_q = digest_queue.pop(0)

            # Update current window size
            current_window_size = current_window_size - (int(old_digest_from_q.split('\t')[2]) - int(old_digest_from_q.split('\t')[1]))

            # Directed
            if old_digest_from_q in ddad_key_set_bl:
                ddad_cnt = ddad_cnt - 1
            ddad_density = ddad_cnt/current_window_size

            # Undirected
            if old_digest_from_q in ud_key_set_bl:
                ud_cnt = ud_cnt - 1
            ud_density = ud_cnt/current_window_size

            # Blacklisted
            if old_digest_from_q in digest_to_be_excluded:
                bl_cnt = bl_cnt - 1
            bl_density = bl_cnt/current_window_size

            # Directed or undirected
            if old_digest_from_q in cd_key_set_bl:
                cd_cnt = cd_cnt - 1
            cd_density = cd_cnt/current_window_size

            # Not interacting
            if old_digest_from_q not in cd_key_set:
                nd_cnt = nd_cnt - 1
            nd_density = nd_cnt/current_window_size

            # Keep track of current chromosome
            last_chr = fields[0]

        # Center window
        window_center_pos = int(digest_density_range/2) - 1 + digest_density_range%2 # works for uneven and even window sizes

        # Write window centered densities to file
        if chr_d_cnt - 1 < window_center_pos:
            bedgraph_stream_ddad_digests_density_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(0.0) + '\n')
            bedgraph_stream_ud_digests_density_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(0.0) + '\n')
            bedgraph_stream_cd_digests_density_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(0.0) + '\n')
            bedgraph_stream_bl_digests_density_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(0.0) + '\n')
            bedgraph_stream_nd_digests_density_output.write(fields[0] + '\t' + str(int(fields[1]) - 1) + '\t' + str(int(fields[2])) + '\t' + str(0.0) + '\n')
        elif digest_density_range <= chr_d_cnt:
            bedgraph_stream_ddad_digests_density_output.write(str(digest_queue[window_center_pos]) + '\t' + str(ddad_density) + '\n')
            bedgraph_stream_ud_digests_density_output.write(str(digest_queue[window_center_pos]) + '\t' + str(ud_density) + '\n')
            bedgraph_stream_cd_digests_density_output.write(str(digest_queue[window_center_pos]) + '\t' + str(cd_density) + '\n')
            bedgraph_stream_bl_digests_density_output.write(str(digest_queue[window_center_pos]) + '\t' + str(bl_density) + '\n')
            bedgraph_stream_nd_digests_density_output.write(str(digest_queue[window_center_pos]) + '\t' + str(nd_density) + '\n')

        line = fp.readline()

bedgraph_stream_ddad_digests_density_output.close()
bedgraph_stream_ud_digests_density_output.close()
bedgraph_stream_cd_digests_density_output.close()
bedgraph_stream_bl_digests_density_output.close()
bedgraph_stream_nd_digests_density_output.close()

print("[INFO] ... done.")
