"""
This script is to prepare BED files containing interacting digests and promoter regions used for motif analysis.
It iterates a file with interactions, gene symbols,
strand and TSS information, etc. created with the script 'get_gene_symbols_of_interactions_script.py'.

Six BED files are created:

   1. BED files with regions of unique exclusive digests (see schematic representation with venn diagrams)

      a. of directed interactions ('S' or 'T' in column 3 of file with gene symbols and interactions)
         '<OUT_PREFIX>_directed_digests.bed'

      b. of undirected reference interactions ('URAA' in column 3 of file with gene symbols and interactions)
         '<OUT_PREFIX>_undirected_digests.bed'

   2. BED files for promoter regions

      a. on directed interacting digests
         i. and on the '-' strand
            '<OUT_PREFIX>_directed_minus_tss.bed'

         ii. and on the '+' strand
            '<OUT_PREFIX>_directed_plus_tss.bed'

      b. on undirected interacting digests
         i. and on the '-' strand
            '<OUT_PREFIX>_undirected_minus_tss.bed'

         ii. and on the '+' strand
            '<OUT_PREFIX>_undirected_plus_tss.bed'

Use copy these lines to a text file:

chr9:128781853-128788644;chr9:128855420-128861085	66776	S	TBC1D13;KYAT1	12:0	AA	8.32	+/-	chr9:128787253:+,chr9:128787207:+;chr9:128860532:-
chr17:39397580-39405334;chr17:39461300-39463357	55966	T	FBXL20;CDK12	33:11	AA	7.37	-/+	chr17:39402556:-,chr17:39401626:-;chr17:39461486:+
chr3:101509230-101513979;chr3:101842446-101853540	328467	URAA	SENP7;NFKBIZ	6:3	AA	1.37	-/+	chr3:101513212:-,chr3:101513241:-;chr3:101849514:+
chr21:30736306-30747469;chr21:30819004-30819702	71535	URAI	KRTAP21-2;	8:4	AI	1.64	-/-1	chr21:30747259:-;
chr5:74778161-74782279;chr5:74867764-74868683	85485	URII	;	17:6	II	4.05	-1/-1	;
chr5:88676106-88682018;chr5:88927332-88929569	245314	NA	LINC00461,MEF2C-AS2;	5:1	II	2.21	d/-1	chr5:88676298:-,chr5:88678448:-,chr5:88676218:+;

gzip the test file and use this file as input for this script.
Digest and promoter will be extracted only for the first two lines (directed interactions) and the third line (undirected reference interaction).

The BED files with digest regions are used for motif discovery using DREME.
The BED files with promoter regions are used the analysis of motif occurrences relative to TSS using CentriMo.

We use 'bedtools getfasta' for sequence extraction and gsub for repeat masking.

We use 'gsub(/a|c|g|t/,"N",$1)' for repeat masking.

"""

import argparse
import gzip
import os
import numpy


### Parse command line
######################

parser = argparse.ArgumentParser(description='Extract regions for motif analysis from file with interactions and gene symbols.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--interaction-gs-file', help='Interaction file with gene symbols created with \'get_gene_symbols_of_interactions_script.py\'.')
parser.add_argument('--chrom-info-file', help='File with chromosome sizes <CHROMOSOME NAME>\\t<CHROMOSOME SIZE>.')
parser.add_argument('--up-dist', help='Number of bases upstream of TSS.', default=1000)
parser.add_argument('--down-dist', help='Number of bases upstream of TSS.', default=1000)
parser.add_argument('--allowed-enrichment-pair-tags', help='Set of allowed enrichment pair tags for digests (\'AA\', \'AI\',\'IA\',\'II\', ...) separated by \';\'.', default="AA")
parser.add_argument('--allowed-strand-pair-tags', help='Set of allowed strand pair tags for digests (\'-/-\', \'+/+\', \'+/-\',\'-/+\',\'-/d\', ...) separated by \';\'.', default="-/-;+/+")
parser.add_argument('--allowed-interaction-categories-directed', help='Set categories of directed interactions (\'S\' and \'T\') separated by \';\'.', default="S;T")
parser.add_argument('--allowed-interaction-categories-undirected', help='Set categories of directed interactions (\'NA\', \'U\', \'URAA\', \'URAI\' and \'URII\') separated by \';\'.', default="URAA")
parser.add_argument('--max-num-of-tss-per-digest', help='Use interactions for which both digests have at most this number of TSS only.', default=-1)
parser.add_argument('--bedtools-path', help='Path to bedtools executable.', default="bedtools")
parser.add_argument('--genome-fasta-path', help='Path to genome FASTA file.', default="/Users/hansep/data/hg38/hg38.fa")

args = parser.parse_args()
out_prefix = args.out_prefix
interaction_gs_file = args.interaction_gs_file # central file that contains a lot of information about interactions
chrom_info_file = args.chrom_info_file # needed to avoid promoter coordinates outside chromosome regions
up_dist = int(args.up_dist) # used to determine start coordinates of promoter regions
down_dist = int(args.down_dist) # used to determine end coordinates of promoter regions
allowed_enrichment_pair_tags = args.allowed_enrichment_pair_tags # only interactions with these enrichment pair tags will be used
allowed_strand_pair_tags = args.allowed_strand_pair_tags # only interactions with these strand pair tags will be used
allowed_interaction_categories_directed = args.allowed_interaction_categories_directed # only directed interactions with these categories will be used
allowed_interaction_categories_undirected = args.allowed_interaction_categories_undirected # only undirected with these categories will be used
bedtools_path = args.bedtools_path
genome_fasta_path = args.genome_fasta_path

print("[INFO] " + "Input parameters")
print("\t[INFO] Analysis for: " + out_prefix)
print("\t[INFO] --interaction-gs-file: " + interaction_gs_file)
print("\t[INFO] --chrom-info-file: " + chrom_info_file)
print("\t[INFO] --up-dist: " + str(up_dist))
print("\t[INFO] --down-dist: " + str(down_dist))
print("\t[INFO] --enrichment-strand-pair-tags: " + allowed_enrichment_pair_tags)
print("\t[INFO] --allowed-strand-pair-tags: " + allowed_strand_pair_tags)
print("\t[INFO] --allowed-interaction-categories-directed: " + allowed_interaction_categories_directed)
print("\t[INFO] --allowed-interaction-categories-undirected: " + allowed_interaction_categories_undirected)
print("\t[INFO] --bedtools-path: " + bedtools_path)
print("\t[INFO] --genome-fasta-path: " + genome_fasta_path)

# convert allowed pair tag strings into lists
allowed_strand_pair_tags = allowed_strand_pair_tags.split(";")
#print(allowed_strand_pair_tags)
allowed_enrichment_pair_tags = str(allowed_enrichment_pair_tags).split(";")
#print(allowed_enrichment_pair_tags)
allowed_interaction_categories_directed = str(allowed_interaction_categories_directed).split(";")
#print(allowed_interaction_categories_directed)
allowed_interaction_categories_undirected = str(allowed_interaction_categories_undirected).split(";")
#print(allowed_interaction_categories_undirected)

### Define auxiliary functions
##############################

def convert_bed_to_fasta(bedtools_path, genome_fasta_path, bed_file_name):
    bed_file_name_base = str(bed_file_name).split(".b")[0]
    sys_cmd = bedtools_path + ' getfasta -name -fi ' + genome_fasta_path + ' -bed  ' + bed_file_name + ' > ' + bed_file_name_base + '.fasta'
    #print("Will execute: " + sys_cmd)
    os.system(sys_cmd)
    return bed_file_name_base + '.fasta'

def mask_repeats(fasta_file_name):
    fasta_file_name_base = str(fasta_file_name).split(".f")[0]
    sys_cmd = 'awk \'{if($1 !~ /^>/){gsub(/a|c|g|t/,"N",$1)};print}\' ' + fasta_file_name + ' > ' + fasta_file_name_base + '_masked.fasta'
    #print("Will execute: " + sys_cmd)
    os.system(sys_cmd)
    return fasta_file_name_base + '_masked.fasta'

def get_base_frequencies(fasta_file_name):
    a_abs = 0
    c_abs = 0
    g_abs = 0
    t_abs = 0
    A_abs = 0
    C_abs = 0
    G_abs = 0
    T_abs = 0
    N_abs = 0
    with open(fasta_file_name, 'rt') as fp:
        line = fp.readline()
        while line:
            if line.count('>') == 0:
                line = line.rstrip()
                a_abs = a_abs + line.count('a')
                c_abs = c_abs + line.count('c')
                g_abs = g_abs + line.count('g')
                t_abs = t_abs + line.count('t')
                A_abs = A_abs + line.count('A')
                C_abs = C_abs + line.count('C')
                G_abs = G_abs + line.count('G')
                T_abs = T_abs + line.count('T')
                N_abs = N_abs + line.count('N')
            line = fp.readline()
    fp.close()

    # calculate repeat content (fraction of lower case)
    acgt_abs = a_abs + c_abs + g_abs + t_abs
    acgtACGTN_abs = a_abs + c_abs + g_abs + t_abs + A_abs + C_abs + G_abs + T_abs + N_abs
    repeat_content = "{:.2f}".format(100*acgt_abs/acgtACGTN_abs)

    # determine GC content within repeat regions
    cg_abs = c_abs + g_abs
    if 0 < acgt_abs:
        gc_content_repeat = "{:.2f}".format(100 * cg_abs / acgt_abs)
    else:
        gc_content_repeat = "{:.2f}".format(0.0)

    # determine GC content for non repeat regions
    CG_abs = C_abs + G_abs
    ACGT_abs = A_abs + C_abs + G_abs + T_abs
    gc_content_non_repeat = "{:.2f}".format(100 * CG_abs / ACGT_abs)

    # determine overall GC content
    cgCG_abs = cg_abs + CG_abs
    acgtACGT_abs = a_abs + c_abs + g_abs + t_abs + A_abs + C_abs + G_abs + T_abs
    gc_content_total = "{:.2f}".format(100 * cgCG_abs / acgtACGT_abs)

    header_line = "fasta_file_name" + "\t" + \
                  "a_occ_abs" + "\t" + \
                  "c_occ_abs" + "\t" + \
                  "g_occ_abs" + "\t" + \
                  "t_occ_abs" + "\t" + \
                  "A_occ_abs" + "\t" + \
                  "C_occ_abs" + "\t" + \
                  "G_occ_abs" + "\t" + \
                  "T_occ_abs" + "\t" + \
                  "N_occ_abs" + "\t" + \
                  "repeat_content" + "\t" + \
                  "gc_content_repeat" + "\t" + \
                  "gc_content_non_repeat" + "\t" + \
                  "gc_content_total"

    value_line = fasta_file_name + "\t" + \
                 str(a_abs) + "\t" + \
                 str(c_abs) + "\t" + \
                 str(g_abs) + "\t" + \
                 str(t_abs) + "\t" + \
                 str(A_abs) + "\t" + \
                 str(C_abs) + "\t" + \
                 str(G_abs) + "\t" + \
                 str(T_abs) + "\t" + \
                 str(N_abs) + "\t" + \
                 str(repeat_content) + "\t" + \
                 str(gc_content_repeat) + "\t" + \
                 str(gc_content_non_repeat) + "\t" + \
                 str(gc_content_total)

    return header_line, value_line, repeat_content, gc_content_repeat, gc_content_non_repeat, gc_content_total


### Prepare output files
########################

# create two BED files for regions of segmented directed and undirected interactions
bed_file_name_directed_digests = out_prefix + "_directed_digests.bed"
bed_stream_name_directed_digests = open(bed_file_name_directed_digests, 'wt')
bed_file_name_undirected_digests = out_prefix + "_undirected_digests.bed"
bed_stream_name_undirected_digests = open(bed_file_name_undirected_digests, 'wt')

# create two BED files for promoters on digests of directed iteractions, one for TSS on the plus and another one for TSS on the minus strand
bed_file_name_directed_digests_tss_minus = out_prefix + "_directed_digests_tss_minus.bed"
bed_stream_name_directed_digests_tss_minus = open(bed_file_name_directed_digests_tss_minus, 'wt')
bed_file_name_directed_digests_tss_plus = out_prefix + "_directed_digests_tss_plus.bed"
bed_stream_name_directed_digests_tss_plus = open(bed_file_name_directed_digests_tss_plus, 'wt')

# create two BED files for promoters on digests of undirected iteractions, one for TSS on the plus and another one for TSS on the minus strand
bed_file_name_undirected_digests_tss_minus = out_prefix + "_undirected_digests_tss_minus.bed"
bed_stream_name_undirected_digests_tss_minus = open(bed_file_name_undirected_digests_tss_minus, 'wt')
bed_file_name_undirected_digests_tss_plus = out_prefix + "_undirected_digests_tss_plus.bed"
bed_stream_name_undirected_digests_tss_plus = open(bed_file_name_undirected_digests_tss_plus, 'wt')

### Prepare variables and data structures
#########################################

chrom_sizes_dict = {} # read chromosome sizes to hash
with open(chrom_info_file, 'rt') as fp:
    line = fp.readline()
    while line:
        chr_name = line.split("\t")[0]
        chr_size = line.split("\t")[1]
        chrom_sizes_dict[chr_name] = int(chr_size)
        line = fp.readline()

directed_digests_set = set()
undirected_digests_set = set()

directed_interaction_num = 0
undirected_interaction_num = 0

first_digest_without_tss_num = 0
second_digest_without_tss_num = 0

directed_digests_symbols_set = set()
undirected_digests_symbols_set = set()

interaction_partners_per_digest_directed_dict = {}
interaction_partners_per_digest_undirected_dict = {}


### Iterate interaction file with gene symbols
##############################################

print("[INFO] Iterating interaction file with gene symbols " + interaction_gs_file + " ...")
with gzip.open(interaction_gs_file, 'rt') as fp:

    n_interaction_total = 0 # counter to track progress
    line = fp.readline()
    cnt = 0
    while line:

        line = line.rstrip()

        n_interaction_total += 1
        if n_interaction_total % 1000000 == 0:
            print("\t[INFO]", n_interaction_total, "interactions processed ...")

        field = line.split("\t")

        if field[5] not in allowed_enrichment_pair_tags: # use only interactions with specified digest enrichment pair tags (II, AI, AA)
            line = fp.readline()
            continue

        if field[7] not in allowed_strand_pair_tags and 'All' not in allowed_strand_pair_tags: # use only interactions with specified digest strand pair tags (-/-, +/+, -/+, +/-, -/d, ...) or 'All'
            line = fp.readline()
            continue

        if field[8].split(";")[0] == '':
            print("Warning: No TSS for first digest!")
            first_digest_without_tss_num += 1
            line = fp.readline()
            continue

        if field[8].split(";")[1] == '':
            print("Warning: No TSS for second digest!")
            second_digest_without_tss_num += 1
            line = fp.readline()
            continue

        # parse interaction info ToDo: write parsing function (saves 20 lines)
        coords = field[0].split(";")
        coord_a = coords[0]
        coord_b = coords[1]
        chr_a = coord_a.split(":")[0]
        chr_b = coord_b.split(":")[0]
        sta_a = int(coord_a.split(":")[1].split("-")[0])
        sta_b = int(coord_b.split(":")[1].split("-")[0])
        end_a = int(coord_a.split(":")[1].split("-")[1])
        end_b = int(coord_b.split(":")[1].split("-")[1])
        syms = field[3].split(";")
        syms_a = syms[0]
        syms_b = syms[1]

        base_name = field[0] + "|" + field[1] + "|" + field[2] + "|" + field[3] + "|" + field[4]+ "|" + field[5] + "|" + field[6] + "|" + field[7]

        # get category of interacting digest
        interaction_category = field[2]

        if interaction_category in allowed_interaction_categories_directed:

            directed_interaction_num += 1

            # add digest regions to set for directed interactions
            directed_digests_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
            directed_digests_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))

            # add associated gene symbols to set for directed interactions
            for symbol in syms_a.split(","):
                directed_digests_symbols_set.add(symbol)
            for symbol in syms_b.split(","):
                directed_digests_symbols_set.add(symbol)

            # count interaction partners for each directed digest
            if chr_a + "\t" + str(sta_a) + "\t" + str(end_a) in interaction_partners_per_digest_directed_dict.keys():
                interaction_partners_per_digest_directed_dict[chr_a + "\t" + str(sta_a) + "\t" + str(end_a)] += 1
            else:
                interaction_partners_per_digest_directed_dict[chr_a + "\t" + str(sta_a) + "\t" + str(end_a)] = 1

            if chr_b + "\t" + str(sta_b) + "\t" + str(end_b) in interaction_partners_per_digest_directed_dict.keys():
                interaction_partners_per_digest_directed_dict[chr_b + "\t" + str(sta_b) + "\t" + str(end_b)] += 1
            else:
                interaction_partners_per_digest_directed_dict[chr_b + "\t" + str(sta_b) + "\t" + str(end_b)] = 1

        if interaction_category in allowed_interaction_categories_undirected:

            undirected_interaction_num += 1

            # add digest regions to set for undirected interactions
            undirected_digests_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
            undirected_digests_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))

            # add associated gene symbols to set for undirected interactions
            for symbol in syms_a.split(","):
                undirected_digests_symbols_set.add(symbol)
            for symbol in syms_b.split(","):
                undirected_digests_symbols_set.add(symbol)

            # count interaction partners for each undirected digest
            if chr_a + "\t" + str(sta_a) + "\t" + str(end_a) in interaction_partners_per_digest_undirected_dict.keys():
                interaction_partners_per_digest_undirected_dict[chr_a + "\t" + str(sta_a) + "\t" + str(end_a)] += 1
            else:
                interaction_partners_per_digest_undirected_dict[chr_a + "\t" + str(sta_a) + "\t" + str(end_a)] = 1

            if chr_b + "\t" + str(sta_b) + "\t" + str(end_b) in interaction_partners_per_digest_undirected_dict.keys():
                interaction_partners_per_digest_undirected_dict[chr_b + "\t" + str(sta_b) + "\t" + str(end_b)] += 1
            else:
                interaction_partners_per_digest_undirected_dict[chr_b + "\t" + str(sta_b) + "\t" + str(end_b)] = 1

        line = fp.readline()

# report interactions with no TSS on one or both digests (should not happen)
if 0 < first_digest_without_tss_num:
    print("\t[WARNING] There were " + str(first_digest_without_tss_num) + " interactions without TSS on the first digest!")
if 0 < second_digest_without_tss_num:
    print("\t[WARNING] There were " + str(second_digest_without_tss_num) + " interactions without TSS on the second digest!")

fp.close()

print("[INFO] ... done.")

### Segment digests and print coordinates for sequence analysis
###############################################################

print("[INFO] Segmenting digests and printing coordinates for sequence analysis ...")

directed_wo_undirected_digests_set = directed_digests_set.difference(undirected_digests_set)
undirected_wo_directed_digests_set = undirected_digests_set.difference(directed_digests_set)
directed_intersect_undirected_digests_set = directed_digests_set.intersection(undirected_digests_set)
directed_union_directed_digests_set = directed_digests_set.union(undirected_digests_set)

cnt = 1
print("\t[INFO] Printing to file: " + bed_file_name_directed_digests)
for digest in directed_wo_undirected_digests_set:
    bed_stream_name_directed_digests.write(digest + "\t" + str(cnt) + "\n")
    cnt += 1
cnt = 1
print("\t[INFO] Printing to file: " + bed_file_name_undirected_digests)
for digest in undirected_wo_directed_digests_set:
    bed_stream_name_undirected_digests.write(digest + "\t" + str(cnt) + "\n")
    cnt += 1

bed_stream_name_directed_digests.close()
bed_stream_name_undirected_digests.close()

print("[INFO] ... done.")

### Segement and print gene symbols for GO analysis
###################################################

print("[INFO] Segmenting and printing genes symbols for GO analysis ...")

directed_wo_undirected_digests_symbols_set = directed_digests_symbols_set.difference(undirected_digests_symbols_set)
tab_file_name_directed_digest_symbols = out_prefix + "_directed_digest_symbols.tab"
tab_stream_name_directed_digest_symbols = open(tab_file_name_directed_digest_symbols, 'wt')
print("\t[INFO] Printing to file: " + tab_file_name_directed_digest_symbols)
for symbol in directed_wo_undirected_digests_symbols_set:
    tab_stream_name_directed_digest_symbols.write(symbol + "\n")
tab_stream_name_directed_digest_symbols.close()

undirected_wo_directed_digests_symbols_set = undirected_digests_symbols_set.difference(directed_digests_symbols_set)
tab_file_name_undirected_digest_symbols = out_prefix + "_undirected_digest_symbols.tab"
tab_stream_name_undirected_digest_symbols = open(tab_file_name_undirected_digest_symbols, 'wt')
print("\t[INFO] Printing to file: " + tab_file_name_undirected_digest_symbols)
for symbol in undirected_wo_directed_digests_symbols_set:
    tab_stream_name_undirected_digest_symbols.write(symbol + "\n")
tab_stream_name_undirected_digest_symbols.close()

print("[INFO] ... done.")

### Determine and print digest sizes
####################################

print("[INFO] Determining and printing digest sizes ...")

print("\t[INFO] Directed digests ...")
tab_file_name_directed_digests_size_distribution = out_prefix + "_directed_digests_size_distribution.tab"
tab_stream_name_directed_digests_size_distribution = open(tab_file_name_directed_digests_size_distribution, 'wt')
print("\t\t[INFO] Printing to file: " + tab_file_name_directed_digests_size_distribution)
directed_digest_size_array = []
for d_coord in directed_wo_undirected_digests_set:
    sta = d_coord.split("\t")[1]
    end = d_coord.split("\t")[2]
    size = int(end) - int(sta)
    directed_digest_size_array.append(size)
    tab_stream_name_directed_digests_size_distribution.write(str(size) + "\n")

print("\t\t[INFO] Q1: " + str(int(numpy.quantile(directed_digest_size_array, 0.25))) + " bp")
print("\t\t[INFO] Q2: " + str(int(numpy.quantile(directed_digest_size_array, 0.5))) + " bp (Median)")
print("\t\t[INFO] Q3: " + str(int(numpy.quantile(directed_digest_size_array, 0.75))) + " bp")
print("\t\t[INFO] Mean: " + str(int(numpy.mean(directed_digest_size_array))) + " bp")
tab_stream_name_directed_digests_size_distribution.close()

print("\t[INFO] Undirected digests ...")
tab_file_name_undirected_digests_size_distribution = out_prefix + "_undirected_digests_size_distribution.tab"
tab_stream_name_undirected_digests_size_distribution = open(tab_file_name_undirected_digests_size_distribution, 'wt')
print("\t\t[INFO] Printing to file: " + tab_file_name_undirected_digests_size_distribution)
undirected_digest_size_array = []
for d_coord in undirected_wo_directed_digests_set:
    sta = d_coord.split("\t")[1]
    end = d_coord.split("\t")[2]
    size = int(end) - int(sta)
    undirected_digest_size_array.append(size)
    tab_stream_name_undirected_digests_size_distribution.write(str(size) + "\n")

print("\t\t[INFO] Q1: " + str(int(numpy.quantile(undirected_digest_size_array, 0.25))) + " bp")
print("\t\t[INFO] Q2: " + str(int(numpy.quantile(undirected_digest_size_array, 0.5))) + " bp (Median)")
print("\t\t[INFO] Q3: " + str(int(numpy.quantile(undirected_digest_size_array, 0.75))) + " bp")
print("\t\t[INFO] Mean: " + str(int(numpy.mean(undirected_digest_size_array))) + " bp")
tab_stream_name_undirected_digests_size_distribution.close()

print("[INFO] ... done.")

### Create file with numbers of interaction partners per digest
###############################################################

print("[INFO] Creating file with numbers of interaction partners per digest ...")

print("\t[INFO] Directed digests ...")
tab_file_name_directed_digests_interaction_partner_per_digest_distribution = out_prefix + "_directed_digests_interaction_partner_per_digest_distribution.tab"
tab_stream_name_directed_digests_interaction_partner_per_digest_distribution = open(tab_file_name_directed_digests_interaction_partner_per_digest_distribution, 'wt')
print("\t\t[INFO] Printing to file: " + tab_file_name_directed_digests_interaction_partner_per_digest_distribution)
directed_inter_partner_num_array = []
for inter_partner_num in interaction_partners_per_digest_directed_dict.values():
    tab_stream_name_directed_digests_interaction_partner_per_digest_distribution.write(str(inter_partner_num) + "\n")
    directed_inter_partner_num_array.append(inter_partner_num)
print("\t\t[INFO] Mean number of partners: " + "{:.2f}".format((numpy.mean(directed_inter_partner_num_array))))
tab_stream_name_directed_digests_interaction_partner_per_digest_distribution.close()

print("\t[INFO] Undirected digests ...")
tab_file_name_undirected_digests_interaction_partner_per_digest_distribution = out_prefix + "_undirected_digests_interaction_partner_per_digest_distribution.tab"
tab_stream_name_undirected_digests_interaction_partner_per_digest_distribution = open(tab_file_name_undirected_digests_interaction_partner_per_digest_distribution, 'wt')
print("\t\t[INFO] Printing to file: " + tab_file_name_undirected_digests_interaction_partner_per_digest_distribution)
undirected_inter_partner_num_array = []
for inter_partner_num in interaction_partners_per_digest_undirected_dict.values():
    tab_stream_name_undirected_digests_interaction_partner_per_digest_distribution.write(str(inter_partner_num) + "\n")
    undirected_inter_partner_num_array.append(inter_partner_num)
print("\t\t[INFO] Mean number of partners: " + "{:.2f}".format((numpy.mean(undirected_inter_partner_num_array))))
tab_stream_name_undirected_digests_interaction_partner_per_digest_distribution.close()

print("[INFO] ... done.")

#######################################################################################################################################



### get TSS associated with digests
###################################

directed_wo_undirected_digests_tss_set = set()
undirected_wo_directed_digests_tss_set = set()

cnt=0
print("[INFO] Iterating interaction file with gene symbols " + interaction_gs_file + " ...")
with gzip.open(interaction_gs_file, 'rt') as fp:

    n_interaction_total = 0 # counter to track progress
    line = fp.readline()
    cnt = 0
    while line:
        line = line.rstrip()

        n_interaction_total += 1
        if n_interaction_total % 1000000 == 0:
            print("\t[INFO]", n_interaction_total, "interactions processed ...")

        field = line.split("\t")

        if field[5] not in allowed_enrichment_pair_tags: # use only interactions with specified digest enrichment pair tags (II, AI, AA)
            line = fp.readline()
            continue

        if field[7] not in allowed_strand_pair_tags and 'All' not in allowed_strand_pair_tags: # use only interactions with specified digest strand pair tags (-/-, +/+, -/+, +/-, -/d, ...) or 'All'
            line = fp.readline()
            continue

        if field[8].split(";")[0] == '':
            print("Warning: No TSS for first digest!")
            first_digest_without_tss_num += 1
            line = fp.readline()
            continue

        if field[8].split(";")[1] == '':
            print("Warning: No TSS for second digest!")
            second_digest_without_tss_num += 1
            line = fp.readline()
            continue

        # get coordinates of interacting digests
        coords = field[0].split(";")
        coord_a = coords[0]
        coord_b = coords[1]
        chr_a = coord_a.split(":")[0]
        chr_b = coord_b.split(":")[0]
        sta_a = int(coord_a.split(":")[1].split("-")[0])
        sta_b = int(coord_b.split(":")[1].split("-")[0])
        end_a = int(coord_a.split(":")[1].split("-")[1])
        end_b = int(coord_b.split(":")[1].split("-")[1])
        syms = field[3].split(";")
        syms_a = syms[0]
        syms_b = syms[1]
        tsss = field[8].split(";")
        tsss_a = tsss[0]
        tsss_b = tsss[1]

        base_name = field[0] + "|" + field[1] + "|" + field[2] + "|" + field[3] + "|" + field[4]+ "|" + field[5] + "|" + field[6] + "|" + field[7]

        d1_coord = str(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
        if d1_coord in directed_wo_undirected_digests_set:
            for tss_a in tsss_a.split(','):
                directed_wo_undirected_digests_tss_set.add(tss_a)
        if d1_coord in undirected_wo_directed_digests_set:
            for tss_a in tsss_a.split(','):
                undirected_wo_directed_digests_tss_set.add(tss_a)

        d2_coord = str(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))
        if d2_coord in directed_wo_undirected_digests_set:
            for tss_b in tsss_b.split(','):
                directed_wo_undirected_digests_tss_set.add(tss_b)
        if d2_coord in undirected_wo_directed_digests_set:
            for tss_b in tsss_b.split(','):
                undirected_wo_directed_digests_tss_set.add(tss_b)

        line = fp.readline()

cnt_minus = 1
cnt_plus = 1
for tss in directed_wo_undirected_digests_tss_set:
    A = tss.split(':')
    chr = A[0]
    sta = int(A[1]) - up_dist
    if sta < 0:
        sta = 0
    end = int(A[1]) + down_dist
    if chrom_sizes_dict[chr] < end:
        end = chrom_sizes_dict[chr] - 1
    strand = A[2]
    if strand == '-':
        bed_stream_name_directed_digests_tss_minus.write(chr + "\t" + str(sta) + "\t" + str(end) + "\t" + str(cnt_minus) + "\n")
        cnt_minus += 1
    elif strand == '+':
        bed_stream_name_directed_digests_tss_plus.write(chr + "\t" + str(sta) + "\t" + str(end) + "\t" + str(cnt_plus) + "\n")
        cnt_plus += 1
    else:
        print("Warning: Strand of TSS was neither \'-\' nor \'+\'!")

cnt_minus = 1
cnt_plus = 1
for tss in undirected_wo_directed_digests_tss_set:
    A = tss.split(':')
    chr = A[0]
    sta = int(A[1]) - up_dist
    if sta < 0:
        sta = 0
    end = int(A[1]) + down_dist
    if chrom_sizes_dict[chr] < end:
        end = chrom_sizes_dict[chr] - 1
    strand = A[2]
    if strand == '-':
        bed_stream_name_undirected_digests_tss_minus.write(chr + "\t" + str(sta) + "\t" + str(end) + "\t" + str(cnt_minus) + "\n")
        cnt_minus += 1
    elif strand == '+':
        bed_stream_name_undirected_digests_tss_plus.write(chr + "\t" + str(sta) + "\t" + str(end) + "\t" + str(cnt_plus) + "\n")
        cnt_plus += 1
    else:
        print("Warning: Strand of TSS was neither \'-\' nor \'+\'!")

bed_stream_name_directed_digests_tss_minus.close()
bed_stream_name_directed_digests_tss_plus.close()
bed_stream_name_undirected_digests_tss_minus.close()
bed_stream_name_undirected_digests_tss_plus.close()


### Create FASTA files with digest and promoter sequences
#########################################################

print("[INFO] Converting BED to FASTA ...")

# get FASTA files with DIGEST sequences of directed interactions
fasta_file_name_directed_digests = convert_bed_to_fasta(bedtools_path, genome_fasta_path, bed_file_name_directed_digests)
fasta_file_name_directed_digests_masked = mask_repeats(fasta_file_name_directed_digests)

# get FASTA files with DIGEST sequences of undirected interactions
fasta_file_name_undirected_digests = convert_bed_to_fasta(bedtools_path, genome_fasta_path, bed_file_name_undirected_digests)
fasta_file_name_undirected_digests_masked = mask_repeats(fasta_file_name_undirected_digests)

# get FASTA files with sequences around TSS on '-' strand that are associated with directed interactions for CentriMo analysis
fasta_file_name_directed_digests_tss_minus = convert_bed_to_fasta(bedtools_path, genome_fasta_path, bed_file_name_directed_digests_tss_minus)

# get FASTA files with sequences around TSS on '+' strand that are associated with directed interactions for CentriMo analysis
fasta_file_name_directed_digests_tss_plus = convert_bed_to_fasta(bedtools_path, genome_fasta_path, bed_file_name_directed_digests_tss_plus)

# concat FASTA files for plus and minus strand for calculation of sequence statistics and motif discovery
fasta_file_name_directed_digests_tss = fasta_file_name_directed_digests_tss_minus.split("_minus")[0] + '.fasta'
sys_cmd = 'cat ' + fasta_file_name_directed_digests_tss_minus + ' ' + fasta_file_name_directed_digests_tss_plus + ' > ' + fasta_file_name_directed_digests_tss
print("Will execute: " + sys_cmd)
os.system(sys_cmd)
fasta_file_name_directed_digests_tss_masked = mask_repeats(fasta_file_name_directed_digests_tss)

# get FASTA files with sequences around TSS on '-' strand that are associated with undirected interactions for CentriMo analysis
fasta_file_name_undirected_digests_tss_minus = convert_bed_to_fasta(bedtools_path, genome_fasta_path, bed_file_name_undirected_digests_tss_minus)

# get FASTA files with sequences around TSS on '+' strand that are associated with undirected interactions for CentriMo analysis
fasta_file_name_undirected_digests_tss_plus = convert_bed_to_fasta(bedtools_path, genome_fasta_path, bed_file_name_undirected_digests_tss_plus)

# concat FASTA files for plus and minus strand for calculation of sequence statistics and motif discovery
fasta_file_name_undirected_digests_tss = fasta_file_name_undirected_digests_tss_minus.split("_minus")[0] + '.fasta'
sys_cmd = 'cat ' + fasta_file_name_undirected_digests_tss_minus + ' ' + fasta_file_name_undirected_digests_tss_plus + ' > ' + fasta_file_name_undirected_digests_tss
print("Will execute: " + sys_cmd)
os.system(sys_cmd)
fasta_file_name_undirected_digests_tss_masked = mask_repeats(fasta_file_name_undirected_digests_tss)


print("[INFO] ... done.")


### Determine base frequencies and write to file
################################################

print("[INFO] Determining base frequencies and write to file ...")

tab_file_name_base_frequencies = out_prefix + "_base_frequencies.tab"
tab_stream_name_base_frequencies = open(tab_file_name_base_frequencies, 'wt')

# directed digests
header_line, value_line, repeat_content, gc_content_repeat, gc_content_non_repeat, gc_content_total = get_base_frequencies(fasta_file_name_directed_digests)
tab_stream_name_base_frequencies.write(header_line + "\n")
tab_stream_name_base_frequencies.write(value_line + "\n")
repeat_content_directed_digests = repeat_content
gc_content_repeat_directed_digests = gc_content_repeat
gc_content_non_repeat_directed_digests = gc_content_non_repeat
gc_content_total_directed_digests = gc_content_total

# undirected digests
header_line, value_line, repeat_content, gc_content_repeat, gc_content_non_repeat, gc_content_total = get_base_frequencies(fasta_file_name_undirected_digests)
tab_stream_name_base_frequencies.write(value_line + "\n")
repeat_content_undirected_digests = repeat_content
gc_content_repeat_undirected_digests = gc_content_repeat
gc_content_non_repeat_undirected_digests = gc_content_non_repeat
gc_content_total_undirected_digests = gc_content_total

# TSS on directed digests
header_line, value_line, repeat_content, gc_content_repeat, gc_content_non_repeat, gc_content_total = get_base_frequencies(fasta_file_name_directed_digests_tss)
tab_stream_name_base_frequencies.write(value_line + "\n")
repeat_content_directed_digests_tss = repeat_content
gc_content_repeat_directed_digests_tss = gc_content_repeat
gc_content_non_repeat_directed_digests_tss = gc_content_non_repeat
gc_content_total_directed_digests_tss = gc_content_total

# TSS on undirected digests
header_line, value_line, repeat_content, gc_content_repeat, gc_content_non_repeat, gc_content_total = get_base_frequencies(fasta_file_name_undirected_digests_tss)
tab_stream_name_base_frequencies.write(value_line + "\n")
repeat_content_undirected_digests_tss = repeat_content
gc_content_repeat_undirected_digests_tss = gc_content_repeat
gc_content_non_repeat_undirected_digests_tss = gc_content_non_repeat
gc_content_total_undirected_digests_tss = gc_content_total

tab_stream_name_base_frequencies.close()

print("[INFO] ... done.")

### Create file with interaction and digest statistics
######################################################

print("[INFO] Creating file with interaction and digest statistics ...")
print("\t[INFO] Number of directed interactions: " + str(directed_interaction_num))
print("\t[INFO] Number of undirected interactions: " + str(undirected_interaction_num))

print("\t[INFO] Number of all unique directed digests: " + str(len(directed_digests_set)))
print("\t[INFO] Number of all unique undirected digests: " + str(len(undirected_digests_set)))

print("\t[INFO] Number of unique directed digests that are not undirected: " + str(len(directed_wo_undirected_digests_set)))
print("\t[INFO] Number of unique undirected digests that are not directed: " + str(len(undirected_wo_directed_digests_set)))

print("\t[INFO] Number of intersecting digests: " + str(len(directed_intersect_undirected_digests_set)))
print("\t[INFO] Size of the union of directed and undirected digests: " + str(len(directed_union_directed_digests_set)))

connectivity_factor_directed = 1-len(directed_digests_set)/(directed_interaction_num*2) # zero, if no digest belongs to more than one interaction
print("\t[INFO] Connectivity factor of directed digests: " + str("{:.2f}".format(connectivity_factor_directed)))

connectivity_factor_undirected = 1-len(undirected_digests_set)/(undirected_interaction_num*2)
print("\t[INFO] Connectivity factor of undirected digests: " + str("{:.2f}".format(connectivity_factor_undirected)))

# Jaccard index
jaccard_index = len(directed_intersect_undirected_digests_set) / len(directed_union_directed_digests_set)
print("\t[INFO] Jaccard index: " + str("{:.2f}".format(jaccard_index)))

# Szymkiewicz–Simpson coefficient
szymkiewicz_simpson_coefficient = len(directed_intersect_undirected_digests_set) / (min(len(directed_digests_set), len(undirected_digests_set)))
print("\t[INFO] Szymkiewicz–Simpson coefficient: " + str("{:.2f}".format(szymkiewicz_simpson_coefficient)))

print("\t[INFO] Repeat content of directed digests: " + str(repeat_content_directed_digests))
print("\t[INFO] Repeat content of undirected digests: " + str(repeat_content_undirected_digests))
print("\t[INFO] Repeat content of TSS on directed digests: " + str(repeat_content_directed_digests_tss))
print("\t[INFO] Repeat content of TSS on undirected digests: " + str(repeat_content_undirected_digests_tss))

print("\t[INFO] GC content within repeat regions of directed digests: " + str(gc_content_repeat_directed_digests))
print("\t[INFO] GC content within repeat regions of undirected digests: " + str(gc_content_repeat_undirected_digests))
print("\t[INFO] GC content within repeat regions of promoters on directed digests: " + str(gc_content_repeat_directed_digests_tss))
print("\t[INFO] GC content within repeat regions of promoters on undirected digests: " + str(gc_content_repeat_undirected_digests_tss))

print("\t[INFO] GC content within non repeat regions of directed digests: " + str(gc_content_non_repeat_directed_digests))
print("\t[INFO] GC content within non repeat regions of undirected digests: " + str(gc_content_non_repeat_undirected_digests))
print("\t[INFO] GC content within non repeat regions of promoters on directed digests: " + str(gc_content_non_repeat_directed_digests_tss))
print("\t[INFO] GC content within non repeat regions of promoters on undirected digests: " + str(gc_content_non_repeat_undirected_digests_tss))

print("\t[INFO] Total GC content of directed digests: " + str(gc_content_total_directed_digests))
print("\t[INFO] Total GC content of undirected digests: " + str(gc_content_total_undirected_digests))
print("\t[INFO] Total GC content of promoters on directed digests: " + str(gc_content_total_directed_digests_tss))
print("\t[INFO] Total GC content of promoters on undirected digests: " + str(gc_content_total_undirected_digests_tss))

tab_file_name_interaction_and_digest_statistics = out_prefix + "_interaction_and_digest_statistics.tab"
tab_file_stream_interaction_and_digest_statistics = open(tab_file_name_interaction_and_digest_statistics, 'wt')

tab_file_stream_interaction_and_digest_statistics.write(
    "out_prefix" + "\t" +

    "directed_interaction_num" + "\t" +
    "undirected_interaction_num" + "\t" +

    "directed_digests_set_size" + "\t" +
    "undirected_digests_set_size" + "\t" +

    "directed_wo_undirected_digests_set_size" + "\t" +
    "undirected_wo_directed_digests_set_size" + "\t" +

    "directed_intersect_undirected_digests_set" + "\t" +
    "directed_union_directed_digests_set_size" + "\t" +

    "connectivity_factor_directed" + "\t" +
    "connectivity_factor_undirected" + "\t" +

    "jaccard_index" + "\t" +
    "szymkiewicz_simpson_coefficient" + "\t" +

    "repeat_content_directed_digests" + "\t" +
    "repeat_content_undirected_digests" + "\t" +
    "repeat_content_directed_digests_tss" + "\t" +
    "repeat_content_undirected_digests_tss" + "\t" +

    "gc_content_repeat_directed_digests" + "\t" +
    "gc_content_repeat_undirected_digests" + "\t" +
    "gc_content_repeat_directed_digests_tss" + "\t" +
    "gc_content_repeat_undirected_digests_tss" + "\t" +

    "gc_content_non_repeat_directed_digests" + "\t" +
    "gc_content_non_repeat_undirected_digests" + "\t" +
    "gc_content_non_repeat_directed_digests_tss" + "\t" +
    "gc_content_non_repeat_undirected_digests_tss" + "\t" +

    "gc_content_total_directed_digests" + "\t" +
    "gc_content_total_undirected_digests" + "\t" +
    "gc_content_total_directed_digests_tss" + "\t" +
    "gc_content_non_repeat_undirected_digests_tss" +
    "\n")
tab_file_stream_interaction_and_digest_statistics.write(
    out_prefix + "\t" +
    str(directed_interaction_num) + "\t" +
    str(undirected_interaction_num) + "\t" +
    str(len(directed_digests_set)) + "\t" +
    str(len(undirected_digests_set)) + "\t" +
    str(len(directed_wo_undirected_digests_set)) + "\t" +
    str(len(undirected_wo_directed_digests_set)) + "\t" +
    str(len(directed_intersect_undirected_digests_set)) + "\t" +
    str(len(directed_union_directed_digests_set)) + "\t" +
    str("{:.2f}".format(connectivity_factor_directed)) + "\t" +
    str("{:.2f}".format(connectivity_factor_undirected)) + "\t" +
    str("{:.2f}".format(jaccard_index)) + "\t" +
    str("{:.2f}".format(szymkiewicz_simpson_coefficient)) + "\t" +

    str(repeat_content_directed_digests) + "\t" +
    str(repeat_content_undirected_digests) + "\t" +
    str(repeat_content_directed_digests_tss) + "\t" +
    str(repeat_content_undirected_digests_tss) + "\t" +

    str(gc_content_repeat_directed_digests) + "\t" +
    str(gc_content_repeat_undirected_digests) + "\t" +
    str(gc_content_repeat_directed_digests_tss) + "\t" +
    str(gc_content_repeat_undirected_digests_tss) + "\t" +

    str(gc_content_non_repeat_directed_digests) + "\t" +
    str(gc_content_non_repeat_undirected_digests) + "\t" +
    str(gc_content_non_repeat_directed_digests_tss) + "\t" +
    str(gc_content_non_repeat_undirected_digests_tss) + "\t" +

    str(gc_content_total_directed_digests) + "\t" +
    str(gc_content_total_undirected_digests) + "\t" +
    str(gc_content_total_directed_digests_tss) + "\t" +
    str(gc_content_total_undirected_digests_tss         ) +
    "\n")

tab_file_stream_interaction_and_digest_statistics.close()