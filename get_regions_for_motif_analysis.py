"""
This script is to prepare BED files containing interacting digests and promoter regions used for motif analysis.
It iterates a file with interactions, gene symbols,
strand and TSS information, etc. created with the script 'get_gene_symbols_of_interactions_script.py'.

Six BED files are created:

   1. BED files with regions of interacting digests

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


### Prepare output files
########################

# convert aallowed pair tag strings into lists
allowed_strand_pair_tags = allowed_strand_pair_tags.split(";")
print(allowed_strand_pair_tags)
allowed_enrichment_pair_tags = str(allowed_enrichment_pair_tags).split(";")
print(allowed_enrichment_pair_tags)
allowed_interaction_categories_directed = str(allowed_interaction_categories_directed).split(";")
print(allowed_interaction_categories_directed)
allowed_interaction_categories_undirected = str(allowed_interaction_categories_undirected).split(";")
print(allowed_interaction_categories_undirected)

# create two BED files for regions of directed and undirected interactions
file_name_directed_digests = out_prefix + "_directed_digests.bed"
directed_digests_output_bed = open(file_name_directed_digests, 'wt')
file_name_undirected_digests = out_prefix + "_undirected_digests.bed"
undirected_digests_output_bed = open(file_name_undirected_digests, 'wt')

# create three BED files for regions of segmented directed and undirected interactions
file_name_directed_digests_exc = out_prefix + "_directed_digests_exc.bed"
directed_digests_output_bed_exc = open(file_name_directed_digests_exc, 'wt')
file_name_undirected_digests_exc = out_prefix + "_undirected_digests_exc.bed"
undirected_digests_output_bed_exc = open(file_name_undirected_digests_exc, 'wt')
file_name_directed_undirected_digests_inter = out_prefix + "_directed_undirected_digests_inter.bed"
directed_undirected_digests_inter = open(file_name_directed_undirected_digests_inter, 'wt')

# create two BED files for promoters on digests of directed iteractions, one for TSS on the plus and another one for TSS on the minus strand
file_name_directed_minus_tss = out_prefix + "_directed_minus_tss.bed"
directed_minus_tss_output_bed = open(file_name_directed_minus_tss, 'wt')
file_name_directed_plus_tss = out_prefix + "_directed_plus_tss.bed"
directed_plus_tss_output_bed = open(file_name_directed_plus_tss, 'wt')

# create two BED files for promoters on digests of undirected iteractions, one for TSS on the plus and another one for TSS on the minus strand
file_name_undirected_minus_tss = out_prefix + "_undirected_minus_tss.bed"
undirected_minus_tss_output_bed = open(file_name_undirected_minus_tss, 'wt')
file_name_undirected_plus_tss = out_prefix + "_undirected_plus_tss.bed"
undirected_plus_tss_output_bed = open(file_name_undirected_plus_tss, 'wt')


### Prepare variables and data structures
#########################################

chrom_sizes = {} # read chromosome sizes to hash
with open(chrom_info_file, 'rt') as fp:
    line = fp.readline()
    while line:
        chr_name = line.split("\t")[0]
        chr_size = line.split("\t")[1]
        chrom_sizes[chr_name] = int(chr_size)
        line = fp.readline()

directed_minus_tss = set()
directed_plus_tss = set()
undirected_minus_tss = set()
undirected_plus_tss = set()

directed_digests = []
directed_digests_num = 0
undirected_digests = []

cnt_first_digest_without_tss = 0
cnt_second_digest_without_tss = 0

directed_digests_set = set() # digest coordinates
undirected_digests_set = set() # digest coordinates

directed_interaction_num = 0
undirected_interaction_num = 0


### Iterate interaction file with gene symbols
##############################################

cnt=0
print("[INFO] Iterating interaction file with gene symbols " + interaction_gs_file + " ...")
with gzip.open(interaction_gs_file, 'rt') as fp:

    n_interaction_total = 0 # counter to track progress
    line = fp.readline()
    cnt = 0
    while line:

        line = line.rstrip()

        n_interaction_total += 1
        if n_interaction_total % 100000 == 0:
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
            cnt_first_digest_without_tss += 1
            line = fp.readline()
            continue

        if field[8].split(";")[1] == '':
            print("Warning: No TSS for second digest!")
            cnt_second_digest_without_tss += 1
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

        base_name = field[0] + "|" + field[1] + "|" + field[2] + "|" + field[3] + "|" + field[4]+ "|" + field[5] + "|" + field[6] + "|" + field[7]

        # get category of interacting digest
        interaction_category = field[2]

        # add digest regions to array for directed interactions
        if interaction_category in allowed_interaction_categories_directed:
            directed_digests.append(chr_a + "\t" + str(sta_a) + "\t" + str(end_a) + "\tD1|" + str(cnt) + "\n")
            directed_digests.append(chr_b + "\t" + str(sta_b) + "\t" + str(end_b) + "\tD2|" + str(cnt) + "\n")
            cnt += 1
            directed_digests_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
            directed_digests_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))
            directed_interaction_num += 1

        # add digest regions to array for undirected interactions
        if interaction_category in allowed_interaction_categories_undirected:
            undirected_digests.append(chr_a + "\t" + str(sta_a) + "\t" + str(end_a) + "\tD1|" + str(cnt) + "\n")
            undirected_digests.append(chr_b + "\t" + str(sta_b) + "\t" + str(end_b) + "\tD2|" + str(cnt) + "\n")
            cnt += 1
            undirected_digests_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
            undirected_digests_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))
            undirected_interaction_num += 1

        # get TSS associated with the first interacting digest
        tss_d1 = field[8].split(";")[0].split(",")
        tss_d1_num = len(tss_d1)

        for tss in tss_d1:
            chromosome = tss.split(":")[0]
            sta = int(tss.split(":")[1]) - up_dist
            strand = tss.split(":")[2]
            if sta < 0:
                sta = 0
            end = int(tss.split(":")[1]) + down_dist
            if chrom_sizes[chromosome] < end:
                end = chrom_sizes[chromosome]

            if interaction_category in allowed_interaction_categories_directed:
                if strand == "-":
                    directed_minus_tss.add(chromosome + "\t" + str(sta) + "\t" + str(end) + "\tD1|" + str(tss_d1_num) + "|" + base_name)
                elif strand == "+":
                    directed_plus_tss.add(chromosome + "\t" + str(sta) + "\t" + str(end) + "\tD1|" + str(tss_d1_num) + "|" + base_name)
                else:
                    print("Warning: Strand symbol of TSS was neither \'-\' not \'+\'!") # should never happen

            if interaction_category in allowed_interaction_categories_undirected:
                if strand == "-":
                    undirected_minus_tss.add(chromosome + "\t" + str(sta) + "\t" + str(end) + "\tD1|" + str(tss_d1_num) + "|" + base_name)
                elif strand == "+":
                    undirected_plus_tss.add(chromosome + "\t" + str(sta) + "\t" + str(end) + "\tD1|" + str(tss_d1_num) + "|" + base_name)
                else:
                    print("Warning: Strand symbol of TSS was neither \'-\' not \'+\'!") # should never happen

        # get TSS associated with the second interacting digest
        tss_d2 = field[8].split(";")[1].split(",")
        tss_d2_num = len(tss_d2)

        for tss in tss_d2:
            chromosome = tss.split(":")[0]
            sta = int(tss.split(":")[1]) - up_dist
            strand = tss.split(":")[2]
            if sta < 0:
                sta = 0
            end = int(tss.split(":")[1]) + down_dist
            if chrom_sizes[chromosome] < end:
                end = chrom_sizes[chromosome]

            if interaction_category in allowed_interaction_categories_directed:
                if strand == "-":
                    directed_minus_tss.add(chromosome + "\t" + str(sta) + "\t" + str(end) + "\tD2|" + str(tss_d2_num) + "|" + base_name)
                elif strand == "+":
                    directed_plus_tss.add(chromosome + "\t" + str(sta) + "\t" + str(end) + "\tD2|" + str(tss_d2_num) + "|" + base_name)
                else:
                    print("Warning: Strand symbol of TSS was neither \'-\' not \'+\'!") # should never happen

            if interaction_category in allowed_interaction_categories_undirected:
                if strand == "-":
                    undirected_minus_tss.add(chromosome + "\t" + str(sta) + "\t" + str(end) + "\tD2|" + str(tss_d2_num) + "|" + base_name)
                elif strand == "+":
                    undirected_plus_tss.add(chromosome + "\t" + str(sta) + "\t" + str(end) + "\tD2|" + str(tss_d2_num) + "|" + base_name)
                else:
                    print("Warning: Strand symbol of TSS was neither \'-\' not \'+\'!") # should never happen

        line = fp.readline()

# report interactions with no TSS on one or both digests (may happen because of imperfect enrichment annotation)
if 0 <= cnt_first_digest_without_tss:
    print("WARNING: There were " + str(cnt_first_digest_without_tss) + " interactions without TSS on the first digest!")
if 0 <= cnt_second_digest_without_tss:
    print("WARNING: There were " + str(cnt_second_digest_without_tss) + " interactions without TSS on the second digest!")

fp.close()


### Segment digests
###################

directed_wo_undirected_digests_set = directed_digests_set.difference(undirected_digests_set)
undirected_wo_directed_digests_set = undirected_digests_set.difference(directed_digests_set)
directed_intersect_undirected_digests_set = directed_digests_set.intersection(undirected_digests_set)
directed_union_directed_digests_set = directed_digests_set.union(undirected_digests_set)

print("")
print("Digest set statistics")
print("=====================")
print("Number of directed interactions: " + str(directed_interaction_num))
print("Number of undirected interactions: " + str(undirected_interaction_num))
print("Number of all unique directed digests: " + str(len(directed_digests_set)))
print("Number of all unique undirected digests: " + str(len(undirected_digests_set)))
print("Number of directed digests that are not undirected: " + str(len(directed_wo_undirected_digests_set)))
print("Number of undirected digests that are not directed: " + str(len(undirected_wo_directed_digests_set)))
print("Number of intersecting digests: " + str(len(directed_intersect_undirected_digests_set)))
print("Size of the union of directed and undirected digests: " + str(len(directed_union_directed_digests_set)))

connectivity_factor_directed = 1-len(directed_digests_set)/(directed_interaction_num*2) # zero, if no digest belongs to more than one interaction
print("Connectivity factor of directed digests: " + str("{:.2f}".format(connectivity_factor_directed)))

connectivity_factor_undirected = 1-len(undirected_digests_set)/(undirected_interaction_num*2)
print("Connectivity factor of undirected digests: " + str("{:.2f}".format(connectivity_factor_undirected)))

# Jaccard index
jaccard_index = len(directed_intersect_undirected_digests_set) / len(directed_union_directed_digests_set)
print("Jaccard index: " + str("{:.2f}".format(jaccard_index)))

# Szymkiewicz–Simpson coefficient
szymkiewicz_Simpson_coefficient = len(directed_intersect_undirected_digests_set) / (min(len(directed_digests_set), len(undirected_digests_set)))
print("Szymkiewicz–Simpson coefficient: " + str("{:.2f}".format(szymkiewicz_Simpson_coefficient)))

# determine average digest length
total_digest_length =0
for digest in directed_wo_undirected_digests_set:
    A = digest.split("\t")
    total_digest_length = total_digest_length + (int(A[2])-int(A[1]))
print("Average length of directed digests: " + str(int(total_digest_length/len(directed_wo_undirected_digests_set))) + " bp")

total_digest_length =0
for digest in undirected_wo_directed_digests_set:
    A = digest.split("\t")
    total_digest_length = total_digest_length + (int(A[2])-int(A[1]))
print("Average length of undirected digests: " + str(int(total_digest_length/len(undirected_wo_directed_digests_set))) + " bp")

total_digest_length =0
for digest in directed_intersect_undirected_digests_set:
    A = digest.split("\t")
    total_digest_length = total_digest_length + (int(A[2])-int(A[1]))
print("Average length of digests in intersect of directed and undirected digests: " + str(int(total_digest_length/len(directed_intersect_undirected_digests_set))) + " bp")

print("")

### write BED files with digest coordinates for DREME analysis
##############################################################

cnt = 1
for digest in directed_wo_undirected_digests_set:
    directed_digests_output_bed_exc.write(digest + "\t" + str(cnt) + "\n")
    cnt += 1

cnt = 1
for digest in undirected_wo_directed_digests_set:
    undirected_digests_output_bed_exc.write(digest + "\t" + str(cnt) + "\n")
    cnt += 1

directed_digests_output_bed_exc.close()
undirected_digests_output_bed_exc.close()
directed_undirected_digests_inter.close()

### Create FASTA files
######################

# get FASTA files with digest sequences of directed interactions
file_name_directed_digests_exc_base = str(file_name_directed_digests_exc).split(".b")[0]
sys_cmd = bedtools_path + ' getfasta -name -fi ' + genome_fasta_path + ' -bed  ' + file_name_directed_digests_exc + ' > ' +  file_name_directed_digests_exc_base + '.fasta'
print(sys_cmd)
os.system(sys_cmd)
sys_cmd = 'awk \'{if($1 !~ /^>/){gsub(/a|c|g|t/,"N",$1)};print}\' ' + file_name_directed_digests_exc_base + '.fasta' + ' > ' + file_name_directed_digests_exc_base + '_masked.fasta'
print(sys_cmd)
os.system(sys_cmd)

# get FASTA files with digest sequences of undirected interactions
file_name_undirected_digests_exc_base = str(file_name_undirected_digests_exc).split(".b")[0]
sys_cmd = bedtools_path + ' getfasta -name -fi ' + genome_fasta_path + ' -bed  ' + file_name_undirected_digests_exc + ' > ' +  file_name_undirected_digests_exc_base + '.fasta'
print(sys_cmd)
os.system(sys_cmd)
sys_cmd = 'awk \'{if($1 !~ /^>/){gsub(/a|c|g|t/,"N",$1)};print}\' ' + file_name_undirected_digests_exc_base + '.fasta' + ' > ' + file_name_undirected_digests_exc_base + '_masked.fasta'
print(sys_cmd)
os.system(sys_cmd)


print("")
print("Digest sequence statistics")
print("==========================")
print("")
print("GC content of directed digests before repeat masking:")
a_occ_total = 0
c_occ_total = 0
g_occ_total = 0
t_occ_total = 0
with open(file_name_directed_digests_exc_base + '.fasta', 'rt') as fp:
    line = fp.readline()
    while line:
        if line.count('>') == 0:
            line = line.rstrip().upper()
            a_occ_total = a_occ_total + line.count('A')
            c_occ_total = c_occ_total + line.count('C')
            g_occ_total = g_occ_total + line.count('G')
            t_occ_total = t_occ_total + line.count('T')
        line = fp.readline()

acgt_occ_total = a_occ_total + c_occ_total + g_occ_total + t_occ_total
print("   A: " + str(a_occ_total) + " (" + str(100*round(a_occ_total / acgt_occ_total,4)) + ")")
print("   C: " + str(c_occ_total) + " (" + str(100*round(c_occ_total / acgt_occ_total,4)) + ")")
print("   G: " + str(g_occ_total) + " (" + str(100*round(g_occ_total / acgt_occ_total,4)) + ")")
print("   T: " + str(t_occ_total) + " (" + str(100*round(t_occ_total / acgt_occ_total,4)) + ")")
print("--------------------------------------------------------------------------")
at_occ_total = a_occ_total + t_occ_total
print("   AT: " + str(at_occ_total) + " (" + str(100*round(at_occ_total / acgt_occ_total,4)) + ")")
gc_occ_total = g_occ_total + c_occ_total
print("   GC: " + str(gc_occ_total) + " (" + str(100*round(gc_occ_total / acgt_occ_total,4)) + ")")
fp.close()

print("")
print("GC content of directed digests after repeat masking:")
a_occ_total = 0
c_occ_total = 0
g_occ_total = 0
t_occ_total = 0
with open(file_name_directed_digests_exc_base + '_masked.fasta', 'rt') as fp:
    line = fp.readline()
    while line:
        if line.count('>') == 0:
            line = line.rstrip().upper()
            a_occ_total = a_occ_total + line.count('A')
            c_occ_total = c_occ_total + line.count('C')
            g_occ_total = g_occ_total + line.count('G')
            t_occ_total = t_occ_total + line.count('T')
        line = fp.readline()

acgt_occ_total = a_occ_total + c_occ_total + g_occ_total + t_occ_total
print("   A: " + str(a_occ_total) + " (" + str(100*round(a_occ_total / acgt_occ_total,4)) + ")")
print("   C: " + str(c_occ_total) + " (" + str(100*round(c_occ_total / acgt_occ_total,4)) + ")")
print("   G: " + str(g_occ_total) + " (" + str(100*round(g_occ_total / acgt_occ_total,4)) + ")")
print("   T: " + str(t_occ_total) + " (" + str(100*round(t_occ_total / acgt_occ_total,4)) + ")")
print("--------------------------------------------------------------------------")
at_occ_total = a_occ_total + t_occ_total
print("   AT: " + str(at_occ_total) + " (" + str(100*round(at_occ_total / acgt_occ_total,4)) + ")")
gc_occ_total = g_occ_total + c_occ_total
print("   GC: " + str(gc_occ_total) + " (" + str(100*round(gc_occ_total / acgt_occ_total,4)) + ")")
fp.close()

print("")
print("GC content of undirected digests before repeat masking:")
a_occ_total = 0
c_occ_total = 0
g_occ_total = 0
t_occ_total = 0
with open(file_name_undirected_digests_exc_base + '.fasta', 'rt') as fp:
    line = fp.readline()
    while line:
        if line.count('>') == 0:
            line = line.rstrip().upper()
            a_occ_total = a_occ_total + line.count('A')
            c_occ_total = c_occ_total + line.count('C')
            g_occ_total = g_occ_total + line.count('G')
            t_occ_total = t_occ_total + line.count('T')
        line = fp.readline()

acgt_occ_total = a_occ_total + c_occ_total + g_occ_total + t_occ_total
print("   A: " + str(a_occ_total) + " (" + str(100*round(a_occ_total / acgt_occ_total,4)) + ")")
print("   C: " + str(c_occ_total) + " (" + str(100*round(c_occ_total / acgt_occ_total,4)) + ")")
print("   G: " + str(g_occ_total) + " (" + str(100*round(g_occ_total / acgt_occ_total,4)) + ")")
print("   T: " + str(t_occ_total) + " (" + str(100*round(t_occ_total / acgt_occ_total,4)) + ")")
print("--------------------------------------------------------------------------")
at_occ_total = a_occ_total + t_occ_total
print("   AT: " + str(at_occ_total) + " (" + str(100*round(at_occ_total / acgt_occ_total,4)) + ")")
gc_occ_total = g_occ_total + c_occ_total
print("   GC: " + str(gc_occ_total) + " (" + str(100*round(gc_occ_total / acgt_occ_total,4)) + ")")
fp.close()

print("")
print("GC content of undirected digests after repeat masking:")
a_occ_total = 0
c_occ_total = 0
g_occ_total = 0
t_occ_total = 0
with open(file_name_undirected_digests_exc_base + '_masked.fasta', 'rt') as fp:
    line = fp.readline()
    while line:
        if line.count('>') == 0:
            line = line.rstrip().upper()
            a_occ_total = a_occ_total + line.count('A')
            c_occ_total = c_occ_total + line.count('C')
            g_occ_total = g_occ_total + line.count('G')
            t_occ_total = t_occ_total + line.count('T')
        line = fp.readline()

acgt_occ_total = a_occ_total + c_occ_total + g_occ_total + t_occ_total
print("   A: " + str(a_occ_total) + " (" + str(100*round(a_occ_total / acgt_occ_total,4)) + ")")
print("   C: " + str(c_occ_total) + " (" + str(100*round(c_occ_total / acgt_occ_total,4)) + ")")
print("   G: " + str(g_occ_total) + " (" + str(100*round(g_occ_total / acgt_occ_total,4)) + ")")
print("   T: " + str(t_occ_total) + " (" + str(100*round(t_occ_total / acgt_occ_total,4)) + ")")
print("--------------------------------------------------------------------------")
at_occ_total = a_occ_total + t_occ_total
print("   AT: " + str(at_occ_total) + " (" + str(100*round(at_occ_total / acgt_occ_total,4)) + ")")
gc_occ_total = g_occ_total + c_occ_total
print("   GC: " + str(gc_occ_total) + " (" + str(100*round(gc_occ_total / acgt_occ_total,4)) + ")")
fp.close()

print("")
print("Repeat content of directed digests:")
n_occ_total = 0
acgtn_occ_total = 0
with open(file_name_directed_digests_exc_base + '_masked.fasta', 'rt') as fp:
    line = fp.readline()
    while line:
        if line.count('>') == 0:
            line = line.rstrip()
            n_occ_total = n_occ_total + line.count('N')
            acgtn_occ_total = acgtn_occ_total + len(line)
        line = fp.readline()

print(100*round(n_occ_total/acgtn_occ_total,4))
fp.close()
print("")
print("Repeat content of undirected digests:")
n_occ_total = 0
acgtn_occ_total = 0
with open(file_name_undirected_digests_exc_base + '_masked.fasta', 'rt') as fp:
    line = fp.readline()
    while line:
        if line.count('>') == 0:
            line = line.rstrip()
            n_occ_total = n_occ_total + line.count('N')
            acgtn_occ_total = acgtn_occ_total + len(line)
        line = fp.readline()

print(100*round(n_occ_total/acgtn_occ_total,4))
fp.close()

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
        if n_interaction_total % 100000 == 0:
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
            cnt_first_digest_without_tss += 1
            line = fp.readline()
            continue

        if field[8].split(";")[1] == '':
            print("Warning: No TSS for second digest!")
            cnt_second_digest_without_tss += 1
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
    if chrom_sizes[chr] < end:
        end = chrom_sizes[chr] - 1
    strand = A[2]
    if strand == '-':
        directed_minus_tss_output_bed.write(chr + "\t" + str(sta) + "\t" + str(end) + "\t" + str(cnt_minus) + "\n")
        cnt_minus += 1
    elif strand == '+':
        directed_plus_tss_output_bed.write(chr + "\t" + str(sta) + "\t" + str(end) + "\t" + str(cnt_plus) + "\n")
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
    if chrom_sizes[chr] < end:
        end = chrom_sizes[chr] - 1
    strand = A[2]
    if strand == '-':
        undirected_minus_tss_output_bed.write(chr + "\t" + str(sta) + "\t" + str(end) + "\t" + str(cnt_minus) + "\n")
        cnt_minus += 1
    elif strand == '+':
        undirected_plus_tss_output_bed.write(chr + "\t" + str(sta) + "\t" + str(end) + "\t" + str(cnt_plus) + "\n")
        cnt_plus += 1
    else:
        print("Warning: Strand of TSS was neither \'-\' nor \'+\'!")

directed_minus_tss_output_bed.close()
undirected_minus_tss_output_bed
undirected_minus_tss_output_bed.close()
undirected_plus_tss_output_bed.close()

# convert BED to FASTA
# get FASTA files with sequences around TSS on '+' strand that are associated with directed interactions
file_name_directed_digests_base = str(file_name_directed_plus_tss).split(".b")[0]
sys_cmd = bedtools_path + ' getfasta -name -fi ' + genome_fasta_path + ' -bed  ' + file_name_directed_plus_tss + ' > ' + file_name_directed_digests_base + '.fasta'
print(sys_cmd)
os.system(sys_cmd)

# get FASTA files with sequences around TSS on '+' strand that are associated with undirected interactions
file_name_undirected_digests_base = str(file_name_undirected_plus_tss).split(".b")[0]
sys_cmd = bedtools_path + ' getfasta -name -fi ' + genome_fasta_path + ' -bed  ' + file_name_undirected_plus_tss + ' > ' + file_name_undirected_digests_base + '.fasta'
print(sys_cmd)
os.system(sys_cmd)

# get FASTA files with sequences around TSS on '-' strand that are associated with directed interactions
file_name_directed_digests_base = str(file_name_directed_minus_tss).split(".b")[0]
sys_cmd = bedtools_path + ' getfasta -name -fi ' + genome_fasta_path + ' -bed  ' + file_name_directed_minus_tss + ' > ' + file_name_directed_digests_base + '.fasta'
print(sys_cmd)
os.system(sys_cmd)

# get FASTA files with sequences around TSS on '-' strand that are associated with undirected interactions
file_name_undirected_digests_base = str(file_name_undirected_minus_tss).split(".b")[0]
sys_cmd = bedtools_path + ' getfasta -name -fi ' + genome_fasta_path + ' -bed  ' + file_name_undirected_minus_tss + ' > ' + file_name_undirected_digests_base + '.fasta'
print(sys_cmd)
os.system(sys_cmd)

exit(0)

# END# END# END# END# END# END# END# END# END# END# END# END# END# END# END# END

# write BED file with digest regions for motif analysis of directed interactions
directed_digests_num = 0
for digest in directed_digests:
    directed_digests_output_bed.write(digest)
    directed_digests_num += 1
undirected_digests_num = 0
for digest in undirected_digests:
    if undirected_digests_num <= directed_digests_num:
        undirected_digests_output_bed.write(digest)
        undirected_digests_num += 1


### write BED files for CentriMo analysis
#########################################

total_number_of_directed_promoters_minus = 0
for promoter in directed_minus_tss:
    field = promoter.split("\t")
    directed_minus_tss_output_bed.write(promoter + "|" + str(total_number_of_directed_promoters_minus) + "\n")
    total_number_of_directed_promoters_minus += 1
total_number_of_directed_promoters_plus = 0
for promoter in directed_plus_tss:
    field = promoter.split("\t")
    directed_plus_tss_output_bed.write(promoter + "|" + str(total_number_of_directed_promoters_plus) + "\n")
    total_number_of_directed_promoters_plus += 1

total_number_of_undirected_promoters_minus = 0
for promoter in undirected_minus_tss:
    field = promoter.split("\t")
    undirected_minus_tss_output_bed.write(promoter + "|" + str(total_number_of_undirected_promoters_minus) + "\n")
    total_number_of_undirected_promoters_minus += 1
    if total_number_of_directed_promoters_minus <= total_number_of_undirected_promoters_minus:
        break
total_number_of_undirected_promoters_plus = 0
for promoter in undirected_plus_tss:
    field = promoter.split("\t")
    undirected_plus_tss_output_bed.write(promoter + "|" + str(total_number_of_undirected_promoters_plus) + "\n")
    total_number_of_undirected_promoters_plus += 1
    if total_number_of_directed_promoters_plus <= total_number_of_undirected_promoters_plus:
        break


### Close output BED files
##########################

directed_digests_output_bed.close() # for motif discovery with DREME
undirected_digests_output_bed.close() # for motif discovery with DREME

directed_minus_tss_output_bed.close() # for motif enrichment analysis with CentriMo
directed_plus_tss_output_bed.close() # for motif enrichment analysis with CentriMo

undirected_minus_tss_output_bed.close() # for motif enrichment analysis with CentriMo
undirected_plus_tss_output_bed.close() # for motif enrichment analysis with CentriMo


### Create FASTA files
######################

# get FASTA files with digest sequences of directed interactions
file_name_directed_digests_base = str(file_name_directed_digests).split(".b")[0]
sys_cmd = bedtools_path + ' getfasta -name -fi ' + genome_fasta_path + ' -bed  ' + file_name_directed_digests + ' > ' +  file_name_directed_digests_base + '.fasta'
print(sys_cmd)
os.system(sys_cmd)
sys_cmd = 'awk \'{if($1 !~ /^>/){gsub(/a|c|g|t/,"N",$1)};print}\' ' + file_name_directed_digests_base + '.fasta' + ' > ' + file_name_directed_digests_base + '_masked.fasta'
print(sys_cmd)
os.system(sys_cmd)

# get FASTA files with digest sequences of undirected interactions
file_name_undirected_digests_base = str(file_name_undirected_digests).split(".b")[0]
sys_cmd = bedtools_path + ' getfasta -name -fi ' + genome_fasta_path + ' -bed  ' + file_name_undirected_digests + ' > ' +  file_name_undirected_digests_base + '.fasta'
print(sys_cmd)
os.system(sys_cmd)
sys_cmd = 'awk \'{if($1 !~ /^>/){gsub(/a|c|g|t/,"N",$1)};print}\' ' + file_name_undirected_digests_base + '.fasta' + ' > ' + file_name_undirected_digests_base + '_masked.fasta'
print(sys_cmd)
os.system(sys_cmd)

# get FASTA files with sequences around TSS on '+' strand that are associated with directed interactions
file_name_directed_digests_base = str(file_name_directed_plus_tss).split(".b")[0]
sys_cmd = bedtools_path + ' getfasta -name -fi ' + genome_fasta_path + ' -bed  ' + file_name_directed_plus_tss + ' > ' + file_name_directed_digests_base + '.fasta'
print(sys_cmd)
os.system(sys_cmd)

# get FASTA files with sequences around TSS on '+' strand that are associated with undirected interactions
file_name_undirected_digests_base = str(file_name_undirected_plus_tss).split(".b")[0]
sys_cmd = bedtools_path + ' getfasta -name -fi ' + genome_fasta_path + ' -bed  ' + file_name_undirected_plus_tss + ' > ' + file_name_undirected_digests_base + '.fasta'
print(sys_cmd)
os.system(sys_cmd)

# get FASTA files with sequences around TSS on '-' strand that are associated with directed interactions
file_name_directed_digests_base = str(file_name_directed_minus_tss).split(".b")[0]
sys_cmd = bedtools_path + ' getfasta -name -fi ' + genome_fasta_path + ' -bed  ' + file_name_directed_minus_tss + ' > ' + file_name_directed_digests_base + '.fasta'
print(sys_cmd)
os.system(sys_cmd)

# get FASTA files with sequences around TSS on '-' strand that are associated with undirected interactions
file_name_undirected_digests_base = str(file_name_undirected_minus_tss).split(".b")[0]
sys_cmd = bedtools_path + ' getfasta -name -fi ' + genome_fasta_path + ' -bed  ' + file_name_undirected_minus_tss + ' > ' + file_name_undirected_digests_base + '.fasta'
print(sys_cmd)
os.system(sys_cmd)
