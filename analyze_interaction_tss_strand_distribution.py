import argparse
import gzip
import diachrscripts_classes as dclass
import numpy

### Parse command line
######################

parser = argparse.ArgumentParser(description='Determine expression category levels for interacting digest pairs.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--ref-gene-file', help='UCSC refGene file (must be gzipped and the same version that was used to create the digest map for Diachromatic).')
parser.add_argument('--interaction-file', help='Diachromatic interaction file.')

args = parser.parse_args()
out_prefix = args.out_prefix
ref_gene_file = args.ref_gene_file
diachromatic_interaction_file = args.interaction_file


### Define Functions
####################

def assign_strand_to_digest(digest, ref_gene_tss_map):
    strand = -1
    for i in range(digest.get_start(), digest.get_end()): # iterate digest from left to right
        key = digest.get_chromosome() + ":" + str(i)
        if ref_gene_tss_map.get_coord_strand(key) == -1: # no TSS at this position
            continue
        elif ref_gene_tss_map.get_coord_strand(key) == 'd': # TSS for '-' and '+' at this position (rare)
            return 'd'
        else:
            if strand == -1: # first TSS on digest
                strand = ref_gene_tss_map.get_coord_strand(key) # must be '-' or '+'
            elif strand != ref_gene_tss_map.get_coord_strand(key):
                return 'd'
    return strand

def categorize_interaction_strand_pair(interaction, ref_gene_tss_map):
    d1_strand = assign_strand_to_digest(interaction.get_first_digest(), ref_gene_tss_map)
    d2_strand = assign_strand_to_digest(interaction.get_second_digest(), ref_gene_tss_map)
    return str(d1_strand) + "/" + str(d2_strand)

def get_fraction(index, values):
    if 0 < sum(values):
        return float(1.0*values[index]/sum(values))
    else:
        return 0.0


### Start Execution
###################

ref_gene_tss_map = dclass.TSSCoordinateMap(ref_gene_file, "refGene") # parse refGene file with TSS

strand_simple = dclass.PairKeyDict(['+', '-', 'd', '-1'])
strand_twisted = dclass.PairKeyDict(['+', '-', 'd', '-1'])
strand_undirected = dclass.PairKeyDict(['+', '-', 'd', '-1'])
strand_undefined = dclass.PairKeyDict(['+', '-', 'd', '-1'])

n_interaction = 0
print("[INFO] Determining pair category for each interaction in " + diachromatic_interaction_file)
with gzip.open(diachromatic_interaction_file, 'r' + 't') as fp:

    line = fp.readline()

    while line:

        # parse line representing one interaction
        values = line.split("\t")

        digest_1 = dclass.Digest(values[0], int(values[1]), int(values[2]))
        if values[3] == 'A':
            digest_1.set_active()
        digest_2 = dclass.Digest(values[4], int(values[5]), int(values[6]))
        if values[7] == 'A':
            digest_2.set_active()

        if ':' in values[8]: # regular Diachromatic file
            values2 = values[8].split(":")
            n_simple = int(values2[0])
            n_twisted = int(values2[1])
            i_type = "TBD" # will be determined using binomial P-value
        else: # from LRT script
            n_simple = values[8]
            n_twisted = values[9]
            i_type = values[13].rstrip()
        interaction = dclass.Interaction(digest_1, digest_2, n_simple, n_twisted)
        interaction.set_interaction_type(i_type)

        # restrict analysis to cis long range interactions
        if not(interaction.is_cis_long_range(10000)):
            line = fp.readline()
            continue

        # restrict analysis to interactions between targeted promoters
        if interaction.get_status_pair_flag() != "AA":
            line = fp.readline()
            continue

        # assign expression level category to digest using max approach
        pair_key = categorize_interaction_strand_pair(interaction, ref_gene_tss_map)

        if interaction.get_interaction_type() == None:
            raise Exception("[FATAL] Interaction type is 'None'. This should never happen.")
        elif interaction.get_interaction_type() == "S":
            strand_simple.pair_dict[pair_key] = strand_simple.pair_dict[pair_key] + 1
        elif interaction.get_interaction_type() == "T":
            strand_twisted.pair_dict[pair_key] = strand_twisted.pair_dict[pair_key] + 1
        elif interaction.get_interaction_type() == "U":
            strand_undirected.pair_dict[pair_key] = strand_undirected.pair_dict[pair_key] + 1
        elif interaction.get_interaction_type() == "NA":
            strand_undefined.pair_dict[pair_key] = strand_undefined.pair_dict[pair_key] + 1
        else:
            line = fp.readline()
            print interaction.get_interaction_type()
            raise Exception("[FATAL] Invalid interaction type. Should be either 'S', 'T', 'U' or 'NA' but was " + interaction.get_interaction_type() + ".")

        n_interaction += 1
        if n_interaction%100 == 0:
            print "\t[INFO]", n_interaction, "interactions processed..."

        line = fp.readline()

    print("...done.")

fp.close()


### Print results to screen
###########################

print out_prefix
print "PAIR\tSIMPLE\tTWISTED\tSIMPLE+SIMPLE\tUNDIRECTED\tUNDEFINED" # absolute frequencies
for i in strand_simple.pair_dict.keys():
    print i + "\t" + str(strand_simple.pair_dict[i]) + "\t" + str(strand_twisted.pair_dict[i]) + "\t" + str(strand_simple.pair_dict[i] + strand_twisted.pair_dict[i]) + "\t" + str(strand_undirected.pair_dict[i]) + "\t" + str(strand_undefined.pair_dict[i])

print strand_undefined.pair_dict.values()
print sum(strand_undefined.pair_dict.values())
print str(float(1.0 * strand_undefined.pair_dict[i] / sum(strand_undefined.pair_dict.values())))

print "PAIR\tSIMPLE\tTWISTED\tSIMPLE+SIMPLE\tUNDIRECTED\tUNDEFINED" # relative frequencies within simple, twisted, undirected and undefined
for i in strand_simple.pair_dict.keys():

    fraction_simple = get_fraction(i,strand_simple.pair_dict)
    fraction_twisted = get_fraction(i, strand_twisted.pair_dict)
    directed = map(sum, zip(strand_simple.pair_dict,strand_simple.pair_dict))
    fraction_directed = get_fraction(i, directed)
    fraction_undirected = get_fraction(i, strand_undirected.pair_dict)
    fraction_undefined = get_fraction(i, strand_undefined.pair_dict)

    print i + "\t" + str(float(1.0 * strand_simple.pair_dict[i] / sum(strand_simple.pair_dict.values()))) + "\t" + str(float(1.0 * strand_twisted.pair_dict[i] / sum(strand_twisted.pair_dict.values()))) + "\t" \
          + str(float((1.0 * strand_simple.pair_dict[i] + strand_twisted.pair_dict[i]) / (sum(strand_simple.pair_dict.values()) + sum(strand_twisted.pair_dict.values())))) + "\t" \
          + str(float(1.0 * strand_undirected.pair_dict[i] / sum(strand_undirected.pair_dict.values()))) + "\t"# \
         # + str(float(1.0 * strand_undefined.pair_dict[i] / sum(strand_undefined.pair_dict.values())))
