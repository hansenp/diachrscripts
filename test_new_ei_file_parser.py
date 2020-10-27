import os
import sys
import gzip
import math
PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))
from diachr2 import EnhancedInteraction, EnhancedInteractionParser
import diachrscripts_toolkit

# The purpose of this file is to demonstrate that we can switch to using the EnhancedInteractionParser in the diachr2 package
# instead of the "line.split" solution in scripts such as diachrscripts/06_select_uir_from_uie_and_uii.py
# We first use the class to parse out a list of lines -- each line is represented as an EnhancedInteraction object
# 
# To run the test, just pass the path to an enhancer interaction file
# python test_new_ei_file_parser.py <path/to/ie-file> 


if len(sys.argv) < 2:
    print('You need to specify the path to the interaction file')
    sys.exit()

input_path = sys.argv[1]


## With the class
ie_parser = EnhancedInteractionParser(input_path)
ei_list = ie_parser.parse()
print("[INFO] The EnhancedInteractionParser extracted %d interaction lines" % len(ei_list))

## With the diachrscripts_toolkit.parse_enhanced_interaction_line_with_gene_symbols(line) function

def raise_error(line, msg, item_n, item_o):
    print("[ERROR]", line)
    print("New '%s' Old: '%s'" % (item_n, item_o))
    raise ValueError(msg)


def check_dif(line1, line2):
    fields1 = line1.split('\t')
    fields2 = line2.split('\t')
    if len(fields1) != len(fields2):
        raise ValueError("Fields1 len=x and field2 len=y")
    for i in range(len(fields2)):
        if fields1[i] != fields2[i]:
            if isinstance(fields1[i], float):
                if float(fields1[i]) == float(fields2[i]):
                    continue
                else:
                    raise ValueError("i=%d f1=%s f2=%s" % (i, fields1[i], fields2[i]))

i = 0
with gzip.open(input_path, 'rt') as fp:
    for line in fp:
        line = line.rstrip()
        # Parse enhanced interactions line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total, i_dist = \
            diachrscripts_toolkit.parse_enhanced_interaction_line_with_gene_symbols(line)
        ei = ei_list[i]
        my_newline = ei.get_line()
        check_dif(line, my_newline)
        
            #raise ValueError("Lines did not match")
        i += 1
        if ei.chr_a != chr_a:
            raise ValueError("chr_a did not match at (%d)" % i)
        if ei.sta_a != sta_a:
            raise ValueError("sta_a did not match at (%d)" % i)
        if ei.end_a != end_a:
            raise ValueError("end_a did not match at (%d)" % i)
        if ei.syms_a != syms_a:
            raise ValueError("syms_a did not match at (%d)" % i)
        if ei.tsss_a != tsss_a:
            raise ValueError("tsss_a did not match at (%d)" % i)
        if ei.chr_b != chr_b:
            raise ValueError("chr_b did not match at (%d)" % i)
        if ei.sta_b != sta_b:
            raise ValueError("sta_b did not match at (%d)" % i)
        if ei.end_b != end_b:
            raise ValueError("end_b did not match at (%d)" % i)
        if ei.tsss_b != tsss_b:
            raise_error(line, "tsss_b did not match at (%d)" % i, ei.tsss_b, tsss_b)
        if ei.enrichment_pair_tag != enrichment_pair_tag:
            raise ValueError("enrichment_pair_tag did not match at (%d)" % i)
        if ei.strand_pair_tag != strand_pair_tag:
            raise ValueError("strand_pair_tag did not match at (%d)" % i)
        if ei.interaction_category != interaction_category:
            raise ValueError("interaction_category did not match at (%d)" % i)
        if ei.neg_log_p_value != neg_log_p_value:
            raise ValueError("neg_log_p_value did not match at (%d)" % i)
        if ei.rp_total != rp_total:
            raise ValueError("rp_total did not match at (%d)" % i)
        if ei.i_dist != i_dist:
            raise ValueError("i_dist did not match at (%d)" % i)


print("[INFO] We compared %d lines and found NO DIFFERENCE" % i)
