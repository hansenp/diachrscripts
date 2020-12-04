"""
In this script, the P-values for directionality of interactions are calculated and, based on these P-values,
the interactions are categorized into directed (DI) und undirected interactions (UI).
Undirected reference interactions (UIR) are then selected from UI, which are comparable to DI with respect to
to the total number of read pairs (n).

The input consists of a file in Diachromatic interaction format and a P-value threshold that was determined using our
FDR procedure. The output is again a file in Diachromatic interaction format, but there are two additional columns
on the right for the calculated P-values and interaction categories.
"""

from diachr.binomial_model import BinomialModel
import time

# PARSE COMMAND LINE

# TEST MODULE FOR P-VALUE CALCULATION

p_vaule_factory = BinomialModel()

start_time = time.time()
n_max = 1000
for n in range(1, n_max + 1):
    for k in range(0, n + 1):
        p_vaule_factory.get_binomial_logsf_p_value(k, n - k)
end_time = time.time()
pval_num = p_vaule_factory._get_pval_dict_size()

print("It took " + "{:.2f}".format(end_time-start_time) + " seconds to calculate " + str(pval_num) + " P-values.")

start_time = time.time()
n_max = 1000
for n in range(1, n_max + 1):
    for k in range(0, n + 1):
        p_vaule_factory.get_binomial_logsf_p_value(k, n - k)
end_time = time.time()
pval_num = p_vaule_factory._get_pval_dict_size()

print("It took " + "{:.2f}".format(end_time-start_time) + " seconds to get " + str(pval_num) + " P-values from the dictionary.")


# LOAD INTERACTIONS WITH 'DiachromaticInteractionparser'

# FIRST PASS: CALCULATE P-VALUES AND DEFINE DIRECTED AND UNDIRECTED INTERACTIONS

# SECOND PASS: SELECT UNDIRECTED REFERENCE INTERACTIONS

# WRITE FILE IN DIACHROMATIC INTERACTION FORMAT WITH TWO ADDITIONAL COLUMNS


