#!/usr/bin/env python

"""
With this script all plots can be created that are also created in the following Jupyter notebook:

       jupyter_notebooks/binomialModel.ipynb

All functionality is in class BinomialInteractionModel that is also used by the notebook.

In between, interactive plot windows open that have to closed for the script to continue.
"""

import os
import sys
import argparse
from diachr import BinomialInteractionModel
from diachr import EnhancedInteraction
from diachr import EnhancedInteractionParser
import matplotlib.pyplot as plt
import scipy, scipy.stats, numpy

parser = argparse.ArgumentParser(description='Explore binomial model using simuated and real interactions.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--diachromatic-interaction-file', help='Enhanced interaction file.')

args = parser.parse_args()
out_prefix = args.out_prefix
diachromatic_interaction_file = args.enhanced_interaction_file

print("[INFO] " + "Input parameters")
print("\t[INFO] --out_prefix: " + out_prefix)
print("\t[INFO] --diachromatic-interaction-file: " + str(diachromatic_interaction_file))

# Create BinomialInteractionModel object
bim = BinomialInteractionModel()
bim._out_prefix = out_prefix

# Verification of the implementation of the P-value calculation through simulation
# --------------------------------------------------------------------------------

# Simulation for n=0,...,400 and a P-value threshold of 0.05
n_max = 400
i_num = 20000000
pvt = 0.05
n_list, n_sim_list, n_sig_list = bim.simulate_interactions(n_max=n_max, i_num=i_num, pvt=pvt)
bim.simulate_interactions_plot(n_max=n_max, i_num=i_num, pvt=pvt,
                               n_list=n_list,n_sim_list=n_sim_list,n_sig_list=n_sig_list, CREATE_PDF=True)


# Simualtion for n=0,...,400 and a P-value threshold of 0.0038
n_max = 400
i_num = 100000000
pvt = 0.0038
n_list, n_sim_list, n_sig_list = bim.simulate_interactions(n_max=n_max, i_num=i_num, pvt=pvt)
bim.simulate_interactions_plot(n_max=n_max, i_num=i_num, pvt=pvt,
                               n_list=n_list,n_sim_list=n_sim_list,n_sig_list=n_sig_list, CREATE_PDF=True)


plt.close('all')

# Binomial distributions for different n and a fixed P-value threshold
# --------------------------------------------------------------------

# Binomial distributions for n=1,...,10 and a P-value threshold of 0.05
pvt=0.05
N=10
bim.analyze_N_binomial_distributions_with_fixed_p_thresh(N=N, pvt=pvt, CREATE_DIST_PLOTS=False, CREATE_PDF=True)

# Comparison with the results from the simualtion with n=1,...10 and a P-value threshold of 0.05
n_max = 10
i_num = 100000000
pvt = 0.05
n_list, n_sim_list, n_sig_list = bim.simulate_interactions(n_max=n_max, i_num=i_num, pvt=pvt)
bim.simulate_interactions_plot(n_max=n_max, i_num=i_num, pvt=pvt,
                               n_list=n_list,n_sim_list=n_sim_list,n_sig_list=n_sig_list, CREATE_PDF=True)

plt.close('all')

# Binomial distributions for n=1,...,40 and a P-value threshold of 0.05
pvt=0.05
N=40
bim.analyze_N_binomial_distributions_with_fixed_p_thresh(N=N, pvt=pvt, CREATE_DIST_PLOTS=False, CREATE_PDF=True)



# Comparison with the results from the simualtion with n=1,...40 and a P-value threshold of 0.05
n_max = 40
i_num = 100000000
pvt = 0.05
n_list, n_sim_list, n_sig_list = bim.simulate_interactions(n_max=n_max, i_num=i_num, pvt=pvt)

bim.simulate_interactions_plot(n_max=n_max, i_num=i_num, pvt=pvt,
                               n_list=n_list,n_sim_list=n_sim_list,n_sig_list=n_sig_list, CREATE_PDF=True)


# Binomial distributions for n=1,...,400 and a P-value threshold of 0.05
pvt=0.05
N=400
bim.analyze_N_binomial_distributions_with_fixed_p_thresh(N=N, pvt=pvt, CREATE_DIST_PLOTS=False, CREATE_PDF=True)



# Binomial distributions for n=1,...,400 and a P-value threshold of 0.0038
pvt=0.0038
N=400
bim.analyze_N_binomial_distributions_with_fixed_p_thresh(N=N, pvt=pvt, CREATE_DIST_PLOTS=False, CREATE_PDF=True)


# Distribution of n for empirical data
# ------------------------------------

if diachromatic_interaction_file == None:
    exit(0)

if not os.path.exists(diachromatic_interaction_file):
    raise FileNotFoundError("Could not find interaction file")
n_list, n_di_list, n_uir_list, n_ui_list = bim.count_di_uir_and_ui_for_each_n(ei_file=diachromatic_interaction_file)

# Plot distribution for n=1,...,400
x_max = 400
bim.count_di_uir_and_ui_for_each_n_plot(n_list=n_list, n_di_list=n_di_list, n_uir_list=n_uir_list, n_ui_list=n_ui_list, x_max = x_max, l_wd=0.0)

# Plot distribution for n=1,...,20
x_max = 20
y_max = 0.04
l_wd=0.2
bim.count_di_uir_and_ui_for_each_n_plot(n_list=n_list, n_di_list=n_di_list, n_uir_list=n_uir_list, n_ui_list=n_ui_list, x_max = x_max, y_max = y_max, l_wd=l_wd)

# Compare empirical and theoretical distribution for n=1,...,20
pvt=0.0038
N=20
bim.analyze_N_binomial_distributions_with_fixed_p_thresh(N=N, pvt=pvt, CREATE_DIST_PLOTS=False, CREATE_PDF=True)
