#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt

# Arguments
fpkm_tracking_file_GSM661638 = sys.argv[1]   # genes.fpkm_tracking file of cuffdiff (Schoenefelder/Ter119+/GSM661638_genes.fpkm_tracking)
fpkm_tracking_file_GSM723776 = sys.argv[2]   # genes.fpkm_tracking file of cuffdiff (Schoenefelder/mESC/GSM723776_genes.fpkm_tracking)
fpkm_tracking_file_SAMEA919746 = sys.argv[3] # genes.fpkm_tracking file of cuffdiff (Schoenefelder/mESC/SAMEA919746_genes.fpkm_tracking)

def readFPKMFileToArray(name, FPKM_file):

    cnt_gene_fpkm = 0
    cnt_zero_gene_fpkm = 0
    fpkm_array = []

    with open(FPKM_file) as fp:
        line = fp.readline()

        while line:
            values = line.split("\t")

            cnt_gene_fpkm += 1

            if values[0] == "tracking_id":
                line = fp.readline()
                continue
            else:
                fpkm = float(values[9])
                if 0 < fpkm:
                    fpkm_array.append(fpkm)
                else:
                    cnt_zero_gene_fpkm += 1

            line = fp.readline()

    fp.close()

    print "Total number of gene FPKM values for", name , ":", cnt_gene_fpkm
    print "Number of zero gene FPKM values:", cnt_zero_gene_fpkm, "\n"
    return fpkm_array

# count and read fpkm values to array
fpkm_array_GSM661638 = readFPKMFileToArray("GSM661638", fpkm_tracking_file_GSM661638)
fpkm_array_GSM723776 = readFPKMFileToArray("GSM723776", fpkm_tracking_file_GSM723776)
fpkm_array_SAMEA919746 = readFPKMFileToArray("SAMEA919746", fpkm_tracking_file_SAMEA919746)

# create plot
bins = 100
f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
f.canvas.set_window_title('Distribution of FPKM values')
plt.xlabel("Distance")
plt.ylabel("Frequency")

ax1.set_title("GSM661638")
n, bins, patches = ax1.hist(np.asarray(fpkm_array_GSM661638, dtype='float'), bins, facecolor='blue', alpha=0.5, range=(1,50))
ax2.set_title("GSM723776")
n, bins, patches = ax2.hist(np.asarray(fpkm_array_GSM723776 , dtype='float'), bins, facecolor='blue', alpha=0.5, range=(1,50))
ax3.set_title("SAMEA919746")
n, bins, patches = ax3.hist(np.asarray(fpkm_array_SAMEA919746, dtype='float'), bins, facecolor='blue', alpha=0.5, range=(1,50))

plt.show()

