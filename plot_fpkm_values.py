#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt

# Arguments
fpkm_tracking_file_GSM661638 = sys.argv[1]   # genes.fpkm_tracking file of cuffdiff (Schoenefelder/Ter119+/GSM661638_genes.fpkm_tracking)
fpkm_tracking_file_GSM723776 = sys.argv[2]   # genes.fpkm_tracking file of cuffdiff (Schoenefelder/mESC/GSM723776_genes.fpkm_tracking)
fpkm_tracking_file_SAMEA919746 = sys.argv[3] # genes.fpkm_tracking file of cuffdiff (Schoenefelder/mESC/SAMEA919746_genes.fpkm_tracking)
fpkm_tracking_file_GSM758572 = sys.argv[4]   # genes.fpkm_tracking file of cuffdiff (Mifsud/GM12878/GSM758572_genes.fpkm_tracking)
fpkm_tracking_file_GSM758559 = sys.argv[5] # genes.fpkm_tracking file of cuffdiff (Mifsud/GM12878/GSM758559_genes.fpkm_tracking)

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
fpkm_array_GSM758572 = readFPKMFileToArray("GSM758572", fpkm_tracking_file_GSM758572    )
fpkm_array_GSM758559 = readFPKMFileToArray("GSM758559", fpkm_tracking_file_GSM758559)

# create plot
bins = 100
f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, sharex=True, sharey=True)
f.canvas.set_window_title('Distribution of FPKM values')
plt.xlabel("Distance")
plt.ylabel("Frequency")
ax1.set_title("GSM661638 (Schoenefelder)")
n, bins, patches = ax1.hist(np.asarray(fpkm_array_GSM661638, dtype='float'), bins, facecolor='blue', alpha=0.5, range=(1,50))
ax2.set_title("GSM723776 (Schoenefelder)")
n, bins, patches = ax2.hist(np.asarray(fpkm_array_GSM723776 , dtype='float'), bins, facecolor='blue', alpha=0.5, range=(1,50))
ax3.set_title("SAMEA919746 (Schoenefelder)")
n, bins, patches = ax3.hist(np.asarray(fpkm_array_SAMEA919746, dtype='float'), bins, facecolor='blue', alpha=0.5, range=(1,50))
ax4.set_title("GSM758572 (Mifsud)")
n, bins, patches = ax4.hist(np.asarray(fpkm_array_GSM758572 , dtype='float'), bins, facecolor='blue', alpha=0.5, range=(1,50))
ax5.set_title("GSM758559 (Mifsud)")
n, bins, patches = ax5.hist(np.asarray(fpkm_array_GSM758559, dtype='float'), bins, facecolor='blue', alpha=0.5, range=(1,50))
plt.show()

# print quartiles to screen
print "Quartiles for GSM661638 (Schoenefelder):", np.quantile(fpkm_array_GSM661638, (.25, .5,.75))
print "Quartiles for GSM723776 (Schoenefelder):", np.quantile(fpkm_array_GSM723776, (.25, .5,.75))
print "Quartiles for SAMEA919746 (Schoenefelder):", np.quantile(fpkm_array_SAMEA919746, (.25, .5,.75))
print "Quartiles for GSM758572 (Mifsud):", np.quantile(fpkm_array_GSM758572, (.25, .5,.75))
print "Quartiles for GSM758559 (Mifsud):", np.quantile(fpkm_array_GSM758559, (.25, .5,.75))


