.. _RST_Unbalanced_interaction_calling:

##############################
Unbalanced interaction calling
##############################

We implemented the calling of unbalanced interactions in the Python script ``DICer.py``.
As input, this script requires a file in Diachromatic interaction format.
If a P-value threshold is specified, then this will be used for the classification of interactions.
Otherwise, a randomization procedure is used to set the P-value threshold such that the FDR
remains below a chosen threshold of 5\%.
The output also consists of a file in Diachromatic interaction format extended by two columns for P-values and
interaction categories. This file contains only interactions that are powered at the P-value threshold
used for classification.

****************
Using the script
****************

- ``--out-prefix <String>``
Common prefix for all generated files, which can also contain a path.
- ``--description-tag <String>``
Short description that appears in generated tables and plots.
- ``--diachromatic-interaction-file <String>``
Input file in Diachromatic interaction format.
- ``--min-inter-dist <Integer>``
Minimal interaction distance
- ``--fdr-threshold <Float>``
The P-value is chosen so that the estimated FDR remains below this threshold.
- ``--nominal-alpha-max <Float>``
Maximum nominal alpha at which iteractions are classified as significant.
- ``--nominal-alpha-step <Float>``
Step size for nominal alphas.
- ``--iter-num <Integer>``
Number of randomizations that will be performed.
- ``--random-seed <Integer>``
Random seed that is used for the first iteration. The random seed is incremented by ``1`` for each further iteration.
- ``--thread-num <Integer>``
Number of processes in which the iterations are performed in batches of the same size.
- ``--p-value-threshold <Float>``

By default, the final P-value threshold is determined via randomization. If a P-value is specified, then this P-value threshold will be used and no randomizations will be performed.

******
Output
******

DICer generates a total of seven files:

- ``DEMO_1_reports.txt``
- ``DEMO_1_randomization_plot.pdf``
- ``DEMO_1_randomization_table.txt``
- ``DEMO_1_randomization_histogram_at_threshold.pdf``
- ``DEMO_1_randomization_histogram_at_001.pdf``
- ``DEMO_1_randomization_histogram_at_005.pdf``
- ``DEMO_1_evaluated_and_categorized_interactions.tsv.gz``

*************
Randomization
*************

XXX.

*************************************
Definition of unbalanced interactions
*************************************

XXX.

********************************************
Selection of balanced reference interactions
********************************************

XXX.

