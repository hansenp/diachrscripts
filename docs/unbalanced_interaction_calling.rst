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

+--------------------------------------------------------------------------------------------+------------------------------------------------------------------------+
| Argument                                                                                   | Meaning                                                                |
+============================================================================================+========================================================================+
| --out-prefix <String>                                                                      | Common prefix for all generated files, which can also include a path.  |
+--------------------------------------------------------------------------------------------+------------------------------------------------------------------------+
| --description-tag <String>                                                                 | Short description tag that appears in generated tables and plots.      |
+--------------------------------------------------------------------------------------------+------------------------------------------------------------------------+
| --diachromatic-interaction-file <String>                                                   | Path to an input file in diachromatic interaction format.              |
+--------------------------------------------------------------------------------------------+------------------------------------------------------------------------+
| --min-inter-dist <Integer>                                                                 | Minimum interaction distance.                                          |
+--------------------------------------------------------------------------------------------+------------------------------------------------------------------------+
| --read-pair-counts-rule <String>                                                           | Rule by which interactions are scored.                                 |
+--------------------------------------------------------------------------------------------+------------------------------------------------------------------------+
| By default, the sum of the two highest counts is compared against                          |                                                                        |
+--------------------------------------------------------------------------------------------+------------------------------------------------------------------------+
| the sum of the two lowest counts using a binomial test (ht).                               |                                                                        |
+--------------------------------------------------------------------------------------------+------------------------------------------------------------------------+
| Alternatively, the sum of counts for type 0 and type 1 read pairs can be compared against  |                                                                        |
+--------------------------------------------------------------------------------------------+------------------------------------------------------------------------+
| the sum of type 2 and type 3 read pairs (st).                                              |                                                                        |
+--------------------------------------------------------------------------------------------+------------------------------------------------------------------------+



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

