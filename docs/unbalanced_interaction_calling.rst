.. _RST_Unbalanced_interaction_calling:

##############################
Unbalanced interaction calling
##############################

DICer will be described here.

****************
Using the script
****************

XXX.

******
Output
******

XXX.

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

*****************************************
Randomization of simple and twisted (old)
*****************************************

To decide whether directed interactions occur more often than expected by chance,
we took an approach in which we randomize the
simple and twisted read pair counts of individual interactions.

For a given dataset,
we first determine the number of interactions
that are significant at a chosen P-value threshold,
using our binomial test for directionality of interactions.
Then we randomize the simple ``s`` and twisted read pair counts ``t``
for each interaction according to our null model,
i.e. we randomly draw a number of simple read pairs ``s'``
from a binomial distribution with ``n=s+t`` and ``p=0.5``.
We then set the number of twisted read pairs to ``t'=n-s'``.
After we randomized the counts of all interactions,
we determine the number of significant interactions.
If we repeat this procedure ``x`` times and never observe more
significant interactions than for the non-randomized data,
we speak of an empirical P-value of ``1/x``.
In addition, we keep track of the number of randomized
significant interactions for each iteration and,
at the end, calculate mean, standard deviation and Z-score.

We implemented this analysis in the following script:

.. code-block:: console

    $ python diachrscripts/02_perform_permutation_analysis.py
       --interaction-file MK/gzdir/JAV_MK_R10.interaction.counts.table.clr_20000.tsv.gz
       --nominal-alpha 0.05
       --iter-num 10000
       --thread-num 20
       --out-prefix MK/JAV_MK_R10

In this example,
we are applying the script to the first replicate for MK (``--interaction-file``).
From this file, trans-interactions, short interactions and interactions with and on ``chrM``
were filtered out.
We use a P-value threshold of 0.05 (``--nominal-alpha``) and perform 10,000 iterations (``--iter-num``).
Multi-threading is supported (``--thread-num``).
In this example,
500 iterations are performed in each of the 20 threads
and the numbers of significant interactions are subsequently combined.

The script generates two files:

.. code-block:: console

    MK/JAV_MK_R10/JAV_MK_R10_permutation_stats.txt
    MK/JAV_MK_R10/JAV_MK_R10_permutation_w_sig_p.txt

The first file contains information the chosen parameters and all calculated values
in form of a table with a header row and another row with corresponding values.

+-------------------------------------+------------------------------------------------+
| Header                              | Values                                         |
+=====================================+================================================+
| OUT_PREFIX                          | MK/JAV_MK_R10                                  |
+-------------------------------------+------------------------------------------------+
| ITER_NUM                            | 10000                                          |
+-------------------------------------+------------------------------------------------+
| NOMINAL_ALPHA                       | 0.05                                           |
+-------------------------------------+------------------------------------------------+
| INDEF_RP_CUTOFF                     | 5                                              |
+-------------------------------------+------------------------------------------------+
| N_INTERACTION                       | 73430640                                       |
+-------------------------------------+------------------------------------------------+
| N_INDEFINABLE_INTERACTION           | 70094483                                       |
+-------------------------------------+------------------------------------------------+
| N_UNDIRECTED_INTERACTION            | 3004191                                        |
+-------------------------------------+------------------------------------------------+
| N_DIRECTED_INTERACTION              | 331966                                         |
+-------------------------------------+------------------------------------------------+
| MEAN_PERMUTED_DIRECTED_INTERACTION  | 177028.25                                      |
+-------------------------------------+------------------------------------------------+
| SD_PERMUTED_DIRECTED_INTERACTION    | 397.56                                         |
+-------------------------------------+------------------------------------------------+
| Z_SCORE                             | 389.72                                         |
+-------------------------------------+------------------------------------------------+
| N_PERMUTED_BETTER_THAN_OBSERVED     | 0                                              |
+-------------------------------------+------------------------------------------------+

The second file contains the numbers of randomized significant interactions,
with each row corresponding to one iteration.
These are the first lines of the file for the example above:

.. code-block:: console

    $ head MK/JAV_MK_R10/JAV_MK_R10_permutation_w_sig_p.txt
    177200
    176372
    177082
    177010

**************************
False discovery rate (old)
**************************

This script uses a simple procedure to estimate the FDR of directed interactions for increasing P-value thresholds.

Initially, a Diachromatic interaction file is ingested and a binomial P-value is calculated for each interaction and
stored in a list. Furthermore, the total numbers of interactions for given read pair numbers ``n`` (the sum of simple
and twisted read pairs) are stored in dictionary with ``n`` as keys and interaction numbers as values.
This dictionary is then used to generate a list of random P-values for all interactions. A random P-value for an
interaction with n read pairs is generated by drawing a number of simple read pairs s' from a binomial distribution
``B(n, p = 0.5)``, setting the number of twisted read pairs to ``t' = n - s'`` and calculating the corresponding binomial
P-value.
Finally, the FDR is estimated for increasing P-value thresholds. For each threshold, the number of significant
interactions is determined for the list of original P-values (``S_o``) and for the lists of randomized P-values
(``S_p``) and ``S_p / S_o`` is used as estimator for the FDR.

We implemented this procedure in the following script:

.. code-block:: console

    $ python diachrscripts/03_perform_fdr_analysis.py
       --diachromatic-interaction-file MK/gzdir/JAV_MK_R10.interaction.counts.table.clr_20000.tsv.gz
       --p-val-c-min 0.00025
       --p-val-c-max 0.05
       --p-val-step-size 0.00025
       --fdr-threshold 0.05
       --out-prefix MK/JAV_MK_RALT/FDR/JAV_MK_RALT

In this example,
we are applying the script to a Diachromatic interaction file for the first replicate for MK
(``--diachromatic-interaction-file``).
The FDR is calculated for increasing P-value thresholds
from ``--p-val-c-min``
to ``--p-val-c-max``
with a step size of ``--p-val-step-size``.
All calculated FDRs and the corresponding P-value thresholds are written to
a tab separated file (``--out-prefix``)
and the line with the desired FDR (``--fdr-threshold``)
is additionally output on the screen.

The script generates one file:

.. code-block:: console

    MK/JAV_MK_RALT/FDR/JAV_MK_RALT_fdr_analysis_results.txt

These are the first lines of this file up to the line
in which the calculated FDR is above the specified FDR threshold (``0.05``).
The last line in which the calculated FDR is still below the threshold
was subsequently marked with an arrow.

.. code-block:: console

    OUT_PREFIX	FDR	PC	NSIG_P	NSIG_O
    MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.003878727859788114	0.0001	376	96939
    MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.006966033903740457	0.0002	782	112259
    MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.011531095177042278	0.00030000000000000003	1457	126354
    MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.013328508533864062	0.0004	1768	132648
    MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.018910363781799288	0.0005	2766	146269
    MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.02002288329519451	0.0006000000000000001	3010	150328
    MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.021868671886084275	0.0007000000000000001	3391	155062
    MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.024751345785056698	0.0008	4014	162173
    MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.025583259451778622	0.0009000000000000001	4213	164678
    MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.03238307526066565	0.001	5724	176759
    MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.03334187538996346	0.0011	5985	179504
    MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.03588880531188952	0.0012000000000000001	6632	184793
    MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.03923472957590568	0.0013000000000000002	7514	191514
    MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.039965068029495815	0.0014000000000000002	7734	193519
    MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.040836216309993977	0.0015	8003	195978
    MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.041607002400016135	0.0016	8252	198332
    MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.04285223179365562	0.0017000000000000001	8601	200713
    MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.04672194378811623	0.0018000000000000002	9630	206113
    -> MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.04791701522502928	0.0019000000000000002	10024	209195
    MK/JAV_MK_RALT/FDR/JAV_MK_RALT	0.05701037416737231	0.002	12590	220837
    ...

The table has four columns:

+----------------+----------------------------------------------------------------------+
| Header         | Explanation                                                          |
+================+======================================================================+
| ``OUT_PREFIX`` | Prefix for output files                                              |
+----------------+----------------------------------------------------------------------+
| ``FDR``        | Calculated FDR at the given P-value threshold (``NSIG_P/NSIG_O``)    |
+----------------+----------------------------------------------------------------------+
| ``PC``         | P-value threshold                                                    |
+----------------+----------------------------------------------------------------------+
| ``NSIG_P``     | Number of randomized significant interactions at given ``PC``        |
+----------------+----------------------------------------------------------------------+
| ``NSIG_O``     | Number of significant interactions actually observed at given ``PC`` |
+----------------+----------------------------------------------------------------------+

Note that ``_P`` stands for permuted and ``_O`` for observed.

******************************************
Evaluate and categorize interactions (old)
******************************************

For each interaction, we calculate a P-value for directionality of interaactions
from a two sided binomial test of the simple and twisted read pair counts.
Based on a given P-value threshold, we then divide the interactions into directed (``DI``) and undirected (``UI``).
Finally, we select undirected reference interactions (``UIR``) from ``UI``, which are comparable to ``DI`` with respect to
to the total number of read pairs (the sum of simple and twisted read pairs).
The input consists of a file in Diachromatic interaction format and a P-value threshold that was determined using our
FDR procedure. The output is again a file in Diachromatic interaction format, but there are two additional columns
on the right for the calculated P-values and interaction categories.

Workflow
========

We have implemented all steps related to the evaluation and categorization of interactions
in class ``DiachromaticInteractionSet``, which is instantiated as follows:

.. code-block:: console

    interaction_set = DiachromaticInteractionSet()

Into this object, interactions can be read in from one or more Diachromatic interaction files.

.. code-block:: console

    interaction_set.parse_file("data/test_04/diachromatic_interaction_file.tsv")

The object ``interaction_set`` then contains a dictionary with objects of class ``DiachromaticInteraction``.
Note that the interactions are unique in terms of their coordinates,
i.e. if the input contains two interactions with identical coordinates,
only one ``DiachromaticInteraction`` object is created,
with the simple and twisted read pair counts added up separately.

The next step is to calculate the P-values of all interactions in ``interaction_set`` and,
based on this and a specified threshold, divide interactions into directed (``DI``) and undirected (``UI``).
This step is carried out in the function ``DiachromaticInteractionSet.rate_and_categorize_interactions``,
which expects the negative natural logarithm of a P-value threshold as input.

.. code-block:: console

    nln_p_val_thresh = -log(0.01)
    rate_and_cat_report_dict = interaction_set.rate_and_categorize_interactions(nln_p_val_thresh)

This step discards interactions that cannot be significant at the given P-value threshold because they
do not have enough read pairs.
All other interactions are assigned either the category ``DI``,
if they have a P-value that is less than or equal to the threshold value,
or the category ``UI``,
if the P-value is greater than the threshold value.
In addition, the function returns a dictionary which contains parameters,
intermediate results and summary statistics.
This dictionary can be used to generate a report on the execution of the function or to test the function.
The dictionary contains the following key-value pairs:

+-----------------------+----------------------------------------------------------+
| Key                   | Value                                                    |
+=======================+==========================================================+
| ``NLN_PVAL_THRESH``   | Chosen P-value threshold                                 |
+-----------------------+----------------------------------------------------------+
| ``MIN_RP``            | Minimum number of read pairs required for significance   |
+-----------------------+----------------------------------------------------------+
| ``MIN_RP_PVAL``       | Largest possible significant P-value                     |
+-----------------------+----------------------------------------------------------+
| ``N_PROCESSED``       | Total number of interactions processed                   |
+-----------------------+----------------------------------------------------------+
| ``N_DISCARDED``       | Number of discarded interactions                         |
+-----------------------+----------------------------------------------------------+
| ``N_UNDIRECTED``      | Number of undirected interactions                        |
+-----------------------+----------------------------------------------------------+
| ``N_DIRECTED``        | Number of directed interactions                          |
+-----------------------+----------------------------------------------------------+

The next step is to select undirected reference interactions that are comparable to
the directed interactions with respect to the number of read pairs per interaction.
This step is carried out in the function ``DiachromaticInteractionSet.select_reference_interactions``.

.. code-block:: console

    select_ref_report_dict = interaction_set.select_reference_interactions()

After this step, interactions that were selected as reference have a ``UIR`` tag
instead of a ``UI`` tag.
The function also returns a dictionary that contains the counts of interactions in the
various categories.
Note that, when selecting reference interactions, we also take into account the enrichment status
of digests involved in an interaction,
where ``E`` stands for enriched and ``N`` for not enriched.
The dictionary conntains the following key-value pairs:

+---------------------------------------------------------+-----------------------------------------------+
| Key                                                     | Value                                         |
+=========================================================+===============================================+
| ``DI_NN``, ``DI_NE``, ``DI_EN``, ``DI_EE``              | Numbers of directed interactions              |
+---------------------------------------------------------+-----------------------------------------------+
| ``UIR_NN``, ``UIR_NE``, ``UIR_EN``, ``UIR_EE``          | Numbers of undirected reference interactions  |
+---------------------------------------------------------+-----------------------------------------------+
| ``M_UIR_NN``, ``M_UIR_NE``, ``M_UIR_EN``, ``M_UIR_EE``  | Numbers of missing reference interactions     |
+---------------------------------------------------------+-----------------------------------------------+
| ``UI_NN``, ``UI_NE``, ``UI_EN``, ``UI_EE``              | Numbers of undirected interactions            |
+---------------------------------------------------------+-----------------------------------------------+

Finally, the function ``DiachromaticInteractionSet.write_diachromatic_interaction_file`` can be used
to write all interactions from a ``DiachromaticInteractionSet`` object into a Diachromatic interaction file.

.. code-block:: console

    interaction_set.write_diachromatic_interaction_file("evaluated_and_categorized_interactions.tsv.gz")

Since the interactions were previously evaluated and categorized,
the created file contains two additional columns for the P-value and interaction category on the left.

Script
======

XXX

Test file
=========


We have prepared a Diachromatic interaction file to test
the selection of undirected reference interactions:

.. code-block:: console

    tests/data/test_04/diachromatic_interaction_file.tsv

This file contains 18 interactions that are clearly directed,
with three interactions in the enrichment category ``NN``,
four in ``NE``, five in ``EN`` and six in ``EE``.
The total number of read pairs per interaction ranges from 101 to 106.

.. code-block:: console

    chr14	43059116	43059494	N	chr14	43101212	43101810	N	100:1
    chr8	129042054	129044258	N	chr8	129121269	129121986	N	100:2
    chr15	73467156	73468652	N	chr15	73526903	73528438	N	100:3

    chr17	72411026	72411616	N	chr17	72712662	72724357	E	100:1
    chr18	38724804	38726198	N	chr18	76794986	76803172	E	100:2
    chr11	114362648	114362686	N	chr11	114396073	114404234	E	100:3
    chr15	56158017	56158267	N	chr15	56462978	56465983	E	100:4

    chr14	34714080	34716362	E	chr14	50135355	50139051	N	100:1
    chr1	91022201	91023797	E	chr1	116561813	116566655	N	100:2
    chr1	15681566	15697108	E	chr1	19411358	19417940	N	100:3
    chr12	13586326	13591414	E	chr12	101116206	101119138	N	100:4
    chr2	113676580	113686263	E	chr2	202796295	202797013	N	100:5

    chr7	25228385	25228778	E	chr7	42234764	42240281	E	100:1
    chr6	117448209	117455412	E	chr6	152981947	152991467	E	100:2
    chr7	123787495	123793134	E	chr7	141015018	141017643	E	100:3
    chr10	100185111	100188716	E	chr10	100911854	100914842	E	100:4
    chr10	73301814	73302637	E	chr10	112212615	112226039	E	100:5
    chr8	11066997	11069837	E	chr8	11473293	11485455	E	100:6

In addition, there are 16 interactions that are clearly undirected and
match the directed interactions in terms of enrichment category and
total number of read pairs per interaction.
These are the undirected reference interactions.

.. code-block:: console

    chr1	68750175	68752699	N	chr1	68795323	68796512	N	50:51
    chr10	24701259	24703782	N	chr10	24993210	24994662	N	50:52
    chr13	77661563	77666713	N	chr13	77716416	77717610	N	50:53

    chr11	86392415	86393959	N	chr11	125833282	125834157	E	50:51
    chr4	143079578	143081309	N	chr4	146945582	146955674	E	50:52
    chr14	97520062	97521143	N	chr14	99641271	99654714	E	50:53

    chr12	92715804	92721658	E	chr12	97575775	97576200	N	50:51
    chr3	149470667	149475062	E	chr3	150148068	150150290	N	50:52
    chr18	31373121	31377878	E	chr18	53036132	53037663	N	50:53
    chr2	16841469	16846902	E	chr2	17018322	17019298	N	50:54
    chr14	92748029	92749630	E	chr14	96578956	96580826	N	50:55

    chr3	185935734	185943372	E	chr3	194132402	194139770	E	50:51
    chr5	65953726	65970373	E	chr5	71582365	71588637	E	50:52
    chr5	138607182	138612232	E	chr5	151368315	151371525	E	50:53
    chr1	23789605	23791509	E	chr1	153951062	153962812	E	50:54
    chr5	112154251	112162055	E	chr5	115769718	115774925	E	50:55

The file also contains 12 interactions,
which are also clearly undirected,
but the number of read pair per interaction is in a different range (50 to 52),
which is why these interactions cannot be used as reference interactions.

.. code-block:: console

    chr8	110169057	110171420	N	chr8	110203244	110203772	N	25:25
    chr21	24706442	24708010	N	chr21	24802647	24806351	N	25:26
    chr1	21038219	21041845	N	chr1	21139524	21140811	N	25:27

    chr7	105237629	105238078	N	chr7	128862293	128872009	E	25:25
    chr6	89716777	89716882	N	chr6	89916736	89920828	E	25:26
    chr9	120719242	120724112	N	chr9	135466848	135469617	E	25:27

    chr14	24167761	24173500	E	chr14	75079406	75087746	N	25:25
    chr12	132700658	132715264	E	chr12	133104166	133106284	N	25:26
    chr8	142668614	142673077	E	chr8	142903855	142904484	N	25:27

    chr17	80806561	80813742	N	chr17	80880355	80886553	N	25:25
    chr2	197200527	197202476	N	chr2	197253179	197256110	N	25:26
    chr2	139005552	139005761	N	chr2	139026311	139029719	N	25:27

For two directed interactions,
there is no matching undirected reference interaction.
One interaction is in enrichment category ``NE`` and has a total of 104 read pairs:

.. code-block:: console

    chr15	56158017	56158267	N	chr15	56462978	56465983	E	100:4

The second interaction is in the enrichment category ``EE`` and has a total of 106 read pairs:

.. code-block:: console

    chr8	11066997	11069837	E	chr8	11473293	11485455	E	100:6
