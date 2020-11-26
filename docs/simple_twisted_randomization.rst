.. _RST_simple_twisted_randomization:

###################################
Randomization of simple and twisted
###################################

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

We implemented this analysis in the follwing script:

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
were filtered out
(:ref:`RST_Interaction_calling`).
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
