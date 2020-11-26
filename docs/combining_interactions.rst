.. _RST_Combining_interactions:

######################
Combining interactions
######################

For the dataset that we analyzed, there are three to four biological replicates
for each cell type.
On the other hand, we need as many read pairs as possible in the individual
interactions in order to decide whether an interaction is directed or not.
We therefore decided on an compromise when combining interactions from different
replicates.
We discard interactions that occur in less than two replicates and,
for the remaining interactions, we add up the read pair counts from
all replicates separately for simple and twisted read pairs.

Using the script
================

This way
of combining interactions from different replicates is implemented in
the following script:

.. code-block:: console

    $ python diachrscripts/01_combine_interactions_from_replicates.py
       --interaction-files-path MK/gzdir
       --required-replicates 2
       --out-prefix MK/JAV_MK_RALT

The script expects a path to a directory that contains gzipped files in Diachromatic's interaction format
(``--interaction-files-path``).
From the files output by Diachromatic,
we have filtered out interactions between different chromosomes (trans),
interactions with a distance of less than 20,000 bp and
interactions with or on chromosome ``chrM``
(see :ref:`RST_Interaction_calling`).
This is the content of the directory with the gzipped files for ``MK``:

.. code-block:: console

    $ ls MK/gzdir
    JAV_MK_R10.interaction.counts.table.clr_20000.tsv.gz
    JAV_MK_R20.interaction.counts.table.clr_20000.tsv.gz
    JAV_MK_R30.interaction.counts.table.clr_20000.tsv.gz
    JAV_MK_R40.interaction.counts.table.clr_20000.tsv.gz

In addition, the required number of replicates must be specified (``--required-replicates``).
All interactions that occur in less replicates
will be discarded.
For the remaining interactions,
the simple and twisted read pair counts from different replicates
will be added up separately.
For this analysis,
we required that an interaction must occur in at least two replicates.
The name of each created file will have the same prefix ``--out-prefix``, which can also contain the path to an already existing directory.
For combined replicates,
we have used the abbreviation ``RALT``,
where ``ALT`` stands for **A**\ t\  **L**\ east **T**\ wo.

The command above will generate the following two files:

.. code-block:: console

    MK/JAV_MK_RALT_at_least_in_2_replicates_summary.txt
    MK/JAV_MK_RALT_at_least_in_2_replicates_interactions.tsv.gz

The first file contains an overview of the numbers of interactions
in the individual files and
the second file contains the combined interactions.


Testing the script
==================

Diachromatic
even outputs interactions that have only a single read pair.
On the other hand, when combining interactions,
the interactions from multiple replicates must be read into memory.
Therefore, the memory consumption can become very high
and we carried out this step on a compute cluster.

We have prepared small input files for testing
so that this step can be followed here.
There are a total of four replicates and four interactions.
The first interaction occurs only in replicate 1,
the second interaction occurs in replicates 1 and 2,
the third interaction occurs in replicates 1,2 and 3 and
the fourth interactions occurs in all replicates.

This is content of the interaction file for replicate 4:

.. code-block:: console

    chr1    46297999   46305684   A   chr1    51777391   51781717   I   2:1
    chr17   72411026   72411616   I   chr17   72712662   72724357   I   3:2
    chr7    69513952   69514636   I   chr7    87057837   87061499   A   4:3
    chr11    9641153    9642657   I   chr11   47259263   47272706   A   5:4

From this file, we created the files for the other three replicates
by deleting interactions from the last line one by one.
By creating the files in this way,
the individual interactions have same simple and twisted read pair counts
for all replicates, which is usually not the case.
However, it simplifies the presentation here,
because we only need to know the content of the file for replicate 4
in order to understand the results.

To get the combined interactions that occur in at least two replicates,
execute the following command:

.. code-block:: console

    $ python diachrscripts/01_combine_interactions_from_replicates.py \
       --interaction-files-path tests/data/test_01/ \
       --required-replicates 2
       --out-prefix TEST \

This is the content of the generated file with the combined interactions:

.. code-block:: console

    chr1    46297999   46305684   A   chr1    51777391   51781717   I   8:4
    chr17   72411026   72411616   I   chr17   72712662   72724357   I   9:6
    chr7    69513952   69514636   I   chr7    87057837   87061499   A   8:6

The interaction on chromosome ``chr11`` does not occur in this file
because it was observed for replicate 4 only.
However, we required that an interaction must have been observed in at least two replicates.
The interaction on chromosome ``chr7`` occurs in the files for replicate 3 and 4.
Since this interaction has the same read pair counts for both replicates,
the counts in the file for combined interactions double
(``4:3`` becomes ``8:6``).
The interaction on chromosome ``chr17`` occurs in the files for replicate 2, 3 and 4
and the counts triple (``3:2`` becomes ``9:6``).
Finally, the interaction on ``chr11`` occurs in the files for all four replicates
and the counts quadruple (``2:1`` becomes ``8:4``).