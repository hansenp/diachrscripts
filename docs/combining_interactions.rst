.. _RST_Combining_interactions:

#################################################
Pooling of interactions from different replicates
#################################################

For the dataset on the hematopoietic cell types, there are three to four biological replicates
for each cell type.
In order to pool interactions from different replicates,
we discard interactions that occur in less than two replicates and,
for the remaining interactions, we add up the read pair counts from
all replicates separately for the four counts.

Using the script
================

We implemented the pooling of interactions from different replicates in the following script:

.. code-block:: console

    $ python diachrscripts/additional_scripts/pooler.py
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
the four read pair counts from different replicates
will be added up separately.
For this analysis,
we required that an interaction must occur in at least two replicates.
The name of each created file will have the same prefix ``--out-prefix``,
which can also contain the path to an already existing directory.
For pooled replicates,
we have used the abbreviation ``RALT``,
where ``ALT`` stands for **A**\ t\  **L**\ east **T**\ wo.

The command above will generate the following two files:

.. code-block:: console

    MK/JAV_MK_RALT_at_least_in_2_replicates_summary.txt
    MK/JAV_MK_RALT_at_least_in_2_replicates_interactions.tsv.gz

The first file contains an overview of the numbers of interactions
in the individual files and
the second file contains the pooled interactions.


Testing the script
==================

Diachromatic
even outputs interactions that have only a single read pair.
On the other hand, when pooling interactions,
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

These are the contents of the four interaction files:

.. code-block:: console

    # REPLICATE 1
    chr1    46297999    46305684    E    chr1    51777391    51781717    N    1:1:1:0

    # REPLICATE 2
    chr1    46297999    46305684    E    chr1    51777391    51781717    N    2:0:1:0
    chr17   72411026    72411616    N    chr17   72712662    72724357    N    3:0:1:1

    # REPLICATE 3
    chr1    46297999    46305684    E    chr1    51777391    51781717    N    0:2:1:0
    chr17   72411026    72411616    N    chr17   72712662    72724357    N    3:0:0:2
    chr7    69513952    69514636    N    chr7    87057837    87061499    E    3:1:1:2

    # REPLICATE 4
    chr1    46297999    46305684    E    chr1    51777391    51781717    N    1:1:1:0
    chr17   72411026    72411616    N    chr17   72712662    72724357    N    3:0:2:0
    chr7    69513952    69514636    N    chr7    87057837    87061499    E    2:2:2:1
    chr11   47259263    47272706    N    chr11   91641153    91642657    E    3:2:1:3

To get the pooled interactions that occur in at least two replicates,
execute the following command:

.. code-block:: console

    $ python diachrscripts/additional_scripts/pooler.py \
       --interaction-files-path tests/data/test_01/ \
       --required-replicates 2
       --out-prefix TEST \

This is the content of the generated file with the pooled interactions:

.. code-block:: console

    chr1    46297999    46305684    E    chr1    51777391    51781717    N    4:4:4:0
    chr17   72411026    72411616    N    chr17   72712662    72724357    N    9:0:3:3
    chr7    69513952    69514636    N    chr7    87057837    87061499    E    5:3:3:3

The interaction on chromosome ``chr11`` does not occur in this file
because it was observed for replicate 4 only,
but we require that an interaction occurs in at least two replicates.

The interaction on chromosome ``chr7`` occurs in the files for replicate 3 and 4.

.. code-block:: console

    chr7    69513952    69514636    N    chr7    87057837    87061499    E    3:1:1:2 (R3)
    chr7    69513952    69514636    N    chr7    87057837    87061499    E    2:2:2:1 (R4)
    ------------------------------------------------------------------------------------------
    chr7    69513952    69514636    N    chr7    87057837    87061499    E    5:3:3:3 (POOLED)

The interaction on chromosome ``chr17`` occurs in the files for replicate 2, 3 and 4.

.. code-block:: console

    chr17   72411026    72411616    N    chr17   72712662    72724357    N    3:0:1:1 (R2)
    chr17   72411026    72411616    N    chr17   72712662    72724357    N    3:0:0:2 (R3)
    chr17   72411026    72411616    N    chr17   72712662    72724357    N    3:0:2:0 (R4)
    ------------------------------------------------------------------------------------------
    chr17   72411026    72411616    N    chr17   72712662    72724357    N    9:0:3:3 (POOLED)

Finally, the interaction on ``chr1`` occurs in the files for all four replicates.

.. code-block:: console

    chr1    46297999    46305684    E    chr1    51777391    51781717    N    1:1:1:0 (R1)
    chr1    46297999    46305684    E    chr1    51777391    51781717    N    2:0:1:0 (R2)
    chr1    46297999    46305684    E    chr1    51777391    51781717    N    0:2:1:0 (R3)
    chr1    46297999    46305684    E    chr1    51777391    51781717    N    1:1:1:0 (R4)
    ------------------------------------------------------------------------------------------
    chr1    46297999    46305684    E    chr1    51777391    51781717    N    4:4:4:0 (POOLED)
