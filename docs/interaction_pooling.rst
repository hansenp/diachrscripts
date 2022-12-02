.. _RST_Interaction_pooling:

###################
Interaction pooling
###################


Write a short text about what the script can be used for


Running the script ``pooler.py``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the following command to run the script: ::

    $ diachrscripts/pooler.py \
       --interaction-files-path gzdir/ \
       --required-replicates 2 \
       --out-prefix out_dir/out_prefix

Available arguments:

+---------------+-------------------------------+-------------------------+-----------+---------------------------------------------------------------------------------------------+-----------------+
| Short option  | Long option                   | Example                 | Required  | Description                                                                                 | Default         |
+===============+===============================+=========================+===========+=============================================================================================+=================+
| ``-i``        | ``--interaction-files-path``  | ``gzdir``               | yes       | Path to directory with Diachromatic interaction files                                       | --              |
+---------------+-------------------------------+-------------------------+-----------+---------------------------------------------------------------------------------------------+-----------------+
| ``-i``        | ``--required-replicates``     | ``3``                   | no        | Interactions that occur in fewer files will be discarded                                    | ``2``           |
+---------------+-------------------------------+-------------------------+-----------+---------------------------------------------------------------------------------------------+-----------------+
| ``-o``        | ``--out-prefix``              | ``out_dir/out_prefix``  | no        | Prefix for output files, which can also contain the path to an already existing directory.  | ``OUT_PREFIX``  |
+---------------+-------------------------------+-------------------------+-----------+---------------------------------------------------------------------------------------------+-----------------+

Output files
~~~~~~~~~~~~

The script produces a Diachromatic interaction file containing the pooled interactions:

    * ``out_dir/out_prefix_at_least_in_2_replicates_interactions.tsv.gz``

In addition, a file is produced that contains summary statistics.

    * ``out_dir/out_prefix_at_least_in_2_replicates_summary.txt``

..
    ******************
    Testing the script
    ******************

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

        $ python diachrscripts/pooler.py \
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
