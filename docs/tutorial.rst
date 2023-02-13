.. _RST_tutorial:

########
Tutorial
########

This flowchart provides an overview of the entire workflow, which we will work through step by step below.

.. image:: img/analysis_flowchart.png

***************************
Downloading paired-end data
***************************

Use the script
`dumpy.sh <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/additional_scripts/dumpy.sh>`_
to download paired-end data from promoter capture Hi-C experiments in GM12878 cells
`(Mifsud et al. 2015) <https://pubmed.ncbi.nlm.nih.gov/25938943/>`_.
For this dataset, ``200 GB`` hard disk space must be available.


.. code-block:: console

    $ diachrscripts/additional_scripts/dumpy.sh MIF_GM12878_CHC_REP1 MIF_GM12878_CHC_REP1 "ERR436029"
    $ diachrscripts/additional_scripts/dumpy.sh MIF_GM12878_CHC_REP2 MIF_GM12878_CHC_REP2 "ERR436028;ERR436030;ERR436033"
    $ diachrscripts/additional_scripts/dumpy.sh MIF_GM12878_CHC_REP3 MIF_GM12878_CHC_REP3 "ERR436026;ERR436031"


A separate directory is created for each of the three replicates.
Because this is paired-end data, there are pairs of forward and reverse FASTQ files with
suffixes ``_1`` and ``_2``.
Both files have the same number of lines and reads with the same line index form a pair.

.. code-block:: console

    $ ls -lh MIF_GM12878_CHC_REP*/
    MIF_GM12878_CHC_REP1/:
    total 28G
    12G MIF_GM12878_CHC_REP1_1.fastq.gz
    12G MIF_GM12878_CHC_REP1_2.fastq.gz
    242 MIF_GM12878_CHC_REP1_md5.txt

    MIF_GM12878_CHC_REP2/:
    total 87G
    36G MIF_GM12878_CHC_REP2_1.fastq.gz
    37G MIF_GM12878_CHC_REP2_2.fastq.gz
    490 MIF_GM12878_CHC_REP2_md5.txt

    MIF_GM12878_CHC_REP3/:
    total 57G
    24G MIF_GM12878_CHC_REP3_1.fastq.gz
    24G MIF_GM12878_CHC_REP3_2.fastq.gz
    366 MIF_GM12878_CHC_REP3_md5.txt

The MD5 checksums are as follows:

.. code-block:: console

    $ cat MIF_GM12878_CHC_REP*/MIF_GM12878_CHC_REP*_md5.txt | grep -v ERR
    c68ab47d9f8cb9b49c1701ff04c7fd41  MIF_GM12878_CHC_REP1/MIF_GM12878_CHC_REP1_1.fastq.gz
    338980ac0a0bc1a3ea26bc6cf55dfb43  MIF_GM12878_CHC_REP1/MIF_GM12878_CHC_REP1_2.fastq.gz
    b6c248b48ded9eee8c4e951f4ea8d4a4  MIF_GM12878_CHC_REP2/MIF_GM12878_CHC_REP2_1.fastq.gz
    a3ff3dfde21d39dc41dc3c1212cbb3e3  MIF_GM12878_CHC_REP2/MIF_GM12878_CHC_REP2_2.fastq.gz
    5ceed4dd8784a467afa35f00b83d132d  MIF_GM12878_CHC_REP3/MIF_GM12878_CHC_REP3_1.fastq.gz
    4260c37ed24ee29c5046ef74096ba1f7  MIF_GM12878_CHC_REP3/MIF_GM12878_CHC_REP3_2.fastq.gz

******************************************
Calling interactions with ``Diachromatic``
******************************************

We have extended our previously published
`Diachromatic software <https://diachromatic.readthedocs.io/en/latest/index.html>`__
for initial processing and quality control of Hi-C
and CHi-C data to report read pair counts of interactions separately for the four relative paired-end orientations.
Diachromatic processes the data in three steps:

1. Truncation of chimeric reads ``>`` Truncated reads
2. Mapping and artifact removal ``>`` Valid mapped read pairs
3. Counting supporting reads pairs of restriction fragment pairs ``>`` Interactions

Truncation of chimeric reads
============================

Paired-end sequencing of Hi-C libraries can result in chimeric reads containing sequences from different
genomic regions. Such reads cannot be mapped to the reference genome. Therefore, they must be truncated so that they
are no longer chimeric. This can be done using ``Diachromatic`` with the subcommand ``truncate``.
The chimeric reads must be cut at the ligation sites,
which is why the restriction enzyme used for the experiment must be specified (``-e``).
The prepared input FASTQ files with the forward and reverse reads are specified using ``-q`` and -``-r``.

.. code-block:: console

    $ java -jar Diachromatic.jar truncate \
       -e HindIII \
       -q MIF_GM12878_CHC_REP1/MIF_GM12878_CHC_REP1_1.fastq.gz \
       -r MIF_GM12878_CHC_REP1/MIF_GM12878_CHC_REP1_2.fastq.gz \
       -o MIF_GM12878_CHC_REP1 \
       -x MIF_GM12878_CHC_REP1

The result files are written to the directory specified using ``-o`` and have the same prefix specified using ``-x``.

Mapping and artifact removal
============================

For mapped Hi-C paired-end reads, no particular distribution of distances (insert sizes) or relative orientation
can be assumed.
However, read mappers typically rely on a minimum and maximum insert size and that mapped read pairs point inwards.
Therefore, the truncated forward and reverse reads must be mapped independently and the mapped reads must be re-paired
afterwards.
In addition, there are certain rules by which artifact read pairs that are specific to Hi-C data can be recognized
and removed.
This can be done using ``Diachromatic`` with the subcommand ``align``. We recommend having ``16`` to ``32 GB``
memory available.
A path to a bowtie2 executable (``-b``) and an index for a corresponding reference sequence (``-i``) must be specified.
If the ``-bsu`` switch is used, reads are considered to be mapped uniquely if they map to only one location.
The ``-p`` option can be used to specify how many CPUs will be used by bowtie2.
For the detection of artifact read pairs, a digest file (``-d``) is required,
which contains all restriction fragments resulting from a complete digestion of the genome.
The FASTQ files with the truncated forward (``-q``) and reverse reads (``-r``) must be specified.

Before the following command can be executed,
the bowtie2 index and the digest map must first be prepared.
How to do this is documented here: :ref:`RST_Diachromatic_input_preparation`.

.. code-block:: console

    $ java -Xmx32000m -jar Diachromatic.jar align \
       -b <BOWTIE2_EXECUTABLE> \
       -i <BOWTIE2_INDEX_PATH>/genome \
       -bsu \
       -p 4 \
       -d <DIGEST_MAP> \
       -q MIF_GM12878_CHC_REP1/MIF_GM12878_CHC_REP1.truncated_R1.fastq.gz \
       -r MIF_GM12878_CHC_REP1/MIF_GM12878_CHC_REP1.truncated_R2.fastq.gz \
       -o MIF_GM12878_CHC_REP1 \
       -x MIF_GM12878_CHC_REP1 \
       -j

All result files from this step are written to the same directory (``-o``)
and have the same prefix (``-x``) as the truncated reads.
The main result from this step is a BAM file with valid mapped read pairs that have not been classified as artifacts.
If the ``-j`` switch is used, then an additional BAM file is created containing all read pairs
that were determined to be invalid and therefore rejected.

Counting supporting read pairs for interacting digest pairs
===========================================================

In Diachromatic, an interaction is defined as any pair of digests having at least one supporting valid mapped read pair.
By using ``Diachromatic`` with the subcommand ``count``,
all interactions with their counts of supporting read pairs can be determined.
The corresponding digest map (``-d``) and
the BAM file from the previous step containing the valid mapped read pairs (``-v``) must be specified.
If the ``-s`` switch is used, the counts of supporting read pairs are reported separately by relative orientation.

.. code-block:: console

    $ java -Xmx32000m -jar Diachromatic.jar count \
       -d <DIGEST_MAP>  \
       -v MIF_GM12878_CHC_REP1/MIF_GM12878_CHC_REP1.valid_pairs.aligned.bam \
       -s \
       -o MIF_GM12878_CHC_REP1 \
       -x MIF_GM12878_CHC_REP1

The interactions are written to the following file:

.. code-block:: console

    MIF_GM12878_CHC_REP1/MIF_GM12878_CHC_REP1.interaction.counts.table.tsv

This file is in the Diachromatic interaction format:

.. code-block:: console

    chr1    46297999   46305684   E   chr1    51777391   51781717   N   2:0:1:0
    chr17   72411026   72411616   N   chr17   72712662   72724357   N   3:0:2:0
    chr7    69513952   69514636   N   chr7    87057837   87061499   E   4:0:3:0
    chr11    9641153    9642657   N   chr11   47259263   47272706   E   5:0:4:0

Each line represents an interaction.
Columns 1 to 3 and 5 to 7 contain the coordinates of the digest pair,
whereby the digest with the smaller coordinates always comes before the other digest.
Columns 4 and 8 indicate the enrichment states of the two digests.
An ``E`` means that the corresponding digest has been selected for target enrichment
and an ``N`` means that it has not been selected.
The last column contains the counts of the supporting read pairs
separated by relative orientations of mapped read pairs (``<Class 0>``:``<Class 1>``:``<Class 2>``:``<Class 3>``).

+--------+--------------------------------------+
| Class  | Relative orientation                 |
+========+============================+=========+
| ``0``  | Reads point inwards        |``-><-`` |
+--------+----------------------------+---------+
| ``1``  | Reads point outwards       |``<-->`` |
+--------+----------------------------+---------+
| ``2``  | Reads both point towards 3'| ``->->``|
+--------+----------------------------+---------+
| ``3``  | Reads both point towards 5'| ``<-<-``|
+--------+----------------------------+---------+

Filtering for cis-chromosomal long range interactions
=====================================================

Interactions between different chromosomes are referred to as trans-chromosomal
and interactions within the same chromosome as cis-chromosomal.
In this tutorial, we restrict our analysis to cis-chromosomal interactions.
Typically, interactions with particularly short distances are excluded from downstream analysis.
We define the distance between the two inner ends of interacting digests (column 3 and 6) as interaction distance
and discard all interactions with a distance smaller than ``20,000 bp``.
We also discard all interactions on chromosome ``chrM``.

.. code-block:: console

    $ mkdir gzdir
    $ awk '{if($1==$5 && $6-$3>=20000){print $0}}' MIF_GM12878_CHC_REP1/MIF_GM12878_CHC_REP1.interaction.counts.table.tsv \
       | grep -v chrM \
       | gzip > gzdir/MIF_GM12878_CHC_REP1.interaction.counts.table.clr_200000.tsv.gz

Do the last four steps for the other two replicates as well.
After that, the directory ``gzdir`` should contain three files.

.. code-block:: console

    $ ls gzdir
    MIF_GM12878_CHC_REP1.interaction.counts.table.clr_200000.tsv.gz
    MIF_GM12878_CHC_REP2.interaction.counts.table.clr_200000.tsv.gz
    MIF_GM12878_CHC_REP3.interaction.counts.table.clr_200000.tsv.gz

*************************************************************
``pooler.py``: Pooling interactions from different replicates
*************************************************************

To pool interactions from different replicates,
we discard those that occur in fewer than a specified number of replicates
and for the remaining, overlapping interactions add up the read pair counts separately by orientation.
For example, if the same interaction occurs in two replicates and has counts ``1:2:3:4`` for the one replicate
and counts ``4:3:2:1`` for the other, then the pooled counts will be ``5:5:5:5``.
We implemented this way of pooling in the pooler.py script,
the usage of which is demonstrated here:
`jupyter_notebooks/usage/usage_of_pooler.ipynb <https://github.com/TheJacksonLaboratory/diachrscripts/blob/master/jupyter_notebooks/usage/usage_of_pooler.ipynb>`_.

As input, the script expects a path to a directory containing gzipped files
in Diachromatic's interaction format (``--interaction-files-path``).
In addition, the minimum number of replicates in which a pooled interaction must occur
must be specified (``--required-replicates``).
In this tutorial, we require that an interaction occurs in at least two replicates.
The name of each file created has the same prefix (``--out-prefix``),
which can also contain a path to a pre-existing directory.

.. code-block:: console

    $ mkdir MIF_GM12878_CHC_REPC
    $ diachrscripts/pooler.py \
       --interaction-files-path gzdir \
       --required-replicates 2 \
       --out-prefix MIF_GM12878_CHC_REPC/MIF_GM12878_CHC_REPC

This command will generate the following two files:

.. code-block:: console

    $ ls MIF_GM12878_CHC_REPC | cat
    MIF_GM12878_CHC_REPC_at_least_in_2_replicates_summary.txt
    MIF_GM12878_CHC_REPC_at_least_in_2_replicates_interactions.tsv.gz

The first file contains summary statistics and the second file contains the pooled interactions.

**Note:** Diachromatic reports each interaction with at least one supporting read pair.
To pool interactions, all interactions must first be read into the main memory.
Depending on the size of the input files, this can lead to high memory requirements.

********************************************
``UICer.py``: Unbalanced Interaction Caller
********************************************

We implemented the calling of unbalanced interactions in the Python script ``UICer.py``.
The script performs the following processing steps:

1. **Randomization:** If no classification threshold is specified,
then a threshold is determined using a randomization procedure so that the FDR remains below 5%.

2. **Classification of interactions:** All interactions that do not have enough read pairs
to be classified as unbalanced at the chosen classification threshold are discarded.
The remaining interactions are classified as unbalanced or balanced and assigned a score.

3. **Selection of comparison sets:** From the unbalanced and balanced interactions,
two comparison sets are selected that are as large as possible and comparable
with respect to their total read pair counts per interaction.

The usage of ``UICer.py`` is demonstrated here:
`jupyter_notebooks/usage/usage_of_UICer.ipynb <https://github.com/TheJacksonLaboratory/diachrscripts/blob/master/jupyter_notebooks/usage/usage_of_UICer.ipynb>`__.

In this tutorial, we will use UICer with an FDR threshold of 5% (``--fdr-threshold``).
This will invoke the randomization procedure with 1,000 iterations (``--iter-num``).
Since the randomization procedure is very computationally intensive,
the iterations are carried out in four batches of equal size (``--thread-num``) with 250 iterations each.
As input, we'll use the pooled interactions file from the previous step (``--diachromatic-interaction-file``).

.. code-block:: console

    $ diachrscripts/UICer.py \
        --out-prefix MIF_GM12878_CHC_REPC/MIF_GM12878_CHC_REPC \
        --description-tag MIF_GM12878_CHC_REPC \
        --diachromatic-interaction-file MIF_GM12878_CHC_REPC/MIF_GM12878_CHC_REPC_at_least_in_2_replicates_interactions.tsv.gz \
        --fdr-threshold 0.05 \
        --iter-num 1000 \
        --random-seed 1 \
        --thread-num 4

``UICer.py`` reports summary statistics on all processing steps and generates corresponding plots.

.. code-block:: console

    $ ls MIF_GM12878_CHC_REPC | cat
    MIF_GM12878_CHC_REPC_at_least_in_2_replicates_summary.txt
    MIF_GM12878_CHC_REPC_at_least_in_2_replicates_interactions.tsv.gz
    MIF_GM12878_CHC_REPC_evaluated_and_categorized_interactions.tsv.gz
    MIF_GM12878_CHC_REPC_randomization_histogram_at_001.pdf
    MIF_GM12878_CHC_REPC_randomization_histogram_at_005.pdf
    MIF_GM12878_CHC_REPC_randomization_histogram_at_010.pdf
    MIF_GM12878_CHC_REPC_randomization_histogram_at_threshold.pdf
    MIF_GM12878_CHC_REPC_randomization_plot.pdf
    MIF_GM12878_CHC_REPC_randomization_table.txt
    MIF_GM12878_CHC_REPC_reports.txt

The main result is a file of classified interactions evaluated in terms of their imbalances in the four counts.
The format of this file corresponds to the Diachromatic interaction format with two additional columns,
one for scores to assess the imbalances in the four read pair counts and one for the interaction categories.
Here is one line for each of the categories for illustration:

.. code-block:: console

    chr1   245051445   245057234   N   chr1   245133022   245136428   E   16:0:0:6   6.62   DIX
    chr21   18333585    18336116   N   chr21   18782489    18791793   E   4:0:0:3    2.11   DI
    chrX   151978880   151979018   N   chrX   152449365   152452950   E   11:3:7:7   1.03   UIR
    chr1    31956115    31963217   N   chr1    32695361    32706402   E   1:2:2:2    0.30   UI

There are four interaction categories:

+-----------+--------------------------------------------------------------+
| Category  | Meaning                                                      |
+===========+==============================================================+
| ``DIX``   | Unbalanced interaction without reference interaction         |
+-----------+--------------------------------------------------------------+
| ``DI``    | Unbalanced interaction with reference interaction            |
+-----------+--------------------------------------------------------------+
| ``UIR``   | Balanced interaction selected as reference interaction       |
+-----------+--------------------------------------------------------------+
| ``UI``    | Balanced interaction not selected as reference interaction   |
+-----------+--------------------------------------------------------------+

**Note:** Depending on the size of the input and the number of iterations,
the randomization procedure can be very computationally intensive.
We provide the ``UICer.py`` output file resulting from the steps in this tutorial for download.

********************************************************
Analysis of unbalanced interactions in Jupyter Notebooks
********************************************************

We implemented a number of analysis and visualization methods in Python modules that we use in Jupyter notebooks
to investigate interactions regarding imbalances in their four read pair counts.

Get started
===========

All input files required for these analyses, including the file generated with ``UICer.py`` in this tutorial,
can be downloaded using this notebook: `jupyter_notebooks/Get_started.ipynb <https://github.com/TheJacksonLaboratory/diachrscripts/blob/master/jupyter_notebooks/Get_started.ipynb>`_

This notebook also contains brief explanations and links to the various analyses.

Analyses
========

There are notebooks for the following analyses:

1. `Frequencies of interaction configurations <https://github.com/TheJacksonLaboratory/diachrscripts/blob/master/jupyter_notebooks/analysis/frequencies_of_interaction_configurations.ipynb>`_

2. `Visualization of configurations at baited fragments <https://github.com/TheJacksonLaboratory/diachrscripts/blob/master/jupyter_notebooks/analysis/visualization_of_configurations.ipynb>`_

3. `Classification of baited fragments <https://github.com/TheJacksonLaboratory/diachrscripts/blob/master/jupyter_notebooks/analysis/baited_fragment_classification.ipynb>`_

4. `Bait analysis <https://github.com/TheJacksonLaboratory/diachrscripts/blob/master/jupyter_notebooks/analysis/bait_analysis.ipynb>`__

5. `Restriction fragment lengths <https://github.com/TheJacksonLaboratory/diachrscripts/blob/master/jupyter_notebooks/analysis/restriction_fragment_lengths.ipynb>`_

6. `Distance-dependent contact frequencies <https://github.com/TheJacksonLaboratory/diachrscripts/blob/master/jupyter_notebooks/analysis/distance_dependent_contact_frequencies.ipynb>`_
