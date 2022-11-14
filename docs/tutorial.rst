.. _RST_tutorial:

########
Tutorial
########

This tutorial describes the entire workflow from downloading the data, calling the interactions with ``Diachromatic``,
pooling interactions from different replicates, calling unbalanced interactions with ``DICer`` up to the various
analyzes of interactions with unbalanced read counts that can be performed in Jupyter Notebooks.

.. image:: img/analysis_flowchart.png

***************************
Downloading paired-end data
***************************

Use the script
`dumpy.sh <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/additional_scripts/dumpy.sh>`__
to download promoter capture Hi-C data on GM12878 cells
`(Mifsud et al. 2015) <https://pubmed.ncbi.nlm.nih.gov/25938943/>`_.
For this dataset, ``200 GB`` hard disk space must be available.


.. code-block:: console

    $ diachrscripts/additional_scripts/dumpy.sh MIF_REP1 MIF_REP1 "ERR436029"
    $ diachrscripts/additional_scripts/dumpy.sh MIF_REP2 MIF_REP2 "ERR436028;ERR436030;ERR436033"
    $ diachrscripts/additional_scripts/dumpy.sh MIF_REP3 MIF_REP3 "ERR436026;ERR436031"

.. code-block:: console

    $ diachrscripts/additional_scripts/dumpy.sh MON_IPSC_REP1 MON_IPSC_REP1 "ERR2365275"
    $ diachrscripts/additional_scripts/dumpy.sh MON_IPSC_REP2 MON_IPSC_REP2 "ERR2365276"
    $ diachrscripts/additional_scripts/dumpy.sh MON_IPSC_REP3 MON_IPSC_REP3 "ERR2365277"

.. code-block:: console

    $ diachrscripts/additional_scripts/dumpy.sh MON_CM_REP1 MON_CM_REP1 "ERR2365269"
    $ diachrscripts/additional_scripts/dumpy.sh MON_CM_REP2 MON_CM_REP2 "ERR2365270"
    $ diachrscripts/additional_scripts/dumpy.sh MON_CM_REP3 MON_CM_REP3 "ERR2365271"

A separate directory is created for each of the three replicates.
Because it is paired-end data, there are pairs of forward and reverse FASTQ files with
suffixes ``_1`` and ``_2``.
Both files have the same number of lines and reads with the same line index form a pair.

.. code-block:: console

    $ ls -lh MIF_REP*/
    MIF_REP1/:
    total 28G
    12G MIF_REP1_1.fastq.gz
    12G MIF_REP1_2.fastq.gz
    242 MIF_REP1_md5.txt

    MIF_REP2/:
    total 87G
    36G MIF_REP2_1.fastq.gz
    37G MIF_REP2_2.fastq.gz
    490 MIF_REP2_md5.txt

    MIF_REP3/:
    total 57G
    24G MIF_REP3_1.fastq.gz
    24G MIF_REP3_2.fastq.gz
    366 MIF_REP3_md5.txt

The MD5 checksums are as follows:

.. code-block:: console

    $ cat MIF_REP*/MIF_REP*_md5.txt | grep -v ERR
    c68ab47d9f8cb9b49c1701ff04c7fd41  MIF_REP1/MIF_REP1_1.fastq.gz
    338980ac0a0bc1a3ea26bc6cf55dfb43  MIF_REP1/MIF_REP1_2.fastq.gz
    b6c248b48ded9eee8c4e951f4ea8d4a4  MIF_REP2/MIF_REP2_1.fastq.gz
    a3ff3dfde21d39dc41dc3c1212cbb3e3  MIF_REP2/MIF_REP2_2.fastq.gz
    5ceed4dd8784a467afa35f00b83d132d  MIF_REP3/MIF_REP3_1.fastq.gz
    4260c37ed24ee29c5046ef74096ba1f7  MIF_REP3/MIF_REP3_2.fastq.gz

******************************************
Calling interactions with ``Diachromatic``
******************************************

We used the Java program
`Diachromatic <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6678864/>`__
to derive interactions from Hi-C and capture Hi-C data.
There is a more
`detailed documentation <https://diachromatic.readthedocs.io/en/latest/index.html>`__
on this program.
Diachromatic processes the data in three steps:

1. Truncation of chimeric reads ``>`` Truncated reads
2. Mapping and artifact removal ``>`` Valid mapped read pairs
3. Counting supporting read pairs for interacting digest pairs ``>`` Interactions

Truncation of chimeric reads
============================

The sequencing of Hi-C libraries can result in chimeric reads containing sequences from different regions.
Such reads cannot be mapped.
Therefore, they must first be truncated so that they are no longer chimeric.
This can be done with ``Diachromatic`` using the subcommand ``truncate``.
The chimeric reads must be cut at the ligation sites, which is why the restriction enzyme used for the experiment must
be specified (``-e``).
The prepared FASTQ files with the forward and reverse reads are specified using the ``-q`` and ``-r`` options.

.. code-block:: console

    $ java -jar Diachromatic.jar truncate \
       -e HindIII \
       -q MIF_REP1/MIF_REP1_1.fastq.gz \
       -r MIF_REP1/MIF_REP1_2.fastq.gz \
       -o MIF_REP1 \
       -x MIF_REP1

All result files are written to the directory specified by the option ``-o`` and have the same prefix specified by the
option ``-x``.

Mapping and artifact removal
============================

For Hi-C data, no distribution particular of distances between reads of mapped pairs can be assumed (insert size).
However, for paired-end data, read mappers rely on a minimum and maximum insert size.
Therefore, the truncated forward and reverse reads must be mapped independently, like single-end data, and the mapped
reads must be re-paired afterwards.
In addition, there are certain rules by which artifacts that are specific to Hi-C data can be recognized and removed.
This can be done with ``Diachromatic`` using the subcommand ``align`` for which we recommend having ``16`` to ``32 GB``
memory available.
For the single-end mappings, paths to ``bowtie2`` (``-b``) and to an index for the matching reference sequence (``-i``)
must be specified. If the ``-bsu`` is used, then reads are considered to be mapped uniquely if they map to only one
location. The ``-p`` option specifies how many CPUs can be used by ``bowtie2``.
For the detection of artifacts, a digest file is required, which contains all restriction fragments resulting from a
complete digestion of the genome and must be specified via the option ``-d``.
The FASTQ files with the truncated forward and reverse reads are specified using the ``-q`` and ``-r`` options.

In order to execute the following command, the ``bowtie2`` index and the digest map must first be prepared.
How to do this is documented here: :ref:`RST_Diachromatic_input_preparation`.

.. code-block:: console

    $ java -Xmx32000m -jar Diachromatic.jar align \
       -b <BOWTIE2_EXECUTABLE> \
       -i <BOWTIE2_INDEX_PATH>/genome \
       -bsu \
       -p 4 \
       -d <DIGEST_MAP> \
       -q MIF_REP1/MIF_REP1.truncated_R1.fastq.gz \
       -r MIF_REP1/MIF_REP1.truncated_R2.fastq.gz \
       -o MIF_REP1 \
       -x MIF_REP1 \
       -j

All result files from this step are written to the same directory (``-o``) and have the same prefix (``-x``) as the
truncated reads.
The main result from this step is a BAM file with valid mapped read pairs that have not been classified as artifacts.
If the ``-j`` option is used, then an additional BAM file is created containing all read pairs that were determined to
be invalid and therefore rejected.

Counting supporting read pairs for interacting digest pairs
===========================================================

In ``Diachromatic``, an interactions is defined as any pair of digests having at least one supporting valid mapped read
pair. Using the subcommand ``count``, the number of supporting read pairs for all interactions can be determined.
To do this, a corresponding digest map (``-d``) and a BAM file containing valid mapped read pairs (``-v``) are required.
The ``-s`` option causes the read pair counts to be reported separately for the four types.

.. code-block:: console

    $ java -Xmx32000m -jar Diachromatic.jar count \
       -d <DIGEST_MAP>  \
       -v MIF_REP1/MIF_REP1.valid_pairs.aligned.bam \
       -s \
       -o MIF_REP1 \
       -x MIF_REP1

The interactions are written to the following file:

.. code-block:: console

    MIF_REP1/MIF_REP1.interaction.counts.table.tsv

This file is in Diachromatic's interaction format:

.. code-block:: console

    chr1    46297999   46305684   E   chr1    51777391   51781717   N   2:0:1:0
    chr17   72411026   72411616   N   chr17   72712662   72724357   N   3:0:2:0
    chr7    69513952   69514636   N   chr7    87057837   87061499   E   4:0:3:0
    chr11    9641153    9642657   N   chr11   47259263   47272706   E   5:0:4:0

Each line represents an interaction.
Columns 1 to 3 and 5 to 7 contain the coordinates of the digest pair,
whereby the digest with the smaller coordinates always comes before the other digest.
Columns 4 and 8 indicate the enrichment states of the digests.
An ``E`` means that the corresponding digest has been selected for target enrichment
and an ``N`` means that it has not been selected.
The last column contains the counts of the supporting read pairs separated by type
(``<Type 0>``:``<Type 1>``:``<Type 2>``:``<Type 3>``).

Filtering for cis-chromosomal long range interactions
=====================================================

Interactions between different chromosomes are referred to as trans-chromosomal and interactions within the same
chromosome cis-chromosomal.
We restricted our analyzes to cis-chromosomal interactions.
Typically, interactions with particularly short distances are excluded from downstream analyzes.
We define the distance between the two inner ends of interacting digests (column 3 and 6) as interaction distance
and discard all interactions with a distance smaller than 20,000 bp.
We also discard all interactions on chromosome ``chrM``.

.. code-block:: console

    $ mkdir gzdir
    $ awk '{if($1==$5 && $6-$3>=20000){print $0}}' MIF_REP1/MIF_REP1.interaction.counts.table.tsv \
       | grep -v chrM \
       | gzip > gzdir/MIF_REP1.interaction.counts.table.clr_200000.tsv.gz

Do the last four steps for the other two replicates as well.
After that, the directory ``gzdir`` should contain three files.

.. code-block:: console

    $ ls gzdir
    MIF_REP1.interaction.counts.table.clr_200000.tsv.gz
    MIF_REP2.interaction.counts.table.clr_200000.tsv.gz
    MIF_REP3.interaction.counts.table.clr_200000.tsv.gz

**********************************************
Pooling interactions from different replicates
**********************************************

This is described here: :ref:`RST_Interaction_pooling`.

.. code-block:: console

    $ mkdir MIF_REPC
    $ diachrscripts/additional_scripts/pooler.py \
       --interaction-files-path gzdir \
       --required-replicates 2 \
       --out-prefix MIF_REPC/MIF_REPC

.. code-block:: console

    $ ls MIF_REPC | cat
    MIF_REPC_at_least_in_2_replicates_summary.txt
    MIF_REPC_at_least_in_2_replicates_interactions.tsv.gz


**********************************************
Calling unbalanced interactions with ``DICer``
**********************************************

So far, this is only described in this
`Jupyter Notebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/Demonstration_of_DICer.ipynb>`__.

.. code-block:: console

    $ diachrscripts/DICer.py \
        --out-prefix MIF_REPC/MIF_REPC \
        --description-tag MIF_REPC \
        --diachromatic-interaction-file MIF_REPC/MIF_REPC_at_least_in_2_replicates_interactions.tsv.gz \
        --fdr-threshold 0.05 \
        --iter-num 1000 \
        --random-seed 1 \
        --thread-num 4

``DICer`` generates a file with the evaluated and categorized interactions and several files with statistics on the
various processing steps.

.. code-block:: console

    $ ls MIF_REPC | cat
    MIF_REPC_at_least_in_2_replicates_summary.txt
    MIF_REPC_at_least_in_2_replicates_interactions.tsv.gz
    MIF_REPC_evaluated_and_categorized_interactions.tsv.gz
    MIF_REPC_randomization_histogram_at_001.pdf
    MIF_REPC_randomization_histogram_at_005.pdf
    MIF_REPC_randomization_histogram_at_010.pdf
    MIF_REPC_randomization_histogram_at_threshold.pdf
    MIF_REPC_randomization_plot.pdf
    MIF_REPC_randomization_table.txt
    MIF_REPC_reports.txt

The format of the interaction file corresponds to the Diachromatic interaction format with two additional columns for
a score to evaluate the imbalances in the four counts and the interaction category.
Here is one line for each category as an example:

.. code-block:: console

    chr1   245051445   245057234   N   chr1   245133022   245136428   E   16:0:0:6   6.62   DIX
    chr21   18333585    18336116   N   chr21   18782489    18791793   E   4:0:0:3    2.11   DI
    chrX   151978880   151979018   N   chrX   152449365   152452950   E   11:3:7:7   1.03   UIR
    chr1    31956115    31963217   N   chr1    32695361    32706402   E   1:2:2:2    0.30   UI

The tags for the interaction categories have the following meanings:

+-----------+--------------------------------------------------------------+
| Category  | Meaning                                                      |
+===========+==============================================================+
| ``DIX``   | Unbalanced counts no reference interaction could be selected |
+-----------+--------------------------------------------------------------+
| ``DI``    | Unbalanced counts reference interaction could be selected    |
+-----------+--------------------------------------------------------------+
| ``UIR``   | Balanced counts selected as reference interaction            |
+-----------+--------------------------------------------------------------+
| ``UI``    | Balanced counts not selected as reference interaction        |
+-----------+--------------------------------------------------------------+

