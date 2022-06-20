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

    $ diachrscripts/additional_scripts/dumpy.sh MIF_R1 MIF_R1 "ERR436029"
    $ diachrscripts/additional_scripts/dumpy.sh MIF_R2 MIF_R2 "ERR436028;ERR436030;ERR436033"
    $ diachrscripts/additional_scripts/dumpy.sh MIF_R3 MIF_R3 "ERR436026;ERR436031"

A separate directory is created for each of the three replicates.
Because it is paired-end data, there are pairs of forward and reverse FASTQ files with
suffixes ``_1`` and ``_2``.
Both files have the same number of lines and reads with the same line index form a pair.

.. code-block:: console

    $ ls -lh MIF_R*/
    MIF_R1/:
    total 28G
    12G MIF_R1_1.fastq.gz
    12G MIF_R1_2.fastq.gz
    242 MIF_R1_md5.txt

    MIF_R2/:
    total 87G
    36G MIF_R2_1.fastq.gz
    37G MIF_R2_2.fastq.gz
    490 MIF_R2_md5.txt

    MIF_R3/:
    total 57G
    24G MIF_R3_1.fastq.gz
    24G MIF_R3_2.fastq.gz
    366 MIF_R3_md5.txt

The MD5 checksums are as follows:

.. code-block:: console

    $ cat MIF_R*/MIF_R*_md5.txt | grep -v ERR
    c68ab47d9f8cb9b49c1701ff04c7fd41  MIF_R1/MIF_R1_1.fastq.gz
    338980ac0a0bc1a3ea26bc6cf55dfb43  MIF_R1/MIF_R1_2.fastq.gz
    b6c248b48ded9eee8c4e951f4ea8d4a4  MIF_R2/MIF_R2_1.fastq.gz
    a3ff3dfde21d39dc41dc3c1212cbb3e3  MIF_R2/MIF_R2_2.fastq.gz
    5ceed4dd8784a467afa35f00b83d132d  MIF_R3/MIF_R3_1.fastq.gz
    4260c37ed24ee29c5046ef74096ba1f7  MIF_R3/MIF_R3_2.fastq.gz

******************************************
Calling interactions with ``Diachromatic``
******************************************

This is described here: :ref:`RST_Interaction_calling`.


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
       -q MIF_R1/MIF_R1_1.fastq.gz \
       -r MIF_R1/MIF_R1_2.fastq.gz \
       -o MIF_R1 \
       -x MIF_R1

All result files are written to the directory specified by the option ``-o`` and have the same prefix specified by the
option ``-x``.

Mapping and artifact removal
============================

For Hi-C data, no distribution particular of distances between reads of mapped pairs can be assumed (insert size).
However, for paired-end data, read mappers rely on a minimum and maximum insert size.
Therefore, the truncated forward and reverse reads must be mapped independently, like single-end data, and the reads
must be re-paired afterwards.
Re-pairing only requires that both reads of a pair have been mapped.
In addition, there are certain rules by which artifacts can be recognized that are specific to Hi-C data.
This can be done with ``Diachromatic`` using the subcommand ``align``.
For the single-end mappings, paths to ``bowtie2`` (``-b``) and to an index for the matching reference sequence (``-i``)
must be specified. If the ``-bsu`` is used, then reads are considered to be mapped uniquely if they map to only one
location. The ``-p`` option specifies how many CPUs can be used by ``bowtie2``.
For the detection of artifacts, a digest file is required, which must be specified via the option ``-d``.
The FASTQ files with the truncated forward and reverse reads are specified using the ``-q`` and ``-r`` options.

.. code-block:: console

    $ java -jar Diachromatic.jar align \
       -b <BOWTIE2_EXECUTABLE> \
       -i <BOWTIE2_INDEX>/hg38 \
       -bsu \
       -p 4 \
       -d <DIGEST_MAP> \
       -q MIF_R1/MIF_R1.truncated_R1.fastq.gz \
       -r MIF_R1/MIF_R1.truncated_R2.fastq.gz \
       -o MIF_R1 \
       -x MIF_R1 \
       -j

All result files from this step are written to the same directory (``-o``) and have the same prefix (``-x``) as the
truncated reads.
The main result from this step is a BAM file with valid read pairs that have not been classified as artifacts.
If the ``-j`` option is used, then an additional BAM file is created containing all read pairs that were determined to
be invalid and therefore rejected.

Counting valid read pairs
=========================

XXX.

.. code-block:: console

    $ java -jar Diachromatic.jar count \
       -d <DIGEST_MAP>  \
       -v MIF_R1/MIF_R1.valid_pairs.aligned.bam \
       -s \
       -o MIF_R1 \
       -x MIF_R1

XXX.

Filtering for cis-chromosomal long range interactions
=====================================================

XXX.

.. code-block:: console

    $ mkdir gzip
    $ awk '{if($1==$5 && $6-$3>=20000){print $0}}' MIF_R1/MIF_R1.interaction.counts.table.tsv \
       | grep -v chrM \
       | gzip > gzdir/MIF_R1.interaction.counts.table.clr_200000.tsv.gz

Do the last four steps for the other two replicates as well.
After that, the directory ``gzdir`` should contain three files.

.. code-block:: console

    $ ls gzdir
    MIF_R1.interaction.counts.table.clr_200000.tsv.gz
    MIF_R2.interaction.counts.table.clr_200000.tsv.gz
    MIF_R3.interaction.counts.table.clr_200000.tsv.gz

**********************************************
Pooling interactions from different replicates
**********************************************

This is described here: :ref:`RST_Combining_interactions`.

.. code-block:: console

    $ mkdir MIF_RALT
    $ diachrscripts/additional_scripts/pooler.py \
       --interaction-files-path gzdir \
       --required-replicates 2 \
       --out-prefix MIF_RALT/MIF_RALT

.. code-block:: console

    $ ls MIF_RALT | cat
    MIF_RALT_at_least_in_2_replicates_summary.txt
    MIF_RALT_at_least_in_2_replicates_interactions.tsv.gz


**********************************************
Calling unbalanced interactions with ``DICer``
**********************************************

So far, this is only described in this
`Jupyter Notebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/Demonstration_of_DICer.ipynb>`__.

.. code-block:: console

    $ diachrscripts/DICer.py \
        --out-prefix MIF_RALT/MIF_RALT \
        --description-tag MIF_RALT \
        --diachromatic-interaction-file MIF_RALT/MIF_RALT_at_least_in_2_replicates_interactions.tsv.gz \
        --fdr-threshold 0.05 \
        --iter-num 1000 \
        --random-seed 1 \
        --thread-num 4

``DICer``` generates a file with the evaluated and categorized interactions and several files with statistics on the
various processing steps.

.. code-block:: console

    $ ls MIF_RALT | cat
    MIF_RALT_at_least_in_2_replicates_summary.txt
    MIF_RALT_at_least_in_2_replicates_interactions.tsv.gz
    MIF_RALT_evaluated_and_categorized_interactions.tsv.gz
    MIF_RALT_randomization_histogram_at_001.pdf
    MIF_RALT_randomization_histogram_at_005.pdf
    MIF_RALT_randomization_histogram_at_010.pdf
    MIF_RALT_randomization_histogram_at_threshold.pdf
    MIF_RALT_randomization_plot.pdf
    MIF_RALT_randomization_table.txt
    MIF_RALT_reports.txt

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

******************************************************
Performing various analyzes on unbalanced interactions
******************************************************

We have implemented all analyzes following the calling of unbalanced interactions in different Jupyter Notebooks.
The ``DiachromaticInteractionSet`` is the central data structure in all of these analyzes.
It can be created from an interaction file generated with ``DICer``.

.. code-block:: python

    from diachr import DiachromaticInteractionSet
    d11_interaction_set = DiachromaticInteractionSet()
    d11_interaction_set.parse_file(
        i_file = "MIF_RALT/MIF_RALT_evaluated_and_categorized_interactions.tsv.gz",
        verbose = True)

Interaction distances
=====================

See this
`notebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/interaction_frequency_distance_analysis.ipynb>`__
and this
`one <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/interaction_frequency_distance_analysis_2.ipynb>`__.


Frequencies of read types and configurations of interactions
============================================================

See this
`notebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/read_pair_and_interaction_types.ipynb>`__.

Representation of interactions in triangle heatmaps
===================================================

See this
`notebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/dtvis.ipynb>`__.

Classification of baited digests
================================

See this
`notebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/interactions_at_baited_digests_select_baited_digests.ipynb>`__.

TAD boundaries
==============

See this
`notebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/tad_boundaries.ipynb>`__.


