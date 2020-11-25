.. _RST_Interaction_calling:

###################
Interaction calling
###################


We used the Java program Diachromatic to derive interactions from
unprocessed capture Hi-C reads.
There is already a
`publication <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6678864/>`_
and a
`detailed documentation <https://diachromatic.readthedocs.io/en/latest/index.html>`_
on this program.
This section is intended to give a brief overview
and to describe the details of the analysis of the
`capture Hi-C dataset for 17 hematopoietic cell types <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5123897/>`_.

The input of Diachromatic essentially consists of the following:

1. Paired-end read data in FASTQ format
2. A bowtie2 index for the reference sequence to which the reads will be mapped
3. The restriction enzyme that was used to generate the data
4. A file that contains the restriction fragments (or digests) of the entire genome that result from digestion with the restriction enzyme

Diachromatic processes the data in three steps:

1. Truncation of reads
2. Mapping and removal of artifact read pairs
3. Counting of read pairs that map to the same digest pairs

In the following subsection you will be guided through
the analysis with Diachromatic so that every step can be followed exactly.

**********
Input data
**********

Paired-end data in FASTQ format
===============================

For paired-end data, the reads of the read pairs are typically in two different files,
one for the forward reads (``R1``) and one for the reverse reads (``R2``).
For the data on the hematopoietic cell types,
the downloadable reads for given replicates are divided into multiple FASTQ files
(beyond the division into ``R1`` and ``R2``).
For each replicate, we have merged the files by writing their contents
into two new files, one for ``R1`` and one for ``R2``.
We made sure that the files for ``R1`` and ``R2`` are in sync,
i.e. that two reads with the same row index form a pair.
For example, we have the following file pairs for replicate 1 on megakaryocytes:

.. code-block:: console

    EGAF00001274989.fq, EGAF00001274990.fq
    EGAF00001274991.fq, EGAF00001274992.fq
    EGAF00001274993.fq, EGAF00001274994.fq

In each line, the first file belongs to ``R1`` and the second file to ``R2``.
We wrote these files into two new files as follows:

.. code-block:: console

    cat EGAF00001274989.fq EGAF00001274991.fq EGAF00001274993.fq > MK/JAV_MK_R10_1.fastq
    cat EGAF00001274990.fq EGAF00001274992.fq EGAF00001274994.fq > MK/JAV_MK_R10_2.fastq

In table TODO.xls,
we have documented which original files were concatenated in which order.


bowtie2 index for reference sequence
====================================

Diachromatic uses bowtie2 to map reads to a reference sequence,
which requires an index for the reference sequence.
Such an index can be created with bowtie2 or a pre-calculated index can be downloaded.
For the data on the hematopoietic cell types,
we used the reference sequence hg38 and downloaded a pre-calculated index from here:

.. code-block:: console

    ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz

Note that the index consists of several files:

.. code-block:: console

    hg38.1.bt2
    hg38.2.bt2
    hg38.3.bt2
    hg38.4.bt2
    hg38.rev.1.bt2
    hg38.rev.2.bt2

Only the path together with the prefix (``hg38``) of these files are passed to Diachromatic.


Restriction enzyme an digest map
================================

For the dataset on the 17 hematopoietic cell types,
the enzyme *HindIII* with the recognition sequence ``AAGCTT`` was used.
Diachromatic reports individual interactions as a pair of restriction fragments (or digest pairs)
along with the number of read pairs that map to this digest pair.
This requires the coordinates of all digests in the reference genome
that are passed to Diachromatic in form of a text file,
which we refer to as digest map.
For a given restriction enzyme and a reference genome,
a corresponding digest map can be created with the GOPHER software
`as described in the documentation <https://diachromatic.readthedocs.io/en/latest/digest.html>`_.
These are the first few line of such a digest map:

.. code-block:: console

    Chromosome      Fragment_Start_Position Fragment_End_Position   Fragment_Number 5'_Restriction_Site     3'_Restriction_Site     Length  5'_GC_Content   3'_GC_Content   5'_Repeat_Content       3'_Repeat_Content       Selected        5'_Probes       3'_Probes
    chr1    1       16007   1       None    HindIII 16007   0.000   0.000   0.000   0.003   F       0       0
    chr1    16008   24571   2       HindIII HindIII 8564    0.018   0.018   0.000   0.015   F       0       0
    chr1    24572   27981   3       HindIII HindIII 3410    0.046   0.046   0.000   0.044   F       0       0
    chr1    27982   30429   4       HindIII HindIII 2448    0.035   0.035   0.047   0.043   F       0       0

To ensure consistency,
we recommend creating the digest map from the same FASTA file that was used
to create the bowtie2 index.


*******************
Truncation of reads
*******************

Use Diachromatic to truncate the read pairs given in FASTQ format as follows:

.. code-block:: console

    java -jar Diachromatic.jar truncate \
       -e HindIII \
       -q MK/JAV_MK_R10_R1.fastq.gz \
       -r MK/JAV_MK_R10_R2.fastq.gz \
       -o MK \
       -x JAV_MK_R10


Diachromatic has an internal list of common restriction enzymes
and will use the appropriate recognition sequence and cutting positions
for ``-e HindIII``.
We use the previously downloaded and concatenated FASTQ files
for the forward (R1, ``-q``) and reverse (R2, ``-r``) as input.
An already existing directory for the output (``-o``) and a prefix
for all generated files (``-x``) can also be specified.
For capture Hi-C data, we don't use the ``--sticky-ends`` option,
i.e. we assume that the sticky ends resulting from the restriction
have been filled in.
More details on the truncation of reads can be found in the
`relevant section of the Diachromatic documentation <https://diachromatic.readthedocs.io/en/latest/truncate.html>`_.


******************************************
Mapping and removal of artifact read pairs
******************************************

Use Diachromatic to map the the truncated read pairs to the reference sequence as follows:

.. code-block:: console

    java -jar Diachromatic.jar align \
       -bsu \
       -d <DIGEST_MAP> \
       -q MK/JAV_MK_R10.truncated_R1.fastq.gz \
       -r MK/JAV_MK_R10.truncated_R2.fastq.gz \
       -b <BOWTIE2_EXECUTABLE> \
       -i <BOWTIE2_INDEX>/hg38 \
       -p 32 \
       -j \
       -o MK \
       -x JAV_MK_R10

In addition to mapping, Diachromatic removes duplicated read pairs and
keeps track of the number of read pairs for different duplication levels.
Depending on the size of the input and the actual duplication rate,
this can take up a lot of memory.
We therefore recommend having 16 to 32 GB memory available.

We use the more stringent mode of Diachromatic to define uniquely mapped reads,
i.e. reads that map to only one location (``-bsu``).
In order to determine artifact read pairs,
for example pairs mapped to the same digest,
the previously created digest map is required (``-d``).
We map the truncated reads from the previous step (``-q,-r``) to ``hg38``.

Diachromatic uses bowtie2 to map the reads to the reference genome.
To do this,
an executable bowtie2 file and an index for the reference must be specified (``-b``, ``-i``).
We use 32 threads for the maapping with bowtie2 (``-p``).

For possible subsequent investigation,
we write the rejected artifact read pairs to an extra BAM file (``-j``).
The valid read pairs are always written to a BAM file
with the suffix ``.valid_pairs.aligned.bam``.
We note that these files do not contain any read pairs that have
been mapped to non-canonical chromosomes
(e.g. ``chrUn_GL000216v2``).
The reads of a pair are mapped independently to all chromosomes,
but a pair for which at least one read is mapped to a non-canonical
chromosome cannot be re-paired.
This is the relevant section in the
`Diachromatic source code <https://github.com/TheJacksonLaboratory/diachromatic/blob/master/src/main/java/org/jax/diachromatic/align/ReadPair.java>`_.

.. code-block:: java

    // check if both reads are not on random chromosomes or EBV for hg38
    if (R1.getReferenceName().contains("_") || R2.getReferenceName().contains("_") || R1.getReferenceName().contains("EBV")|| R2.getReferenceName().contains("EBV")) {
        this.isPaired = false;
    }

The output can be redirected and given prefixes as with the ``truncate`` command.
More details on the mapping and removal of artifact read pairs can be found in the
`relevant section of the Diachromatic documentation <https://diachromatic.readthedocs.io/en/latest/mapping.html>`_.

***************************************************
Counting of valid read pairs mapped to digest pairs
***************************************************

Use Diachromatic to count valid read pairs between interacting digest pairs as follows:

.. code-block:: console

    java -jar Diachromatic.jar count \
       -d <DIGEST_MAP>  \
       -v JAV_MK_R10.valid_pairs.aligned.bam \
       -s \
       -o MK \
       -x JAV_MK_R10

In Diachromatic, interactions are defined as digest pairs that have at least
one supporting read pair.
In this step, the supporting read pairs for individual interactions are counted.
To do this, the digest map is required (``-d``).
We use the unique valid pairs from the previous step as input (``-v``),
i.e. duplicates and artifact read pairs have been removed.
We use the ``-s`` option so that the simple and twisted read pairs counts
of individual interactions are reported separately.
Note that the ``-s`` option is currently only available on the ``develop`` branch
of the GitHub repository for Diachromatic.

More details on counting read pairs between interacting digest regions can be found in the
`relevant section of the Diachromatic documentation <https://diachromatic.readthedocs.io/en/latest/count.html>`_.

The interactions with their read pair counts are written to the following file:

.. code-block:: console

    MK/JAV_MK_R10.interaction.counts.table.tsv


These are the first few lines from such a file:

.. code-block:: console

    chr1    46297999   46305684   A   chr1    51777391   51781717   I   2:1
    chr17   72411026   72411616   I   chr17   72712662   72724357   I   3:2
    chr7    69513952   69514636   I   chr7    87057837   87061499   A   4:3
    chr11    9641153    9642657   I   chr11   47259263   47272706   A   5:4

Each line represents one interaction.
Columns 1 to 3 and 5 to 7 contain the coordinates of the digest pair,
whereby the smaller coordinates are always in columns 1 to 3.

In column 4 and 8 there is either an ``A`` or an ``I``,
where column 4 belongs to the first and column 8 belongs to the second digest.
An ``A`` means that the corresponding digest was selected for target enrichment
and an ``I`` means that it was not selected.
The information about digests that were selected for enrichment
is taken from the digest map that was generated with GOPHER.
When generating the digest map with GOPHER,
we used the shortcut option *All protein coding genes*
because the analyzed data comes from a whole-promoter capture Hi-C
experiment.
This has the effect that all digests which contain at least on TSS
and are potentially suited for enrichment are marked with an ``A``
and all others with an ``I``.
The digests marked with an ``A`` do not exactly correspond to the digests
that were actually selected for the experiment.
This is due to different annotations as well as different criteria
for the selection of enrichable digests.
Therefore, we only used the markings with ``A`` and ``I`` at the beginning
to roughly asses the how well the enrichment worked.
For the following,
we did not use these markings with ``A`` and ``I``,
but instead a list of enriched digests from the original publication
of Javierre et al. 2016 (see below).
To make it easier to distinguish,
we denote enriched enriched digests with an ``E`` and non-enriched
digests with an ``N``, when using the annotation from this list.

The last columnn in a Diachromatic interaction file
shows the counts of simple and twisted read pairs
separated by a colon.
For example, ``5:4`` means that five simple and four twisted
read pairs were counted for an interaction.

************************************
Subsequent filtering of interactions
************************************

We filtered out interactions between different chromosomes (trans).
From the remaining interactions (cis),
we have filtered out short interactions with a distance 20,000 bp
and interactions with and on chromosome ``chrM``.
We implemented the filtering with the
`command line tool AWK <https://en.wikipedia.org/wiki/AWK>`_:

.. code-block:: console

    awk '{if($1==$5 && $6-$3>=20000){print $0}}' MK/JAV_MK_R10.interaction.counts.table.tsv \
       | grep -v chrM \
       | gzip > MK/gzdir/JAV_MK_R10.interaction.counts.table.clr_200000.tsv.gz

The filtered interactions are written to a gzip-compressed file
in a directory ``MK/gzdir``.
This is where the interactions for all replicates are written
and from there they are read in by the script ``01_combine_interactions_from_replicates.py``
(see :ref:`RST_Combining_interactions`).
