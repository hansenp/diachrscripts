.. _RST_Diachromatic_input_preparation:

##############################
Diachromatic input preparation
##############################

*************
bowtie2 index
*************

Diachromatic uses bowtie2 to map reads to a reference sequence,
which requires an index for the reference sequence.
Such an index can be created with bowtie2 or a pre-calculated index can e.g. be downloaded from the
`iGenomes website <https://support.illumina.com/sequencing/sequencing_software/igenome.html>`_

.. code-block:: console

    $ wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
    $ mkdir Homo_sapiens_UCSC_hg38
    $ tar -xf Homo_sapiens_UCSC_hg38.tar.gz -C Homo_sapiens_UCSC_hg38

The index consists of several files:

.. code-block:: console

    $ ls Homo_sapiens_UCSC_hg38/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/ | cat
    genome.1.bt2
    genome.2.bt2
    genome.3.bt2
    genome.4.bt2
    genome.fa
    genome.rev.1.bt2
    genome.rev.2.bt2

Only the path together with the file prefix of these files are passed to Diachromatic.

.. code-block:: console

    -i <PATH>/genome

In addition to the bowtie2 index, the downloaded directory contains other data and is about ``30 GB`` in size.

***********
Digest file
***********

Diachromatic reports individual interactions as a pairs of restriction fragments (or digest pairs)
along with the number of supporting read pairs.
This requires the coordinates of all possible digests in the entire reference genome.
These are passed to Diachromatic in form of a text file, which we refer to as digest file.
For a given restriction enzyme and a reference genome,
a corresponding digest map can be created with the GOPHER software
`as described in the documentation <https://diachromatic.readthedocs.io/en/latest/digest.html>`__.
For example, here are the first few lines from a digest file:

.. code-block:: console

    Chromosome      Fragment_Start_Position Fragment_End_Position   Fragment_Number 5'_Restriction_Site     3'_Restriction_Site     Length  5'_GC_Content   3'_GC_Content   5'_Repeat_Content       3'_Repeat_Content       Selected        5'_Probes       3'_Probes
    chr1    1       16007   1       None    HindIII 16007   0.000   0.000   0.000   0.003   F       0       0
    chr1    16008   24571   2       HindIII HindIII 8564    0.018   0.018   0.000   0.015   F       0       0
    chr1    24572   27981   3       HindIII HindIII 3410    0.046   0.046   0.000   0.044   F       0       0
    chr1    27982   30429   4       HindIII HindIII 2448    0.035   0.035   0.047   0.043   F       0       0

Each line represents one digest.
The first three columns contain the digest coordinates.
In the ``Selected`` column, enriched digests are marked with a ``T`` and all others with an ``F``.
Diachromatic passes the information about enriched digests through to the reported interactions.
In Diachromatic interaction files,
an ``E`` corresponds to an ``T`` and an ``N`` to an ``F``.

.. code-block:: console

    chr11    9641153    9642657   E   chr11   47259263   47272706   N   5:4


Selecting enriched digests
==========================

If you have created the capture Hi-C design for your own experiment with GOPHER,
then all digests for which baits were ordered will be marked with a ``T``.
If this is not the case,
then with promoter capture Hi-C experiments for which all promoters have been selected for target enrichment,
there is the possibility to create a design with the shortcut ``All protein coding genes``
before exporting the digest file from GOPHER.
However, this is inaccurate because the selection of digests depends on the software used and the parameter settings.
For cases in which the coordinates of enriched digests are known,
we provide a script that can be used to create an appropriate digest file.

.. code-block:: console

    $ python diachrscripts/additional_scripts/create_diachromatic_digest_file.py
        --enriched-digests-file <ENRICHED_DIGEST_COORDINATES>.bed
        --diachromatic-digest-file <DIGEST_FILE_TEMPLATE>.txt
        --out-prefix <CUSTOM_DIGEST_FILE_PREFIX>

This script is passed a file containing the coordinates of the digests that have actually been selected for enrichment.
Such information can be found, for example, in the supplementary material of the corresponding publication
(see examples below).
In addition, a digest file for the appropriate reference genome and restriction enzyme must be passed.
You can export a digest file from GOPHER before creating a design or you can  take any existing digest file.
It is only used as a template and the information related to enrichment will be completely rewritten.

Selecting enriched digests for Mifsud et al. 2015
=================================================

Supplementary Table 4 of the work published by
`Mifsud et al. 2015 <https://pubmed.ncbi.nlm.nih.gov/25938943/>`__
contains the coordinates and sequences of the baits used.
Save this table in text format and extract the coordinates.

.. code-block:: console

    $ cat mifsud_supplementary_table_4.txt | \
        awk '{if($1 ~ /^>/){split($0,a," ");split(a[1],b,":");gsub(/>C/,"c",b[1]);split(b[2],c,"-");print b[1]"\t"c[1]"\t"c[2]}}' \
        > mifsud_bait_coords_hg19.bed

Use
`UCSC's LiftOver tool <https://genome.ucsc.edu/cgi-bin/hgLiftOver>`_
to convert the coordinates to from ``hg19`` to ``hg38`` and save the resulting file as ``mifsud_bait_coords_hg38.bed``.
Generate a digest file for ``hg38`` with GOPHER and extract the coordinates of all digests in the genome.

.. code-block:: console

    $ tail -n+2 template_digest_file_hg38_HindIII.txt \
    | awk '{print $1"\t"$2"\t"$3}' > all_hg38_digests.bed

Use
`bedtools <https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html>`_
to extract all digests that contain at least one bait completely.

.. code-block:: console

    $ intersectBed -wa -u -F 1.00 -a all_hg38_digests.bed -b mifsud_bait_coords_hg38.bed \
    > mifsud_baited_digests_hg38.bed

Use our script to create a digest file in which digests that Mifsud et al. have selected for enrichment are marked with
a ``T`` and all others with an ``F``.

.. code-block:: console

    $ python diachrscripts/additional_scripts/create_diachromatic_digest_file.py \
        --enriched-digests-file mifsud_baited_digests_hg38.bed \
        --diachromatic-digest-file template_digest_file_hg38_HindIII.txt \
        --out-prefix mifsud_hg38_HindIII

This will produce the file ``mifsud_hg38_HindIII_diachromatic_digest_file.txt`` that can be used as input for Diachromatic.

Selecting enriched digests for Javierre et al. 2016
===================================================

For the work published by
`Javierre et al. 2016 <https://pubmed.ncbi.nlm.nih.gov/27863249/>`__,
the ``hg19`` coordinates of the baited digests can be downloaded from
`OFS <https://osf.io/e594p/>`__.

Download an archive that expands into a *design folder* that can be provided to the interaction caller
`CHiCAGO <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4908757/>`_.

.. code-block:: console

    $ wget -O human_PCHiC_hg19_HindIII_design.tar.gz https://osf.io/e594p/download
    $ tar -xf human_PCHiC_hg19_HindIII_design.tar.gz

Along with other files, this folder contains the
`CHiCAGO's bait map file <https://bioconductor.org/packages/devel/bioc/vignettes/Chicago/inst/doc/Chicago.html>`_
that consists of the following columns:
``chr``, ``start``, ``end``, ``fragmentID``, ``geneName``.

.. code-block:: console

    $ head -n 4 Human_hg19/Digest_Human_HindIII_baits_e75_ID.baitmap
        1	831895	848168	218	RP11-54O7.16;RP11-54O7.1
        1	848169	850618	219	RP11-54O7.2
        1	850619	874081	220	AL645608.1;RP11-54O7.3;SAMD11
        1	889424	903640	223	KLHL17;NOC2L;PLEKHN1

Convert the bait map file into BED format.

.. code-block:: console

    $ awk '{print "chr"$1"\t"$2"\t"$3}' Human_hg19/Digest_Human_HindIII_baits_e75_ID.baitmap \
    > mifsud_baited_digests_hg19.bed

Use
`UCSC's LiftOver tool <https://genome.ucsc.edu/cgi-bin/hgLiftOver>`_
to convert the coordinates to from ``hg19`` to ``hg38`` and save the resulting file as
``javierre_baited_digests_hg38.bed``.

Use our script to create a digest file in which digests that Javierre et al. have selected for enrichment are marked
with a ``T`` and all others with an ``F``.

.. code-block:: console

    $ python diachrscripts/additional_scripts/create_diachromatic_digest_file.py \
        --enriched-digests-file javierre_baited_digests_hg38.bed \
        --diachromatic-digest-file template_digest_file_hg38_HindIII.txt \
        --out-prefix javierre_hg38_HindIII

This will produce the file ``javierre_hg38_HindIII_diachromatic_digest_file.txt`` that can be used as input for Diachromatic.

XXXXXXXXXXX.

Coordinates are avaiable for a total of 22,076 digests.

These coordinates refer to the genome build ``hg19``.
We used
`UCSC's LiftOver tool <https://genome.ucsc.edu/cgi-bin/hgLiftOver>`_
to convert the coordinates to ``hg38``.
22,056 digests were successfully converted  to ``hg38``.
The conversion failed for 20 digests
because ``hg19`` coordinates in ``hg38``
are either split or partially deleted.
The resulting file in BED format can be found here:

.. code-block:: console

    additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed

hg38 digest file
================

We wrote a Python script to overwrite the values in the ``Selected`` column
of a digest file:

.. code-block:: console

    $ python diachrscripts/additional_scripts/create_diachromatic_digest_file.py
       --enriched-digests-file Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed
       --diachromatic-digest-file no_digests_selected_HindIII_hg38_DigestedGenome.txt
       --out-prefix /JAV_hg38_HindIII


It is important that the coordinates in the two files refer to the same genome build,
e.g. ``hg19`` or ``hg38``.

For each line of the digest file, the script checks
whether there is a digest with matching coordinates in the BED file.
If this is the case, the ``Selected`` field is overwritten with a ``T`` and otherwise with an ``F``.
Furthermore, the fields ``5'_Probes`` and ``3'_Probes`` are set to ``1``.

We applied the script to the prepared enriched digest BED file for the Javierre data
and the digest file for ``hg38`` and ``HindIII`` in which no digest is marked as selected.
For the command above,
the created digest file has the following name:

.. code-block:: console

    JAV_HindIII_hg38_diachromatic_digest_file.txt

The script reports that for 22,008 of the 22,056 enriched digests
no matching coordinates were found in the digest file,
i.e. no matching coordinates were found for 48 digests.
The coordinates of these digests are written to the following file:

.. code-block:: console

    JAV_HindIII_hg38_digests_not_found.bed

The script has an option ``--verbose`` that can be used to examine such cases
more closely by printing the associated lines from the digest file.
Three categories of error were responsible for the 48 cases in which a digest could not be mapped.
In 34 cases, the enriched digest is shifted three positions to the right
with respect to the corresponding digest in the Diachromatic digest file.
In 10 cases, the enriched digest spans a restriction site
(i.e., overlaps two or more digests in the Diachromatic digest file).
And in four cases, the enriched digest is completely contained in a digest
from the Diachromatic digest file.
We assumed that these cases result from the LiftOver from ``hg19`` to ``hg38``
and repeated the same procedure for ``hg19``.
In this case, all enriched digests are found in the Diachromatic digest file,
which confirms our assumption.

