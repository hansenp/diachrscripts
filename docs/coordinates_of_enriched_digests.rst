.. _RST_coordinates_of_enriched_digests:

###############################
Coordinates of enriched digests
###############################

When analyzing data from capture Hi-C experiments,
it is important to know which restriction fragments or digests were selected for enrichment.
In this section,
we describe how to introduce this information into the analysis with
``Diachromatic`` and ``diachrscripts``.

GOPHER's digest file as input for Diachromatic
==============================================

If the probes for a capture Hi-C experiment were created with GOPHER,
then the digest file, which
`can be exported from the corresponding GOPHER design <https://diachromatic.readthedocs.io/en/latest/digest.html>`_,
contains precise information about all enriched digests.
For example, here are the first few lines from a GOPHER digest file:

.. code-block:: console

    Chromosome      Fragment_Start_Position Fragment_End_Position   Fragment_Number 5'_Restriction_Site     3'_Restriction_Site     Length  5'_GC_Content   3'_GC_Content   5'_Repeat_Content       3'_Repeat_Content       Selected        5'_Probes       3'_Probes
    chr1    1       16007   1       None    HindIII 16007   0.000   0.000   0.000   0.003   F       0       0
    chr1    16008   24571   2       HindIII HindIII 8564    0.018   0.018   0.000   0.015   F       0       0
    chr1    24572   27981   3       HindIII HindIII 3410    0.046   0.046   0.000   0.044   F       0       0
    chr1    27982   30429   4       HindIII HindIII 2448    0.035   0.035   0.047   0.043   F       0       0

Each line represents one digest,
and the file contains the digests for the entire genome.
In the ``Selected`` column,
enriched digests are marked with a ``T`` and all others with an ``F``.
Diachromatic uses the information about enriched digests for quality control only
and passes it through to the interactions reported in the interaction file.
In this file,
an ``E`` corresponds to an ``T`` (enriched) and an ``N`` to an ``F`` (not enriched).

.. code-block:: console

    chr11    9641153    9642657   E   chr11   47259263   47272706   N   5:4

Digest file used for Javierre data
==================================

hg38 coordinates of enriched digests
------------------------------------

We used a
`list of enriched digests <https://osf.io/u8tzp/>`_
that can be downloaded as part of the publication by Javierre et al. 2016.
This list can be found in the *Capture design* in the archive named:

.. code-block:: console

    human_PCHiC_hg19_HindIII_design.tar.gz

that expands into a folder that can be provided to the interaction caller
`CHiCAGO <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4908757/>`_
as the *design folder*.
Along with four other files, this folder contains the following file:

.. code-block:: console

    Digest_Human_HindIII_baits_e75_ID.baitmap

This is
`CHiCAGO's bait map file <http://regulatorygenomicsgroup.org/resources/Chicago_vignette.html#input-files-required>`_
that contains the following columns:
``chr``, ``start``, ``end``, ``fragmentID``, ``geneName``.
Here are the first four lines of this file:

.. code-block:: console

    1	831895	848168	218	RP11-54O7.16;RP11-54O7.1
    1	848169	850618	219	RP11-54O7.2
    1	850619	874081	220	AL645608.1;RP11-54O7.3;SAMD11
    1	889424	903640	223	KLHL17;NOC2L;PLEKHN1

We extracted the coordinates of enriched digests from this file with 
the following awk command:

.. code-block:: console

    $ awk '{print "chr"$1"\t"$2"\t"$3}' Digest_Human_HindIII_baits_e75_ID.baitmap
    chr1	831895	848168
    chr1	848169	850618
    chr1	850619	874081
    chr1	889424	903640
    (...)

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
----------------

In order to create a Diachromatic digest file for the analysis of the Javierre data,
we first created a GOPHER project with the name ``no_digests_selected_HindIII``.
Then we set up the project for ``hg38``
(no need to download ``Transcripts`` and ``Alignability Map``)
and only selected the restriction enzyme ``HindIII`` for the design parameters.
Finally, we exported the following digest file:

.. code-block:: console

    no_digests_selected_HindIII_hg38_DigestedGenome.txt

Because GOPHER was not used to select capture probes,
no digest is marked as selected in this file.
We wrote a Python script to overwrite the values in the ``Selected`` column
of a digest file:

.. code-block:: console

    $ python diachrscripts/create_diachromatic_digest_file.py
       --enriched-digests-file Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed
       --diachromatic-digest-file no_digests_selected_HindIII_hg38_DigestedGenome.txt
       --out-prefix /JAV_hg38_HindIII

This script takes a BED file with coordinates of digests selected for enrichment (``--enriched-digests-file``)
and a Diachromatic digest file (``--diachromatic-digest-file``).
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

We also tried to correct the 48 cases.
To do this, we extracted a BED file from from the Diachromatic digest file as follows:

.. code-block:: console

    $ awk '{print $1"\t"$2"\t"$3}' no_digests_selected_HindIII_hg38_DigestedGenome.txt | tail -n+2 > no_digests_selected_HindIII_hg38_DigestedGenome.bed

Then we used BedTools to get all digests from the Diachromatic interaction file
that overlap at least 90% with an enriched digest.

.. code-block:: console

    $ bedtools intersect -f 0.90 -r -a no_digests_selected_HindIII_hg38_DigestedGenome.bed -b JAV_HindIII_hg38_digests_not_found.bed -wa

With the corrected BED file for enriched digests,
we get a Diachromatic digest file in which 22,045 are marked as enriched.
In our analysis of the Javierre data,
we used this file as input for Diachromatic.
