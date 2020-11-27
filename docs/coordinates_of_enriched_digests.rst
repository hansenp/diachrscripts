.. _RST_coordinates_of_enriched_digests:

###############################
Coordinates of enriched digests
###############################

When analyzing data from capture Hi-C experiments,
it is important to know which restriction fragments or digests were selected for enrichment.
In this section,
we describe two ways of introducing this information into the analysis with
``Diachromatic`` and ``diachrscripts``.

GOPHER's digest file
====================

If the probes for a capture Hi-C experiment have been created with GOPHER,
then the digest file, which
`can be exported from the corresponding GOPHER design <https://diachromatic.readthedocs.io/en/latest/digest.html>`_,
contains precise information about all enriched digests.
These are the first few line of such a digest file:

.. code-block:: console

    Chromosome      Fragment_Start_Position Fragment_End_Position   Fragment_Number 5'_Restriction_Site     3'_Restriction_Site     Length  5'_GC_Content   3'_GC_Content   5'_Repeat_Content       3'_Repeat_Content       Selected        5'_Probes       3'_Probes
    chr1    1       16007   1       None    HindIII 16007   0.000   0.000   0.000   0.003   F       0       0
    chr1    16008   24571   2       HindIII HindIII 8564    0.018   0.018   0.000   0.015   F       0       0
    chr1    24572   27981   3       HindIII HindIII 3410    0.046   0.046   0.000   0.044   F       0       0
    chr1    27982   30429   4       HindIII HindIII 2448    0.035   0.035   0.047   0.043   F       0       0

Each line stands represents one digest,
and the entire file contains the digests for the entire genome.
In the ``Selected`` column,
enriched digests are marked with a ``T`` and all others with an ``F``.
Diachromatic uses the information about enriched digests for quality control only
and passes it through to the interactions reported in the interaction file.
In this file,
an ``A`` corresponds to an ``T`` (enriched) and an ``I`` to an ``F`` (not enriched).

Digest file used for Javierre data
==================================

For analyzing the data from Javierre 2016,
we used the GOPHER's shortcut option *All protein coding genes*
to generate a digest map
because the analyzed data comes from a whole-promoter capture Hi-C
experiment.
This has the effect that all digests which contain at least on TSS
and are potentially suited for enrichment are marked with an ``A``
and all others with an ``I``.
In this case,
the digests marked with an ``A`` do not exactly correspond to the digests
that were actually selected for the experiment.
This is due to different annotations as well as different criteria
for the selection of enrichable digests.
Therefore, we only used the markings with ``A`` and ``I`` at the beginning
to roughly asses the how well the enrichment worked.
For further analyzes,
we did not use these markings with ``A`` and ``I``,
but instead a list of enriched digests from the original publication
of Javierre et al. 2016 (see below).


hg38 coordinates for Javierre 2016
==================================

For all subsequent analyzes,
we used a
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

From this file, we extracted the coordinates of enriched digests
by adding a ``chr`` to the chromosome numbers in column 1 and
writing them them together with the start and end positions to a new file.
These are the first four lines of the file with the extracted digest coordinates:

.. code-block:: console

    $ awk '{print "chr"$1"\t"$2"\t"$3}' Digest_Human_HindIII_baits_e75_ID.baitmap | head -n 4
    chr1	831895	848168
    chr1	848169	850618
    chr1	850619	874081
    chr1	889424	903640

In total, there are coordinates for 22,076 digests.
These coordinates refer to the genome build ``hg19``.
We used
`UCSC's LiftOver tool <https://genome.ucsc.edu/cgi-bin/hgLiftOver>`_
to convert the coordinates to ``hg38``.
22,056 digest could be converted successfully to ``hg38``.
The conversion failed for 20 digests
because ``hg19`` coordinates in ``hg38``
are either split or partially deleted.
The resulting file in BED format can be found here:

.. code-block:: console

    diachrscripts/additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed

If an analysis with ``diachrscripts`` needs information about enriched digests,
we use this BED file as input.
In order to distinguish better from the Diachromatic noation with A and I,
we denote enriched enriched digests with an ``E`` and non-enriched
digests with an ``N``, when using the annotation from this list.
