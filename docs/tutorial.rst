.. _RST_tutorial:

########
Tutorial
########

This tutorial describes the entire workflow from downloading the data, calling the interactions with ``Diachromatic``,
pooling interactions from different replicates, calling unbalanced interactions with ``DICer`` up to the various
analyzes of interactions with unbalanced read counts that can be performed in Jupyter Notebooks.

.. image:: img/analysis_flowchart.png

************************************************
Downloading paired-end Hi-C or capture Hi-C data
************************************************

Use the script
`dumpy.sh <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/additional_scripts/dumpy.sh>`__
to download promoter capture Hi-C data on GM12878 cells
`(Mifsud et al. 2015) <https://pubmed.ncbi.nlm.nih.gov/25938943/>`_.

.. code-block:: console

    $ dumpy.sh MIF_R1 MIF_R1 "ERR436029"
    $ dumpy.sh MIF_R2 MIF_R2 "ERR436028 ERR436030 ERR436033"
    $ dumpy.sh MIF_R3 MIF_R3 "ERR436031 ERR436026"

A separate directory is created for each of the three replicates.
Because it is paired-end data, there are pairs of forward and reverse FASTQ files with
suffixes ``_1`` and ``_2``.
The lines of the two files are in sync, meaning that the forward read from line ``i`` in the first file
and the reverse read from line i in the second file form a pair.


.. code-block:: console

    $ ls -l MIF_R*/
    XXX


Both files contain the same number of reads.

******************************************
Calling interactions with ``Diachromatic``
******************************************

This is described here: :ref:`RST_Interaction_calling`.

**********************************************
Pooling interactions from different replicates
**********************************************

This is described here: :ref:`RST_Combining_interactions`.

**********************************************
Calling unbalanced interactions with ``DICer``
**********************************************

So far, this is only described in this
`Jupyter Notebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/Demonstration_of_DICer.ipynb>`__.


******************************************************
Performing various analyzes on unbalanced interactions
******************************************************

We have implemented all analyzes following the calling of unbalanced interactions in different Jupyter Notebooks.

Interaction distances
=====================

See this
`motebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/interaction_frequency_distance_analysis.ipynb>`__
and this
`one <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/interaction_frequency_distance_analysis_2.ipynb>`__.


Frequencies of read types and configurations of interactions
============================================================

See this
`motebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/read_pair_and_interaction_types.ipynb>`__.

Representation of interactions in triangle heatmaps
===================================================

See this
`motebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/dtvis.ipynb>`__.

Classification of baited digests
================================

See this
`motebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/interactions_at_baited_digests_select_baited_digests.ipynb>`__.

TAD boundaries
==============

See this
`motebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/tad_boundaries.ipynb>`__.


