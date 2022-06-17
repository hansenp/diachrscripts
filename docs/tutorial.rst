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

Use the following code to download promoter capture Hi-C data on GM12878 cells
`(Mifsud et al. 2015) <https://pubmed.ncbi.nlm.nih.gov/25938943/>`_.
There are three replicates. Because it is paired-end data, there are pairs of forward and reverse FASTQ files with
suffixes ``_1`` and ``_2``. Both files contain the same number of reads.

.. code-block:: console

    #!/bin/bash

    #OUT_DIR=$1
    #OUT_PREFIX=$2
    #SRR_LIST=$3

    OUT_DIR="MIF_R1"
    OUT_PREFIX="MIF_R1"
    SRR_LIST="ERR436029"

    #OUT_DIR="MIF_R2"
    #OUT_PREFIX="MIF_R2"
    #SRR_LIST="ERR436028 ERR436030 ERR436033"

    #OUT_DIR="MIF_R3"
    #OUT_PREFIX="MIF_R3"
    #SRR_LIST="ERR436031 ERR436026"

    mkdir $OUT_DIR
    > $OUT_DIR/$OUT_PREFIX\_1.fastq # Forward
    > $OUT_DIR/$OUT_PREFIX\_2.fastq # Reverse
    for SRR in $SRR_LIST;
    do
            # Forward
            wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR436/$SRR/$SRR\_1.fastq.gz -O $OUT_DIR/$SRR\_1.fastq.gz
            zcat $OUT_DIR/$SRR\_1.fastq.gz >> $OUT_DIR/$OUT_PREFIX\_1.fastq
            md5sum $OUT_DIR/$SRR\_1.fastq.gz

            # Reverse
            wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR436/$SRR/$SRR\_2.fastq.gz -O $OUT_DIR/$SRR\_2.fastq.gz
            zcat $OUT_DIR/$SRR\_2.fastq.gz >> $OUT_DIR/$OUT_PREFIX\_2.fastq
            md5sum $OUT_DIR/$SRR\_2.fastq.gz
    done
    gzip $OUT_DIR/$OUT_PREFIX\_1.fastq # Forward
    md5sum $OUT_DIR/$OUT_PREFIX\_1.fastq
    gzip $OUT_DIR/$OUT_PREFIX\_2.fastq # Reverse
    md5sum $OUT_DIR/$OUT_PREFIX\_2.fastq


+-----------------------+----------------------------------------------------------+
| Dataset               | MD5 checksum                                             |
+=======================+==========================================================+
| ``ERR436029_1``       | XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                         |
+-----------------------+----------------------------------------------------------+
| ``ERR436029_2``       | XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                         |
+-----------------------+----------------------------------------------------------+
| ``ERR436028_1``       | XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                         |
+-----------------------+----------------------------------------------------------+
| ``ERR436030_1``       | XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                         |
+-----------------------+----------------------------------------------------------+
| ``ERR436033_1``       | XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                         |
+-----------------------+----------------------------------------------------------+
| ``ERR436028_2``       | XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                         |
+-----------------------+----------------------------------------------------------+
| ``ERR436030_2``       | XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                         |
+-----------------------+----------------------------------------------------------+
| ``ERR436033_2``       | XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                         |
+-----------------------+----------------------------------------------------------+
| ``ERR436031_1``       | XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                         |
+-----------------------+----------------------------------------------------------+
| ``ERR436026_1``       | XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                         |
+-----------------------+----------------------------------------------------------+
| ``ERR436031_2``       | XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                         |
+-----------------------+----------------------------------------------------------+
| ``ERR436026_2``       | XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                         |
+-----------------------+----------------------------------------------------------+


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


