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

For our publication we used data for which you first have to sign a European Genome-phenome Archive (EGA) data usage
agreement. But in the `Sequence Read Archive (SRA) <https://www.ncbi.nlm.nih.gov/sra/docs/>`_ many datasets are freely
available. In order to download data from there, you must first
`download, install and configure the SRA Toolkit <https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit>`_.

Alternatively, data can also be downloaded from the `European Nucleotide Archive (ENA) <https://www.ebi.ac.uk/ena/browser/home>`_.

Replicate 1 - Forward:

.. code-block:: console

    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR436/ERR436029/ERR436029_1.fastq.gz -O ERR436029_1.fastq.gz

Replicate 1 - Reverse:

.. code-block:: console

    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR436/ERR436029/ERR436029_2.fastq.gz -O ERR436029_2.fastq.gz

.. code-block:: console

    #OUT_PREFIX="MIFSUD_R10"
    #SRR_LIST="ERR436029"

    #OUT_PREFIX="MIFSUD_R20"
    #SRR_LIST="ERR436028 ERR436030 ERR436033"

    OUT_PREFIX="MIFSUD_R30"
    SRR_LIST="ERR436031 ERR436026"

    > $OUT_DIR/$OUT_PREFIX\_1.fastq
    > $OUT_DIR/$OUT_PREFIX\_2.fastq
    for SRR in $SRR_LIST;
    do
            wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR436/$SRR/$SRR\_1.fastq.gz -O $OUT_DIR/$SRR\_1.fastq.gz
            wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR436/$SRR/$SRR\_2.fastq.gz -O $OUT_DIR/$SRR\_2.fastq.gz
            zcat $OUT_DIR/$SRR\_1.fastq.gz >> $OUT_DIR/$OUT_PREFIX\_1.fastq
            zcat $OUT_DIR/$SRR\_2.fastq.gz >> $OUT_DIR/$OUT_PREFIX\_2.fastq
    done
    gzip $OUT_DIR/$OUT_PREFIX\_1.fastq
    gzip $OUT_DIR/$OUT_PREFIX\_2.fastq

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


