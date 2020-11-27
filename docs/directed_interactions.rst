.. _RST_c:

#####################
Directed interactions
#####################

Remodeling plans
================

Until now, the P-values (logsf) are still calculated in the script:

.. code-block:: console

    diachrscripts/04_extract_gene_symbols_and_tss.py

In addition,
the scripts generates the Enhanced Interaction (EI) format.
The script is complicated and takes a long time to assign gene annotations
to digests of interactions.
At the moment, we don't need any gene annotations and want to do without
the EI-Format.

Based on the P-values, the interactions are divided into directed (``DI``)
and undirected interactions (``UI``).

.. code-block:: console

    diachrscripts/05_define_di_uie_and_uii.py

We further distinguish between ``UI`` that have a digest that is also involved in
a directed interaction (``UII``) and other interactions (``UIE``).

In a next step, we select undirected reference interactions from ``UI``,
whereby we to not distinguish between ``UII`` and ``UIE``.

.. code-block:: console

    diachrscripts/06_select_uir_from_uie_and_uii.py

After that, we have three categories of interactions: ``DI``, ``UIR``, ``UI``.

If we discard the EI format, we could combine the categorization of innteractions
into one single script.

Current implementation
======================

We carried out the first steps of the analysis on a computing cluster.
This include this step, which is implemented in the following script:

.. code-block:: console

    diachrscripts/05_define_di_uie_and_uii.py

This script uses the previously determined P-value threshold
in order to classify interactions as either directed or undirected.
The result is a file in EI format in which each line stands for one interaction
and the third column indicates whether the interaction is directed or undirected.

To reproduce the analysis steps below on your computer,
download the following files:

.. code-block:: console

    JAV_MK_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
    JAV_ERY_RALT_0.0018_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
    JAV_NEU_RALT_0.0011_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
    JAV_MON_RALT_0.0012_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
    JAV_MAC_M0_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
    JAV_MAC_M1_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
    JAV_MAC_M2_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
    JAV_EP_RALT_0.0017_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
    JAV_NB_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
    JAV_TB_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
    JAV_FOET_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
    JAV_NCD4_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
    JAV_TCD4_RALT_0.002_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
    JAV_ACD4_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
    JAV_NACD4_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
    JAV_NCD8_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
    JAV_TCD8_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz


Place these files in the directory:
```
diachrscripts/results/05_define_directed_interactions
```

Each of these files contains the combined interactions for one of the hematopoietic 17 cell types.