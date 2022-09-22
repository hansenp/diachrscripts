`diachrscripts`
===============

Analysis of `Diachromatic` interactions with separate counts for the four types of mapped paired-end reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After mapping of Hi-C and capture Hi-C paired-end reads,
four types of read pairs can be distinguished depending on
their relative orientation and which strand they were
mapped to.
This information was previously used only to remove artifacts.
We have extended our ``Diachromatic`` software to report
the counts of the remaining valid read pairs
for each interaction separately by type.
In ``diachscripts`` we have implemented methods that can be used
to analyze interactions with respect to the four counts.
The tutorial provides an overview of the entire workflow and possible
analyzes is therefore a good starting point for your own analyzes.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   tutorial
   diachromatic_input_preparation
   interaction_pooling
   unbalanced_interaction_calling
   jupyter_notebooks
   binomial_model






