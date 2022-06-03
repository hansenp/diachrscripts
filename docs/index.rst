`diachrscripts`
===============

Analysis of `Diachromatic` interactions with separate counts for the four types of mapped paired-end reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After mapping of Hi-C and capture Hi-C paired-end reads,
four types of read pairs can be distinguished depending on
their relative orientation and which strand they were
mapped to.
This information was previously used only to remove artifacts.
We have extended our Diachromatic software to report
the counts of the remaining valid read pairs
for each interaction separately by type.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   interaction_calling
   coordinates_of_enriched_digests
   combining_interactions
   tutorial
   binomial_model
   simple_twisted_randomization
   false_discovery_rate
   rate_and_categorize_interactions






