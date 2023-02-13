######################
diachrscripts tutorial
######################

***************************************************************************************************
Analysis of interactions with respect to the four relative orientations of mapped paired-end reads
***************************************************************************************************

After mapping of Hi-C or capture Hi-C paired-end reads, four relative orientations of read pairs can be distinguished
depending on the order and strands of the two reads.
We have extended our ``Diachromatic`` software to report the read pair counts for each interaction separately
by orientation.
In ``diachrscripts``, we implemented analysis and visualization methods in Python modules that can be used
in scripts and Jupyter notebooks to study interactions with respect to imbalances in the four read pair counts.
This tutorial provides an overview of the entire workflow, from downloading the paired-end data,
through the initial processing with ``Diachromatic``, pooling interactions from different replicates with ``pooler.py``,
calling of unbalanced interactions with ``UICer.py``,
to various analyses in Jupyter notebooks.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   tutorial
   diachromatic_input_preparation
   resources
