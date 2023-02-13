.. _RST_Resources:

#########
Resources
#########

******
GOPHER
******

`GOPHER <https://gopher.readthedocs.io/en/latest/>`_
is a tool that can be used to design baits for capture Hi-C experiments.
In the context of this work, we used it to prepare digest files as input for Diachromatic.

************
Diachromatic
************

`Diachromatic <https://diachromatic.readthedocs.io/en/latest/>`_ is a tool for quality control and initial processing
of Hi-C and CHi-C data.
We extended Diachromatic so that counts of supporting read pairs of interactions are
reported separately by relative orientation of mapped paired-end reads.

*************
diachrscripts
*************

`diachrscripts <https://github.com/TheJacksonLaboratory/diachrscripts/>`_ is a collection of Python modules,
scripts, and Jupyter notebooks that can be used to analyze interactions with respect to imbalances in the
four counts.
Diachromatic interaction files can be passed to the UICer.py script that produces ``Diachromatic11`` files
with two additional columns, one for imbalance scores and one for interaction categories (balanced or unbalanced).
All further analyzes in diachrscripts are based on ``Diachromatic11`` files.







