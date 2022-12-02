.. _RST_jupyter_notebook:

#################
Jupyter notebooks
#################

******************************************************
Performing various analyzes on unbalanced interactions
******************************************************

We have implemented all analyzes following the calling of unbalanced interactions in different Jupyter notebooks.
The ``DiachromaticInteractionSet`` is the central data structure in all of these analyzes.
It can be created from an interaction file generated with ``UICer``.

.. code-block:: python

    from diachr import DiachromaticInteractionSet
    d11_interaction_set = DiachromaticInteractionSet()
    d11_interaction_set.parse_file(
        i_file = "MIF_REPC/MIF_REPC_evaluated_and_categorized_interactions.tsv.gz",
        verbose = True)

An interaction file generated with ``UICer.py`` for the Mifsud data can be downloaded as follows:

.. code-block:: console

    $ wget https://www.genecascade.org/downloads/diachrscripts/MIF_REPC_evaluated_and_categorized_interactions.tsv.gz


Read type and configuration frequencies
=======================================

See this
`notebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/publication/read_type_and_configuration_frequencies.ipynb>`__.

Triangle interaction maps
=========================

See this
`notebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/publication/triangle_interaction_maps.ipynb>`__.

Classification of baited digests
================================

See this
`notebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/publication/baited_digest_analysis_1.ipynb>`__.

Analysis of BDC0, BDC1 and BDC2 digests
=======================================

See this
`notebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/publication/baited_digest_analysis_2.ipynb>`__.

Distance-dependent interaction frequencies
==========================================

See this
`notebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/publication/interaction_frequency_distance_analysis_1.ipynb>`__
and this
`one <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/publication/interaction_frequency_distance_analysis_2.ipynb>`__.
