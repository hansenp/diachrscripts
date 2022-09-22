.. _RST_jupyter_notebook:

#################
Jupyter notebooks
#################

******************************************************
Performing various analyzes on unbalanced interactions
******************************************************

We have implemented all analyzes following the calling of unbalanced interactions in different Jupyter Notebooks.
The ``DiachromaticInteractionSet`` is the central data structure in all of these analyzes.
It can be created from an interaction file generated with ``DICer``.

.. code-block:: python

    from diachr import DiachromaticInteractionSet
    d11_interaction_set = DiachromaticInteractionSet()
    d11_interaction_set.parse_file(
        i_file = "MIF_REPC/MIF_REPC_evaluated_and_categorized_interactions.tsv.gz",
        verbose = True)


Frequencies of read types and configurations of interactions
============================================================

See this
`notebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/read_pair_and_interaction_types.ipynb>`__.

Representation of interactions in triangle heatmaps
===================================================

See this
`notebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/dtvis.ipynb>`__.

Classification of baited digests
================================

See this
`notebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/interactions_at_baited_digests_select_baited_digests.ipynb>`__.

Interaction distances
=====================

See this
`notebook <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/interaction_frequency_distance_analysis.ipynb>`__
and this
`one <https://github.com/TheJacksonLaboratory/diachrscripts/blob/develop/jupyter_notebooks/interaction_frequency_distance_analysis_2.ipynb>`__.
