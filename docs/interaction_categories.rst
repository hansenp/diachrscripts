.. _RST_interaction_categories:

######################
Interaction categories
######################

To compare directed and undirected interactions, we have decomposed the into different categories.
At the top level, we differentiate between directed (``DI``), undirected reference (``UIR``),
and undirected interactions (``UI``).

We further break down directed interactions into the categories *simple* and *twisted*,
depending which of the two types of read pairs predominates.

At the lowest level, we differentiate interactions with regard to their enrichment status,
i.e. whether both digests (``EE``), the second digest (``NE``), the first digest (``EN``)
or neither of the two digests (``NN``) were selected for target enrichment.
Note that the distinction between ``NE`` and ``EN`` is also a distinction with regard to the direction,
i.e. whether an interaction, viewed from the enriched digest, goes to the left or right.

For undirected interactions, a distinction between *simple* and *twisted* is not possible,
as neither of the two read types predominates.
In order to extend our comparison between *simple* and *twisted* to undirected interactions,
we made an additional analysis at the level of read pairs.
For this purpose, we keep the division into ``DI``, ``UIR`` and ``UI``,
but read pairs between given digests are no longer combined into interactions.
If, for example, an interaction has 29 simple and 27 twisted read pairs,
then the interaction distance goes 29 times into the corresponding distribution of distances for simple
and 27 times into the distance distribution for twisted.

At the level of interactions, the decomposition results in 20 subcategories:

1. **Directed interactions** are decomposed into simple and twisted interactions.
There is also a category that contains all directed interactions,
i.e. simple and twisted interactions.
Within these three categories, we differentiate interactions according to the enrichment
status of the digests involved (``EE``, ``NE``, ``EN`` or ``NN``). This results in 12 categories.

2. **Undirected reference interactions** cannot be decomposed into simple and twisted interactions
and are decomposed according the enrichment status, resulting in four categories.

3. **Undirected interactions** cannot be decomposed into simple annd twisted innteractions either,
and there are only four caategories for the enrichment states of involved digests.

.. image:: img/interaction_categories.png

At the level of read pairs, all three categories can be decomposed into simple and twisted
interactions, which are further differentiated according to enrichment status,
which results in 24 subcategories.

The script:

.. code-block:: console

    diachrscripts/07_analyze_interaction_distances.py
       --enhanced-interaction-file <OUT_PREFIX>_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz
       --out-prefix <OUT_PREFIX>

implements the decomposition of interactions and creates one file for each subcategory
that contains one distance in each line.

The name of each created file will have the same prefix (``--out-prefix``),
which can also contain the path to an already existing directory.

A file in enhanced interaction format (``--enhanced-interaction-file``) that was created with the script ``06_select_uir_from_uie_and_uii.py``.
The interaction distances are taken from the second column of this file.

The third column contains the tag for the interaction category that
is either ``DI``, ``UIR`` or ``UI``.
Whether an interaction is *simple* or *twisted* in decided on the basis of the read
pair counts in the fifth column.
The enrichment states are taken from the sixth column.