.. _RST_Rate_and_categorize_interactions:

################################
Rate and categorize interactions
################################

In this script, the P-values for directionality of interactions are calculated and, based on these P-values,
the interactions are categorized into directed (DI) und undirected interactions (UI).
Undirected reference interactions (UIR) are then selected from UI, which are comparable to DI with respect to
to the total number of read pairs (n).

The input consists of a file in Diachromatic interaction format and a P-value threshold that was determined using our
FDR procedure. The output is again a file in Diachromatic interaction format, but there are two additional columns
on the right for the calculated P-values and interaction categories.

We use a two-sided binomial test to asses the directionality of interactions.

Based on the P-values, the interactions are divided into directed (``DI``)
and undirected interactions (``UI``).

In a next step, we select undirected reference interactions from ``UI``,
whereby we to not distinguish between ``UII`` and ``UIE``.

.. code-block:: console

    diachrscripts/04_rate_and_categorize_interactions.py

After that, we have three categories of interactions: ``DI``, ``UIR``, ``UI``.