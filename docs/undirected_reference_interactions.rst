.. _RST_undirected_reference_interactions:

#################################
Undirected reference interactions
#################################

Our binomial test for directionality rates interactions with fewer
read pairs less often as significant that interactions with many read pairs
(:ref:`RST_binomial_model`).
Therefore, we select undirected reference interactions that are comparable
to directed interactions with respect to the distribution of read pair numbers.





This  is implemented in the  script:

.. code-block:: console

    $ python diachrscripts/06_select_uir_from_uie_and_uii.py \
       --enhanced-interaction-file X \
       --out-prefix X


Once you have downloaded the files from the previous step and placed them in the
directory in the designated directory, you can select the reference interactions for
each of the 17 cell types by executing the following bash script:

.. code-block:: console

    bash_scripts/06_select_reference_interactions/06_select_reference_interactions_run_4_all.sh

This script writes the results to the directory:

.. code-block:: console

    results/06_select_reference_interactions/

whereby a separate subdirectory is created for each cell type.

