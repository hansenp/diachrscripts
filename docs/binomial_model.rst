.. _RST_binomial_model:

##############
Binomial model
##############

We use a binomial test in order to assess imbalances of simple and twisted
read pairs within individual interactions for statistical significance,
whereby the total number of read pairs (n) corresponds to a parameter
of the binomial distribution (*number of trials*).
Therefore, the test is based on different null distributions for interactions with
different n.
In this section, we review our implementation of the P-value calculation
and investigate the consequences of the different null distributions.

For our analyzes, we have implented various functions in the class `BinomialInteractionModel`.

.. code-block:: console

    diachrscripts/diachr/binomial_interaction_model.py

In a Jupyter notebook, these functions are explained step by step and used to carry out the entire analysis:

.. code-block:: console

    diachrscripts/jupyter_notebooks/binomialModel.ipynb

In this notebook, a total 15 plots are generated.
Alternatively, the following script can be executed to generate these plots:

.. code-block:: console

    $ python diachrscripts 00_explore_binomial_model.py
       --out-prefix results/00_explore_binomial_model/EBM

If no enhanced interaction file is passed to the script, only 8 plots are generated.
In order to generate all plots, execute the following command:

.. code-block:: console

    $ python diachrscripts 00_explore_binomial_model.py
       --out-prefix results/00_explore_binomial_model/EBM
       --enhanced-interaction-file MK_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz

Note that column 3 in the interaction file must indicate whether the interactions
are directed (DI), undirected reference (UIR) or undirected interactions (UI).
This classification is made in the script `06_select_uir_from_uie_and_uii.py`.
If you have not yet created interaction files with reference
interactions, the execute the first of the two commands and otherwise
the second.

If the commands are executed as stated above,
all plots will be written to the following directory:

.. code-block:: console

    results/00_explore_binomial_model

