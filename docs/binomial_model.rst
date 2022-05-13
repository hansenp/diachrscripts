.. _RST_binomial_model:

##############
Binomial model
##############

We use a binomial test in order to assess the imbalances
in the four read pair counts of individual interactions for statistical significance,
whereby the total number of read pairs corresponds the number of trials (n) and
the sum of the two highest counts to the number of successes (k).
For our null model, we assume a probability of success of p=0.5.

For different total numbers of reads (n), the test is based on different null
distributions and has different power.
In addition, the binomial distribution is a discrete distribution.
In this section, we examine the implications of this.

For our analyzes, we have implemented some methods in the class `BinomialInteractionModel`.

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

If no Diachromatic interaction file is passed to the script, only 8 plots are generated.
In order to generate all plots, execute the following command:

.. code-block:: console

    $ python diachrscripts 00_explore_binomial_model.py
       --out-prefix results/00_explore_binomial_model/EBM
       --diachromatic-interaction-file MK_evaluated_and_categorized_interactions.tsv.gz

Note that column 11 in the interaction file must indicate whether the interactions
are directed (DI), undirected reference (UIR) or undirected interactions (UI).
If you have not yet created interaction files with P-values and interaction categories,
execute the first of the two commands and otherwise the second.

If the commands are executed as stated above,
all plots will be written to the following directory:

.. code-block:: console

    results/00_explore_binomial_model

