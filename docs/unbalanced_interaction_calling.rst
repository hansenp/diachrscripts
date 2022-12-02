.. _RST_Unbalanced_interaction_calling:

##############################
Unbalanced interaction calling
##############################

We implemented the calling of unbalanced interactions in the Python script ``UICer.py``.
As input, this script requires a file in Diachromatic interaction format.
If a P-value threshold is specified, then this will be used for the classification of interactions.
Otherwise, a randomization procedure is used to set the P-value threshold such that the FDR
remains below a chosen threshold of 5\%.
The output also consists of a file in Diachromatic interaction format extended by two columns for P-values and
interaction categories. This file contains only interactions that are powered at the P-value threshold
used for classification.


Running the script ``UICer.py``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the following command to run the script: ::

    $ diachrscripts/UICer.py \
       --interaction-files-path gzdir/ \
       --out-prefix out_dir/out_prefix \
       --xxxxxx


Available arguments (ToDo: Complete table):

+---------------+-------------------------------+-------------------------+-----------+---------------------------------------------------------------------------------------------+-----------------+
| Short option  | Long option                   | Example                 | Required  | Description                                                                                 | Default         |
+===============+===============================+=========================+===========+=============================================================================================+=================+
| ``-i``        | ``--interaction-files-path``  | ``gzdir``               | yes       | Path to directory with Diachromatic interaction files                                       | --              |
+---------------+-------------------------------+-------------------------+-----------+---------------------------------------------------------------------------------------------+-----------------+

By default, the final P-value threshold is determined via randomization. If a P-value is specified,
then this P-value threshold will be used and no randomizations will be performed.

Output files
~~~~~~~~~~~~

``UICer.py`` generates a total of seven files:

- ``out_dir/out_prefix_reports.txt``
- ``out_dir/out_prefix_randomization_plot.pdf``
- ``out_dir/out_prefix_randomization_table.txt``
- ``out_dir/out_prefix_randomization_histogram_at_threshold.pdf``
- ``out_dir/out_prefix_randomization_histogram_at_001.pdf``
- ``out_dir/out_prefix_randomization_histogram_at_005.pdf``
- ``out_dir/out_prefix_evaluated_and_categorized_interactions.tsv.gz``
