# diachrscripts

This is a collection of Python (3.7) scripts for the analysis of simple and twisted read pairs and interactions.

## Setup

First of all you need to clone this repository and change into the ```diachscripts``` directory with:
```shell script
$ git clone https://github.com/TheJacksonLaboratory/diachrscripts.git
$ cd diachscripts
```

We used ```conda``` to manage the virtual environment.
On a Linux machine, setup the environment as follows:
```shell script
$ conda env create -f environment_linux_p37env.yml
```
On a Mac, setup the environment with a different yml file:
```shell script
$ conda env create -f environment_mac_p37env.yml
```
Active the environment as follows:
```shell script
$ conda activate diachscripts_p37env
```

Install an ipython kernel for the Jupyter notebooks with:
```shell script
python -m ipykernel install --user --name diachscripts_p37env --display-name "Python 3 (diachscripts_p37env)"
```

From this environment start Jupyter notebooks with:
```shell script
(diachscripts_p37env) $ jupyter-notebook
```

Open one of the notebooks in ```diachscripts/jupyter-notebooks``` and choose the right kernel via:

```Kernel -> Change kernel```

Follow the instructions given in the notebook in order to get the required input data and to perform the analyses.

You can select the environment in ```PyCharm``` via:

```Preferences -> Project interpreter -> Add -> Existing environment -> </YOUR/ANACONDA/PATH>/anaconda2/envs/diachrscripts_env/bin/python```

Don't forget to update the YML file after changing the environment.
For this purpose, u         se:
```shell script
(diachscripts_p37env) $ conda env export --no-builds > environment_linux_p37env.yml
```


## Set up (new)

To make a smaller requirements file, we start off like this
```
conda create -n p3dia
```
To check the environment, enter
```
conda activate p3dia
```
To add the packages we need for diachromatic-scripts,
```
pip install --user --requirement requirements.txt
```

We now want to put this kernel into jupypter. We will use the ipykernel package for this.
```
pip install ipykernel
python -m ipykernel install --user
```

## Exploration of the binomial model for directionality of interactions (00)

We use a binomial test in order to assess imbalances of simple and twisted
read pairs within individual interactions for statistical significance,
whereby the total number of read pairs (n) corresponds to a parameter
of the binomial distribution (*number of trials*).
Therefore, the test is based on different null distributions for interactions with
different n.
In this section, we review our implementation of the P-value calculation
and investigate the consequences of the different null distributions.

For our analyzes, we have implented various functions in the class `BinomialInteractionModel`.
```
diachr/binomial_interaction_model.py
```

In a Jupyter notebook, these functions are explained step by step and used to carry out the entire analysis:

```
jupyter_notebooks/binomialModel.ipynb
```
In this notebook, a total 15 plots are generated.
Alternatively, the following script can be executed to generate these plots:
```
PYTHON_PATH="<PATH_TO_PYTHON_3.7>/python3.7"
$PYTHON_PATH 00_explore_binomial_model.py --out-prefix results/00_explore_binomial_model/EBM
```
If no enhanced interaction file is passed to the script, only 8 plots are generated.
In order to generate all plots, execute the following command:
```
$PYTHON_PATH 00_explore_binomial_model.py --out-prefix results/00_explore_binomial_model/EBM \
--enhanced-interaction-file results/06_select_reference_interactions/MK/MK_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz
```
Note that column 3 in the interaction file must indicate whether the interactions
are directed (DI), undirected reference (UIR) or undirected interactions (UI).
This classification is made in the script `06_select_uir_from_uie_and_uii.py`.
If you have not yet created interaction files with reference
interactions, the execute the first of the two commands and otherwise 
the second.

If the commands are executed as stated above,
all plots will be written to the following directory:
```
results/00_explore_binomial_model
```

## Combination of reproducible interactions from different replicates (01)

XXX

## Definition of directed interactions (05)

We carried out the first steps of the analysis on a computing cluster.
This include this step, which is implemented in the following script:
```
diachrscripts/05_define_di_uie_and_uii.py
```
This script uses the previously determined P-value threshold
in order to classify interactions as either directed or undirected.
The result is a file in EI format in which each line stands for one interaction
and the third column indicates whether the interaction is directed or undirected.

To reproduce the analysis steps below on your computer,
download the following files:
```
JAV_MK_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
JAV_ERY_RALT_0.0018_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
JAV_NEU_RALT_0.0011_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
JAV_MON_RALT_0.0012_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
JAV_MAC_M0_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
JAV_MAC_M1_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
JAV_MAC_M2_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
JAV_EP_RALT_0.0017_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
JAV_NB_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
JAV_TB_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
JAV_FOET_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
JAV_NCD4_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
JAV_TCD4_RALT_0.002_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
JAV_ACD4_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
JAV_NACD4_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
JAV_NCD8_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
JAV_TCD8_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz
```

Place these files in the directory:
```
diachrscripts/results/05_define_directed_interactions
```

Each of these files contains the combined interactions for one of the hematopoietic 17 cell types.

## Selection of undirected reference interactions

Our binomial is underpowered for interactions with few read pairs.
Therefore, we select undirected reference interactions that are comparable
to directed interactions with respect to the distribution of read pair numbers.
This  is implemented in the  script:
```
diachrscripts/06_select_uir_from_uie_and_uii.py
```
Once you have downloaded the files from the previous step and placed them in the
directory in the designated directory, you can select the reference interactions for
each of the 17 cell types by executing the following bash script:
```
bash_scripts/06_select_reference_interactions/06_select_reference_interactions_run_4_all.sh
```
This script writes the results to the directory:
```
results/06_select_reference_interactions/
```
whereby a separate subdirectory is created for each cell type.

## Analysis of interaction distances

To compare directed and undirected interactions, we have decomposed the into different categories.
At the top level, we differentiate between directed (*DI*), undirected reference (*UIR*), and undirected interactions (*UI*).
We further break down directed interactions into the categories *simple* and *twisted*,
depending which of the two types of read pairs predominates.
At the lowest level, we differentiate interactions with regard to their enrichment status,
i.e. whether both digests (*EE*), the second digest (*NE*), the first digest (*EN*)
or neither of the two digests (*NN*) were selected for target enrichment.
Note that the distinction between *NE* and *EN* is also a distinction with regard to the direction,
i.e. whether an interaction, viewed from the enriched digest, goes to the left or right.

For undirected interactions, a distinction between *simple* and *twisted* is not possible,
as neither of the two read types predominates.
In order to extend our comparison between *simple* and *twisted* to undirected interactions,
we made an additional analysis at the level of read pairs.
For this purpose, we keep the division into *DI*, *UIR* and *UI*,
but read pairs between given digests are no longer combined into interactions.
If, for example, an interaction has 29 simple and 27 twisted read pairs,
then the interaction distance goes 29 times into the corresponding distribution of distances for simple
and 27 times into the distance distribution for twisted.

At the level of interactions, the decomposition results in 20 subcategories:

1. **Directed interactions** are decomposed into simple and twisted interactions.
There is also a category that contains all directed interactions,
i.e. simple and twisted interactions.
Within these three categories, we differentiate interactions according to the enrichment
status of the digests involved (*EE*, *NE*, *EN* or *NN*). This results in 12 categories.

2. **Undirected reference interactions** cannot be decomposed into simple and twisted interactions
and are decomposed according the enrichment status, resulting in four categories.

3. **Undirected interactions** cannot be decomposed into simple annd twisted innteractions either,
and there are only four caategories for the enrichment states of involved digests.

![Decompose interactions](doc/07_analyze_interaction_distances/interaction_decomposition.png)

At the level of read pairs, all three categories can be decomposed into simple and twisted
interactions, which are further differentiated according to enrichment status,
which results in 24 subcategories.

The script:
```
diachrscripts/07_analyze_interaction_distances.py
```
implements the decomposition of interactions and creates one file for each subcategory
that contains one distance in each line.
 
The script has the following parameters:
```
--out-prefix
```
The name of each created file will have this prefix,
which can also contain the path to an already existing directory.
```
--enhanced-interaction-file <OUT_PREFIX>_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz
```
A file in enhanced interaction format that was created with the script `06_select_uir_from_uie_and_uii.py`.
The interaction distances are taken from the second column of this file.
The third column contains the tag for the interaction category that
is either `DI`, `UIR` or `UI`.
Whether an interaction is *simple* or *twisted* in decided on the basis of the read
pair counts in the fifth column.
The enrichment states are taken from the sixth column.









## Analysis of interaction profiles

xxx




