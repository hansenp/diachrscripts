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


## Definition of directed interactions

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

For each of the 44 subcategories, the script creates aa text file that contains an
interaction distance in each line. For the analysis at the interaction level,
the following files will be created:
```
# Directed interactions
<OUT_PREFIX>_di_ee_dist_array.tab
<OUT_PREFIX>_di_ne_dist_array.tab
<OUT_PREFIX>_di_en_dist_array.tab
<OUT_PREFIX>_di_nn_dist_array.tab

# Undirected reference interactions
<OUT_PREFIX>_uir_ee_dist_array.tab
<OUT_PREFIX>_uir_ne_dist_array.tab
<OUT_PREFIX>_uir_en_dist_array.tab
<OUT_PREFIX>_uir_nn_dist_array.tab

# Undirected  interactions
<OUT_PREFIX>_ui_ee_dist_array.tab
<OUT_PREFIX>_ui_ne_dist_array.tab
<OUT_PREFIX>_ui_en_dist_array.tab
<OUT_PREFIX>_ui_nn_dist_array.tab

# Directed simple interactions
<OUT_PREFIX>_di_ee_s_dist_array.tab
<OUT_PREFIX>_di_ne_s_dist_array.tab
<OUT_PREFIX>_di_en_s_dist_array.tab
<OUT_PREFIX>_di_nn_s_dist_array.tab

# Directed twisted interactions
<OUT_PREFIX>_di_ee_t_dist_array.tab
<OUT_PREFIX>_di_ne_t_dist_array.tab
<OUT_PREFIX>_di_en_t_dist_array.tab
<OUT_PREFIX>_di_nn_t_dist_array.tab
```

For the analysis at the read pair level, the following files will be created:
```
# Simple read pairs from directed interactions
<OUT_PREFIX>_di_ee_s_rp_dist_array.tab
<OUT_PREFIX>_di_ne_s_rp_dist_array.tab
<OUT_PREFIX>_di_en_s_rp_dist_array.tab
<OUT_PREFIX>_di_nn_s_rp_dist_array.tab

# Twisted read pairs from directed interactions
<OUT_PREFIX>_di_ee_t_rp_dist_array.tab
<OUT_PREFIX>_di_ne_t_rp_dist_array.tab
<OUT_PREFIX>_di_en_t_rp_dist_array.tab
<OUT_PREFIX>_di_nn_t_rp_dist_array.tab

# Simple read pairs from undirected reference interactions
<OUT_PREFIX>_uir_ee_s_rp_dist_array.tab
<OUT_PREFIX>_uir_ne_s_rp_dist_array.tab
<OUT_PREFIX>_uir_en_s_rp_dist_array.tab
<OUT_PREFIX>_uir_nn_s_rp_dist_array.tab

# Twisted read pairs from undirected reference interactions
<OUT_PREFIX>_uir_ee_t_rp_dist_array.tab
<OUT_PREFIX>_uir_ne_t_rp_dist_array.tab
<OUT_PREFIX>_uir_en_t_rp_dist_array.tab
<OUT_PREFIX>_uir_nn_t_rp_dist_array.tab

# Simple read pairs from undirected interactions
<OUT_PREFIX>_ui_ee_s_rp_dist_array.tab
<OUT_PREFIX>_ui_ne_s_rp_dist_array.tab
<OUT_PREFIX>_ui_en_s_rp_dist_array.tab
<OUT_PREFIX>_ui_nn_s_rp_dist_array.tab

# Simple read pairs from undirected interactions
<OUT_PREFIX>_ui_ee_t_rp_dist_array.tab
<OUT_PREFIX>_ui_ne_t_rp_dist_array.tab
<OUT_PREFIX>_ui_en_t_rp_dist_array.tab
<OUT_PREFIX>_ui_nn_t_rp_dist_array.tab
``` 

The files generated are read in by four R scripts for further analysis.

### Comparison of DI, UIR and UI interactions

The first R script is for comparing the distances between directed (DI), undirected
reference (UIR) and undirected interactions (UI) and is executed as follows:
```
Rscript --vanilla rscripts/07_analyze_interaction_distances/interaction_distances.r \
<OUT_DIR>/ \
<OUT_PREFIX> \
<OUT_PREFIX>_di_ee_dist_array.tab \
<OUT_PREFIX>_di_ne_dist_array.tab \
<OUT_PREFIX>_di_en_dist_array.tab \
<OUT_PREFIX>_di_nn_dist_array.tab \
<OUT_PREFIX>_uir_ee_dist_array.tab \
<OUT_PREFIX>_uir_ne_dist_array.tab \
<OUT_PREFIX>_uir_en_dist_array.tab \
<OUT_PREFIX>_uir_nn_dist_array.tab \
<OUT_PREFIX>_ui_ee_dist_array.tab \
<OUT_PREFIX>_ui_ne_dist_array.tab \
<OUT_PREFIX>_ui_en_dist_array.tab \
<OUT_PREFIX>_ui_nn_dist_array.tab
```
The first argument (`<OUT_DIR>/`) is the directory to which the results will be written.
The second argument (`<OUT_PREFIX>`) is used as prefix for the names of the generated files.
This is followed by 12 files that contain the interaction distances for DI, UIR, UI
within the enrichment categories `EE`, `NE`, `EN` and `NN`.

The script produces two files, a PDF file with histograms and a TSV file with
summary statistics.

The PDF file
```
<OUT_PREFIX>_i_distance_statistics_ee_ne_en_nn.pdf
```
contains a field of histograms for the distributions of distances
in the 12 subcategories.

![Interaction histograms](doc/07_analyze_interaction_distances/MK_i_distance_histograms_ee_ne_en_nn.png)

The first line of this fileds contains the histograms for DI (orange) and `EE`, `NE`, `EN` and `NN`.
The following two lines have the same structure and are for UIR (light blue) and UI (gray).
The axes of all histograms have the same range. Therefore, the histograms can be compared
with each other. Three summary statistics are determined for each subcategory:

1. Number of interactions (**n**)

2. Median interaction distance (**median**)

3. Interquartile ranges (**IQR**)
 
We use the IQR as a measure for dispersion.
The IQR is defined as the difference between the 75th (Q3) and 25th (Q1) percentiles,
i.e. the range that contains the middle 50% of the distances.
The IQR corresponds to the length of the box in a boxplot.
These summary statistics are shown in the legends of the histograms and are also written to
a TSV file with the name:
```
<OUT_PREFIX>_i_distance_statistics_ee_ne_en_nn.tsv
```

This file contains only one header line and one line with the corresponding 36 values.
```
OUT_PREFIX	DI_EE_N	DI_NE_N	DI_EN_N	DI_NN_N	UIR_EE_N	UIR_NE_N	UIR_EN_N	UIR_NN_N	UI_EE_N	UI_NE_N	UI_EN_N	UI_NN_N	DI_EE_MED	DI_NE_MED	DI_EN_MED	DI_NN_MED	UIR_EE_MED	UIR_NE_MED	UIR_EN_MED	UIR_NN_MED	UI_EE_MED	UI_NE_MED	UI_EN_MED	UI_NN_MED	DI_EE_IQR	DI_NE_IQR	DI_EN_IQR	DI_NN_IQR	UIR_EE_IQR	UIR_NE_IQR	UIR_EN_IQR	UIR_NN_IQR	UI_EE_IQR	UI_NE_IQR	UI_EN_IQR	UI_NN_IQR
MK	9343	97187	100013	2652	9157	99138	97917	2606	190206	2263030	2261998	116260	117502	113607	109861	81096	173721	187142	190400	59266	453700	348054	348125	39831	194838	186641	180064	140207	287543	318649	315731	116242	630739	467046	467140	74126
```
The TSV files will later be used to perform a combined analysis for all cell types.

### Comparison of simple and twisted interactions
XXX




## Analysis of interaction profiles

xxx




