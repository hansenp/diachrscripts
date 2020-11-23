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

## Calling interactions with Diachromatic

We used the Java program Diachromatic to derive interactions from
unprocessed capture Hi-C reads.
There is already a
[publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6678864/)
and a [detailed documentation](https://diachromatic.readthedocs.io/en/latest/index.html)
on this program.
This section is intended to give a brief overview
and to describe the details of the analysis of the capture Hi-C
dataset for 17 hematopoietic cell types.

The input of Diachromatic essentially consists of the following:

1. Paired-end read data in FASTQ format
2. A bowtie2 index for the reference sequence to which the reads will be mapped
3. The restriction enzyme that was used to generate the data
4. A file that contains the restriction fragments (or digests)
of the entire genome that result from digestion with the restriction enzyme

Diachromatic processes the data in three steps:

1. Truncation of reads
2. Mapping and removal of artifact read pairs
3. Counting of read pairs that map to the same digest pairs

In the following subsection you will be guided through
the alalysis with Diachromatic so that every step can be followed exactly.

### Input data

#### Paired-end data in FASTQ format

For paired-end data, the reads of the read pairs are typically in two different files,
one for the forward reads (R1) and one for the reverse reads (R2).
For the data on the hematopoietic cell types,
the downloadable reads for given replicates are divided into multiple FASTQ files
(beyond the division into R1 and R2).
For each replicate, we have merged the files by writing their contents
into two new files, one for R1 and one for R2.
We made sure that the files for R1 and R2 are in sync,
i.e. that two reads with the same row index form a pair.
For example, we have the following file pairs for replicate 1 on megakaryocytes:
```
EGAF00001274989.fq, EGAF00001274990.fq
EGAF00001274991.fq, EGAF00001274992.fq
EGAF00001274993.fq, EGAF00001274994.fq
```
In each line, the first file belongs to R1 annd the second file to R2.
We wrote these files into two new files as follows:
```
cat EGAF00001274989.fq EGAF00001274991.fq EGAF00001274993.fq > MK/JAV_MK_R10_1.fastq
cat EGAF00001274990.fq EGAF00001274992.fq EGAF00001274994.fq > MK/JAV_MK_R10_2.fastq
```
In table TODO.xls,
we have documented which original files were concatenated in which order.

#### bowtie2 index for reference sequence

Diachromatic uses bowtie2 to map reads to a reference sequence,
which requires an index for the reference sequence.
Such an index can be created with bowtie2 or a precalculated index can be downloaded.
For the data on the hematopoietic cell types,
we used the reference sequence hg38 and downloaded a precalculated index from here:
```
ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
```

#### Restriction enzyme an digest map

For the dataset on the 17 hematopoietic cell types,
the enzyme *HindIII* with the recognition sequence `AAGCTT` was used.
Diachromatic reports individual interactions as a pair of restriction fragments (or digest pairs)
along with the number of read pairs that map to this digest pair.
This requires the coordinates of all digests in the reference genome
that are passed to Diachromatic in form of a text file,
which we refer to as digest map.
For a given restriction enzyme and a reference genome,
a corresponding digest map can be created with the GOPHER software
[as described in the documentation](https://diachromatic.readthedocs.io/en/latest/digest.html).
To ensure consistency,
we recommend creating the digest map from the same FASTA file that was used
to create the bowtie2 index.

### Truncation of reads

Use Diachromatic to truncate the read pairs given in FASTQ format as follows:
```
java -jar Diachromatic.jar truncate \
   -e HindIII \
   -q MK/JAV_MK_R10_R1.fastq.gz \
   -r MK/JAV_MK_R10_R2.fastq.gz \
   -o MK \
   -x JAV_MK_R10
```

Diachromatic has an internal list of common restriction enzymes
and will use the appropriate recognition sequence and cutting positions
for `-e HindIII`.
We use the previously downloaded and concatenated FASTQ files
for the forward (R1, `-q`) and reverse (R2, `-r`) as input.
An already existing directory for the output (`-o`) and a prefix
for all generated files (`-x`) can also be specified.
For capture Hi-C data, we don't use the `--sticky-ends` option,
i.e. we assume that the sticky ends resulting from the restriction
have been filled in.
More details on the truncation of reads can be found in the
[relevant section of the Diachromatic documentation](https://diachromatic.readthedocs.io/en/latest/truncate.html).

### Mapping and removal of artifact read pairs

Use Diachromatic to map the the truncated read pairs to the reference sequence as follows:

```
java -jar Diachromatic.jar align \
-bsu \
-d <DIGEST_MAP> \
-q MK/JAV_MK_R10.truncated_R1.fastq.gz \
-r MK/JAV_MK_R10.truncated_R2.fastq.gz \
-b <BOWTIE2_PATH> \
-i <BOWTIE2_INDEX> \
-p 32 \
-j \
-o MK \
-x JAV_MK_R10
```

In addition to mapping, Diachromatic removes duplicated read pairs and
keeps track of the number of read pairs for different duplication levels.
Depending on the size of the input and the actual duplication rate,
this can take up a lot of memory.
We therefore recommend providing 16 to 32 GB memory.

We use the more stringent mode of Diachromatic to define uniquely mapped reads,
i.e. reads that map to only one location (`-bsu`).
In order to determine artifact read pairs,
for example pairs mapped to the same digest,
the previously created digest map is required (`-d`).

Diachromatic uses bowtie2 to map the reads to the reference genome.
To do this,
an executable bowtie2 file and an index for the reference must be specified (`-b`, `-i`).
We use 32 threads for the maapping with bowtie2 (`-p`).

For possible subsequent investigation,
we write the rejected artifact read pairs to an extra BAM file (`-j`).
The output can be redirected and given prefixes as with the `truncate` command.
More details on the mapping and removal of artifact read pairs can be found in the
[relevant section of the Diachromatic documentation](https://diachromatic.readthedocs.io/en/latest/mapping.html).

### Counting of read pairs that map to the same digest pairs

Use Diachromatic to count read pairs between interacting digest regions as follows:

```
java -jar Diachromatic.jar count \
-d <DIGEST_MAP>  \
-v JAV_MK_R10.valid_pairs.aligned.bam \
-s \
-o MK \
-x JAV_MK_R10
```

In Diachromatic, interactions are defined as digest pairs that have at least
one supporting read pair.
In this step, the supporting read pairs for individual interactions are counted,
for which the digest map is required (`-d`).

We use the unique valid pairs from the previous step as input,
i.e. duplicates and artifact read pairs have been removed.
We use the `-s` option so that the simple and twisted read pairs counts
of individual interactions are reported separately.

More details on counting read pairs between interacting digest regions can be found in the
[relevant section of the Diachromatic documentation](https://diachromatic.readthedocs.io/en/latest/count.html).

The interactions with their read pair counts are written to the following file:

```
MK/JAV_MK_R10.interaction.counts.table.tsv
```

These are the first few lines from such a file:
```
chr1    46297999   46305684   A   chr1    51777391   51781717   I   2:1
chr17   72411026   72411616   I   chr17   72712662   72724357   I   3:2
chr7    69513952   69514636   I   chr7    87057837   87061499   A   4:3
chr11    9641153    9642657   I   chr11   47259263   47272706   A   5:4
```
Each line represents one interaction.

Columns 1 to 3 and 5 to 7 contain the coordinates of the digest pair,
whereby the smaller coordinates are always in columns 1 to 3.

In column 4 and 8 there is either an `A` or an `I`,
where column 4 belongs to the first and column 8 belongs to the second digest.
An `A` means that the corresponding digest was selected for target enrichment
and an `I` means that it was not selected.
The information about digests that were selected for enrichment
is taken from the digest map that was generated with GOPHER.
When generating the digest map with GOPHER,
we used the shortcut option *All protein coding genes*
because the analyzed data comes from a whole-promoter capture Hi-C
experiment.
This has the effect that all digests which contain at least on TSS
and are potentially suited for enrichment are marked with an `A`
and all others with an `I`.
The digests marked with an `A` do not exactly correspond to the digests
that were actually selected for the experiment.
This is due to different annotations as well as different criteria
for the selection of enrichable digests.
Therefore, we only used the markings with `A` and `I` at the beginning
to roughly asses the how well the enrichment worked.
For the following,
we did not use these markings with `A` and `I`,
but instead a list of enriched digests from the original publication
of Javierre et al. 2016 (see below).
To make it easier to distinguish,
we denote enriched enriched digests with an `E` and non-enriched
digests with an `N`, when using the annotation from this list.

The last columnn in a Diachromatic interaction file
shows the counts of simple and twisted read pairs
separated by a colon.
For example, `5:4` means that five simple and four twisted
read pairs were counted for an interaction.

### Subsequent filtering of interactions

XXX

## hg38 coordinates of enriched digests

For all subsequent analyzes,
we used a [list of enriched digests](https://osf.io/u8tzp/)
that can be downloaded as part of the publication by Javierre et al. 2016.
This list can be found in the *Capture design* in the archive named:
```
human_PCHiC_hg19_HindIII_design.tar.gz
```
that expands into a folder that can be provided to Chicago as the
*design folder*.

Along with four other files, this folder contains the following file:
```
Digest_Human_HindIII_baits_e75_ID.baitmap
```
This is
[Chicago's bait map file](http://regulatorygenomicsgroup.org/resources/Chicago_vignette.html#input-files-required)
that contains the following columns:
`chr`, `start`, `end`, `fragmentID`, `geneName`.
Here are the first four lines of this file:
```
1	831895	848168	218	RP11-54O7.16;RP11-54O7.1
1	848169	850618	219	RP11-54O7.2
1	850619	874081	220	AL645608.1;RP11-54O7.3;SAMD11
1	889424	903640	223	KLHL17;NOC2L;PLEKHN1
```
From this file, we extracted the coordinates of enriched digests
by adding a `chr` to the chromosome numbers in column 1 and
writing them them together with the start and end positions to a new file.
These are the first four lines of the file with the extracted digest coordinates:
```
$ awk '{print "chr"$1"\t"$2"\t"$3}' Digest_Human_HindIII_baits_e75_ID.baitmap | head -n 4
chr1	831895	848168
chr1	848169	850618
chr1	850619	874081
chr1	889424	903640

```
In total, there are coordinates for 22,076 digests.
These coordinates refer to the genome build hg19.
We used
[UCSC's LiftOver tool](https://genome.ucsc.edu/cgi-bin/hgLiftOver)
to convert the coordinates to hg38.
22,056 digest could be converted successfully to hg38.
The conversion failed for 20 digests
because hg19 coordinates in hg38
are either split or partially deleted.
The resulting file in BED format cann be found here:
```
diachrscripts/additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed
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

For the dataset that we analyzed, there are three to four biological replicates
for each cell type.
On the other hand, we need as many read pairs as possible in the individual
interactions in order to decide whether an interaction is directed or not.
We therefore decided on an compromise when combining interactions from different
replicates.
We discard interactions that occur in less than two replicates and,
for the remaining interactions, we add up the read pair counts from
all replicates separately for simple and twisted read pairs. This way
of combining interactions from different replicates is implemented in
the following script:
```
01_combine_interactions_from_replicates.py
```
This script has the following parameters:
```
--out-prefix <string>
```
The name of each created file will have this prefix, which can also contain the path to an already existing directory.
```
--interaction-files-path <string>
```
Path to a directory that contains gzipped files in Diachromatic's interaction format.
From the files output by Diachromatic,
we have filtered out interactions between different chromosomes,
interactions with a distance of less than 20,000 bp and
interactions with or on chromosome `chrM` (see above).
```
--required-replicates <int>
```
Required number of replicates. All interactions that occur in less replicates
will be discarded.
For the remaining interactions,
the simple and twisted read pair counts from different replicates
will be added up separately.

Diachromatic does not filter interactions
and even outputs interactions that have only a single read pair.
On the other hand, when combining interactions,
the interactions from multiple replicates must be read into memory.
Therefore, the memory consumption can become very high
and we carried out this step on a compute cluster.


We have prepared small input files for testing
so that this step can be followed here.
There are a total of four replicates and four interactions.
The first interaction occurs only in replicate 1,
the second interaction occurs in replicates 1 and 2,
the third interaction occurs in replicates 1,2 and 3 and
the fourth interactions occurs in all replicates.

This is content of the interaction file for replicate 4:
```
chr1    46297999   46305684   A   chr1    51777391   51781717   I   2:1
chr17   72411026   72411616   I   chr17   72712662   72724357   I   3:2
chr7    69513952   69514636   I   chr7    87057837   87061499   A   4:3
chr11    9641153    9642657   I   chr11   47259263   47272706   A   5:4
```
From this file, we created the files for the other three replicates
by deleting interactions from the last line one by one.
By creating the files in this way,
the individual interactions have same simple and twisted read pair counts
for all replicates, which is usually not the case.
However, in simplifies the presentation here,
because we only need to know the content of the file for replicate 4
in order to understand the results.

To get the combined interactions that occur in at least two replicates
execute the following command:
```
$PYTHON_PATH 01_combine_interactions_from_replicates.py \
--out-prefix TEST \
--interaction-files-path tests/data/test_01/ \
--required-replicates 2
```

This will generate two files:
```
TEST_at_least_in_2_replicates_summary.txt
TEST_at_least_in_2_replicates_interactions.tsv.gz
```
The first file contains an overview of the numbers of interactions
in the individual files and the combined interactions.
The second file contains the combined interactions:
```
chr1    46297999   46305684   A   chr1    51777391   51781717   I   8:4
chr17   72411026   72411616   I   chr17   72712662   72724357   I   9:6
chr7    69513952   69514636   I   chr7    87057837   87061499   A   8:6
```
The interaction on chromosome `chr11` does not occur in this file
because it was observed for replicate 4 only.
However, we required that an interaction must have been observed in at least two replicates.
The interaction on chromosome `chr7` occurs in the files for replicate 3 and 4.
Since this interaction has the same read pair counts for both replicates,
the counts in the file for combined interactions double
(`4:3` becomes `8:6`).
The interaction on chromosome `chr17` occurs in the files for replicate 2, 3 and 4
and the counts triple (`3:2` becomes `9:6`).
Finally, the interaction on `chr11` occurs in the files for all four replicates
and the counts quadruple (`2:1` becomes `8:4`).

## Permuation of simple and twisted read pairs (02)

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




