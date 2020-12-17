.. _RST_undirected_reference_interactions:

#################################
Undirected reference interactions
#################################

Our binomial test for directionality rates interactions with fewer
read pairs less often as significant that interactions with many read pairs
(:ref:`RST_binomial_model`).
Therefore, we select undirected reference interactions that are comparable
to directed interactions with respect to the distribution of read pair numbers.

Test file
=========

We have prepared aa Diachromatic interaction file to test
the selection of undirected reference interactions:

.. code-block:: console

    tests/data/test_04/diachromatic_interaction_file.tsv

This file contains 18 interactions that are clearly directed,
with three interactions in the enrichment category ``NN``,
four in ``NE``, five in ``EN`` and six in ``EE``.
The total number of read pairs per interaction ranges from 101 to 106.

.. code-block:: console

    chr14	43059116	43059494	N	chr14	43101212	43101810	N	100:1
    chr8	129042054	129044258	N	chr8	129121269	129121986	N	100:2
    chr15	73467156	73468652	N	chr15	73526903	73528438	N	100:3

    chr17	72411026	72411616	N	chr17	72712662	72724357	E	100:1
    chr18	38724804	38726198	N	chr18	76794986	76803172	E	100:2
    chr11	114362648	114362686	N	chr11	114396073	114404234	E	100:3
    chr15	56158017	56158267	N	chr15	56462978	56465983	E	100:4

    chr14	34714080	34716362	E	chr14	50135355	50139051	N	100:1
    chr1	91022201	91023797	E	chr1	116561813	116566655	N	100:2
    chr1	15681566	15697108	E	chr1	19411358	19417940	N	100:3
    chr12	13586326	13591414	E	chr12	101116206	101119138	N	100:4
    chr2	113676580	113686263	E	chr2	202796295	202797013	N	100:5

    chr7	25228385	25228778	E	chr7	42234764	42240281	E	100:1
    chr6	117448209	117455412	E	chr6	152981947	152991467	E	100:2
    chr7	123787495	123793134	E	chr7	141015018	141017643	E	100:3
    chr10	100185111	100188716	E	chr10	100911854	100914842	E	100:4
    chr10	73301814	73302637	E	chr10	112212615	112226039	E	100:5
    chr8	11066997	11069837	E	chr8	11473293	11485455	E	100:6

In addition, there are 16 interactions that are clearly undirected and
match the directed interactions in terms of enrichment category and
total number of read pairs per interaction.
These are the undirected reference interactions.

.. code-block:: console

    chr1	68750175	68752699	N	chr1	68795323	68796512	N	50:51
    chr10	24701259	24703782	N	chr10	24993210	24994662	N	50:52
    chr13	77661563	77666713	N	chr13	77716416	77717610	N	50:53

    chr11	86392415	86393959	N	chr11	125833282	125834157	E	50:51
    chr4	143079578	143081309	N	chr4	146945582	146955674	E	50:52
    chr14	97520062	97521143	N	chr14	99641271	99654714	E	50:53

    chr12	92715804	92721658	E	chr12	97575775	97576200	N	50:51
    chr3	149470667	149475062	E	chr3	150148068	150150290	N	50:52
    chr18	31373121	31377878	E	chr18	53036132	53037663	N	50:53
    chr2	16841469	16846902	E	chr2	17018322	17019298	N	50:54
    chr14	92748029	92749630	E	chr14	96578956	96580826	N	50:55

    chr3	185935734	185943372	E	chr3	194132402	194139770	E	50:51
    chr5	65953726	65970373	E	chr5	71582365	71588637	E	50:52
    chr5	138607182	138612232	E	chr5	151368315	151371525	E	50:53
    chr1	23789605	23791509	E	chr1	153951062	153962812	E	50:54
    chr5	112154251	112162055	E	chr5	115769718	115774925	E	50:55

The file also contains 12 interactions,
which are also clearly undirected,
but the number of read pair per interaction is in a different range (50 to 52),
which is why these interactions cannot be used as reference interactions.

.. code-block:: console

    chr8	110169057	110171420	N	chr8	110203244	110203772	N	25:25
    chr21	24706442	24708010	N	chr21	24802647	24806351	N	25:26
    chr1	21038219	21041845	N	chr1	21139524	21140811	N	25:27

    chr7	105237629	105238078	N	chr7	128862293	128872009	E	25:25
    chr6	89716777	89716882	N	chr6	89916736	89920828	E	25:26
    chr9	120719242	120724112	N	chr9	135466848	135469617	E	25:27

    chr14	24167761	24173500	E	chr14	75079406	75087746	N	25:25
    chr12	132700658	132715264	E	chr12	133104166	133106284	N	25:26
    chr8	142668614	142673077	E	chr8	142903855	142904484	N	25:27

    chr17	80806561	80813742	N	chr17	80880355	80886553	N	25:25
    chr2	197200527	197202476	N	chr2	197253179	197256110	N	25:26
    chr2	139005552	139005761	N	chr2	139026311	139029719	N	25:27

For two directed interactions,
there is no matching undirected reference interaction.
One interaction is in enrichment category ``NE`` and has a total of 104 read pairs:

.. code-block:: console

    chr15	56158017	56158267	N	chr15	56462978	56465983	E	100:4

The second interaction is in the enrichment category ``EE`` and has a total of 106 read pairs:

.. code-block:: console

    chr8	11066997	11069837	E	chr8	11473293	11485455	E	100:6


Old
===

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

