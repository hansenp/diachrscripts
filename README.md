# diachrscripts

This is a collection of Python (2.7) scripts for the analysis of simple and twisted read pairs and interactions.
These scripts can be executed in a Jupyter notebook that is also included in this repository    .

First of all you need to clone this repository and change into the ```diachscripts``` directory with:
```shell script
$ git clone https://github.com/TheJacksonLaboratory/diachrscripts.git
$ cd diachscripts
```

We used ```conda``` to manage the virtual environment.
Setup and activate the environment as follows:
```shell script
$ conda env create -f environment.yml
$ conda activate diachscripts_env
```

From this environment start the Jupyter notebook with:
```shell script
(diachscripts) $ jupyter-notebook
```

Follow the instructions given in the notebook in order to get the required input data and to perform the analyses.

You can select the environment in ```PyCharm``` via:

```Preferences -> Project interpreter -> Add -> Existing environment -> </YOUR/ANACONDA/PATH>/anaconda2/envs/diachrscripts_env/bin/python```

Don't forget to update the YML file after changing the environment.
For this purpose, use:
```shell script
(diachscripts) $ conda env export > environment.yml
```


# Intallation (new)
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






