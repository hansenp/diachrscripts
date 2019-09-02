# diachrscripts

This is a collection of Python (3.7) scripts for the analysis of simple and twisted read pairs and interactions.
These scripts can be executed in Jupyter notebooks that are also included in this repository    .

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






