# diachrscripts

This is a collection of Python modules, scripts and Jupyter notebooks that can be used to investigate
Diachromatic interactions in terms of the four counts for the different orientations of mapped paired-end reads.

## Setup

First you need to clone this repository with:
```
$ git clone https://github.com/TheJacksonLaboratory/diachrscripts.git
```

Check which Python version you are using:
```
python --version
```
For convenience, you could do something like this: `alias python=python3.9`.
We tested Python versions `3.7`, `3.8`, `3.9` and `3.10` on macOS.
For Python version `3.11` we get an 
[error](https://github.com/pysam-developers/pysam/issues/1154)
when installing the `pysam` package.


Then change into the cloned directory, create a virtual environment and activate it:
```
$ cd diachrscripts
$ python -m venv dscripts
$ source dscripts/bin/activate
```
We use the Python package `pybedtools`,
which requires BedTools to be installed and available through the `$PATH` variable. 
To install BedTools follow
[these instructions](https://bedtools.readthedocs.io/en/latest/content/installation.html).

Install all required packages into this environment:
```
(dscripts)$ pip install -r requirements.txt
```
This may take a few minutes to complete.
Depending on your Python installation,
you may get an error message saying that `Python.h` is missing.
In this case, install `python-devel` and execute the last command again.

## Jupyter notebooks

To make the environment available in Jupyter notebooks create a Jupyter kernel.
```
(dscripts)$ pip install ipykernel
(dscripts)$ python -m ipykernel install --user --name=dscripts
```
When using Jupyter notebooks select the kernel `dscripts`.

Most of the analysis in `diachrscripts` is done in Jupyter notebooks.
To get started with the analyses, change into the directory `jupyter_notebooks` and launch the notebook server.
```
(dscripts)$ cd jupyter_notebooks
(dscripts)$ jupyter-notebook
```
Depending on the configurations, a browser window will then open automatically,
showing the contents of the directory from which the notebook server was started.
If this is not the case, copy the URL that was output to the terminal and paste it into your browser's address bar.
In the browser, select the notebook
[Get_started.ipynb](jupyter_notebooks/Get_started.ipynb).

## Documentation

Change into the `docs` directory and generate
the ReadTheDocs documentation as follows:
```
(dscripts)$ cd docs
(dscripts)$ make html
```
You should now find documentation here: `diachrscripts/docs/_build/html/index.html`

## Testing

Run all tests as follows:
```
(dscripts)$ python -m unittest discover tests/
```
To deactivate the environment use `deactivate`.
