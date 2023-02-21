# diachrscripts

This is a collection of Python (3.8) modules, scripts and Jupyter notebooks that can be used to investigate
Diachromatic interactions in terms of the four counts for the different orientations of mapped paired-end reads.

## Setup

First you need to clone this repository with:
```
$ git clone https://github.com/TheJacksonLaboratory/diachrscripts.git
```

If `virtualenv` is not installed on your system, install it with:
```
$ pip install virtualenv
```

Then change into the cloned directory, create a virtual environment and activate it:
```
$ cd diachrscripts
$ virtualenv dscripts
$ source dscripts/bin/activate
```

Install all required packages into this environment:
```
(dscripts)$ pip install -r requirements.txt
```
This may take a few minutes to complete.

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
the RTD documentation as follows:
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
