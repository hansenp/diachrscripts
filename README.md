# diachrscripts

This is a collection of Python (3.8) modules, scripts and Jupyter notebooks that can be used to investigate
Diachromatic interactions in terms of the four counts for the different types of mapped paired-end reads.

## Setup

First you need to clone this repository with:
```
$ git clone https://github.com/TheJacksonLaboratory/diachrscripts.git
```

If `virtualenv` is not installed on your system, install it with:
```
pip install virtualenv
```

Then change into the cloned directory, create a virtual environment and activate it:
```
$ virtualenv dscripts
$ source dscripts/bin/activate
```

Install all required packages into this environment:
```
(dscripts)$ pip install -r requirements.txt
```

To make the environment available in Jupyter notebooks create a Jupyter kernel.
```
(dscripts)$ pip install ipykernel
(dscripts)$ python -m ipykernel install --user --name=dscripts
```
When using Jupyter notebooks select the kernel `dscripts`.

## Documentation

Change into the `docs` directory and generate
the RTD documentation as follows:
```
(dscripts)$ make html
```
You should now find documentation here: `diachrscripts/docs/_build/html/index.html`

## Testing

Run all tests as follows:
```
(dscripts)$ python -m unittest discover tests/
```
To deactivate the environment use `deactivate`.

