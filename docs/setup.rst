.. _rstsetup

Setup
#####

First of all you need to clone this repository with:

.. code-block:: console

    git clone https://github.com/TheJacksonLaboratory/diachrscripts.git

To make a smaller requirements file, we start off like this

.. code-block:: console

    conda create -n p3dia

To check the environment, enter

.. code-block:: console

    conda activate p3dia

To add the packages we need for diachromatic-scripts,

.. code-block:: console

    pip install --user --requirement requirements.txt


We now want to put this kernel into jupypter. We will use the ipykernel package for this.

.. code-block:: console

    pip install ipykernel
    python -m ipykernel install --user

