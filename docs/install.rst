.. |numpy| raw:: html

  <a href="http://www.numpy.org/" target="_blank"> numpy </a>

.. |scipy| raw:: html

  <a href="http://www.scipy.org/" target="_blank"> scipy </a>

.. |matplotlib| raw:: html

  <a href="http://matplotlib.org/" target="_blank"> matplotlib </a>

.. |h5py| raw:: html

  <a href="http://www.h5py.org/" target="_blank"> h5py </a>

.. |scikit-learn| raw:: html

  <a href="https://scikit-learn.org/" target="_blank"> scikit-learn </a>

.. |Homebrew| raw:: html

  <a href="http://brew.sh/" target="_blank"> Homebrew </a>


Download and Installation
=========================
  

Quick Installation using ``pip``
--------------------------------

It is **recommended** method to install gmx_clusterByFeatrues.

Following Python versions are supported:

* Python 3.9
* Python 3.10
* Python 3.11
* Python 3.12

**No dependency on GROMACS. Install and use it.**

On Linux
~~~~~~~~

Distributions with **glibc version 2.28 or later** are supported. Some of the supported distributions are:

* Debian 10+
* Ubuntu 18.10+
* Fedora 29+
* CentOS/RHEL 8+

Use following steps to install gmx_clusterByFeatrues:

.. code:: bash

    sudo python3 -m pip install gmx-clusterByFeatures


On MacOS
~~~~~~~~

Python3 is available through |Homebrew| package manager.

Currently, MacOS versions **12.0+** and **13.0+** versions are supported.

.. code:: bash

    python3 -m pip install gmx-clusterByFeatures


Updating gmx_clusterByFeatrues
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To update the gmx_clusterByFeatrues package use following command:

.. code-block:: bash

    python3 -m pip install --upgrade --no-deps gmx-clusterByFeatures


``--upgrade`` flag is used to update the package and ``--no-deps`` prevents
update of dependent packages like numpy, scipy, matplotlib etc.


****



Installation from source-code
-----------------------------

Installation from source-code is recommended through `conda` environment.

1. Clone the repository:

.. code:: bash

    git clone --recursive https://github.com/rjdkmr/gmx_clusterByFeatures.git

2. Change directory to ``gmx_clusterByFeatures``:
  
.. code:: bash

  cd gmx_clusterByFeatures

3. Create conda environment and install dependencies:
   
.. code:: bash

  conda env create -y --prefix ./venv --file environment.yaml

4. Activate the environment:

.. code:: bash

  conda activate ./venv

5. Run the following script to install local GROMACS and subsequently install gmx_clusterByFeatures:

.. code:: bash

  bash -i scripts/build_dev_setup_conda.sh