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


Requirements and Installation
=============================


Requirements
------------

It depends on following two packages:
  * **GROMACS** : 2016 and above version
  * **Python** : 3.4 and above version

**Python packages required for installation:** These packages are installed automatically during gmx_clusterByFeatures installation.

* |scikit-learn|
* |numpy|
* |scipy|
* |matplotlib|


****

Installation using ``pip``
--------------------------

It is **recommended** method to install gmx_clusterByFeatrues.

**Not require to install GROMACS**

**Only available on Linux, MacOS-10.12 and MacOS-10.13**

On Linux
~~~~~~~~

1. Python3 is available through package managers such as **yum** (Fedora, CentOs), **YaST** (OpenSuse) and **apt-get**
   (Ubuntu, Linux Mint). For example on ubuntu: run ``sudo apt-get install python3`` command to install Python3.

2. Install **gmx_clusterByFeatrues** by ``sudo pip3 install gmx-clusterByFeatrues`` command.



On MacOS
~~~~~~~~

1. Python3 is available through |Homebrew| package manager. After installing Homebrew, run ``brew install python3`` command to install Python3.

2. Install **gmx_clusterByFeatrues** by ``pip3 install gmx-clusterByFeatrues`` command.


Updating gmx_clusterByFeatrues
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To update the gmx_clusterByFeatrues package use following command:

.. code-block:: bash

    pip3 install --upgrade --no-deps gmx-clusterByFeatrues


``--upgrade`` flag is used to update the package and ``--no-deps`` prevents
update of dependent packages like numpy, scipy, matplotlib etc.


****



Installation from source-code
-----------------------------

It can be downloaded using git as follows

.. code-block:: bash

  git clone -b master https://github.com/rjdkmr/gmx_clusterByFeatures

It can be also downloaded as `zip file <https://github.com/rjdkmr/gmx_clusterByFeatures/archive/master.zip>`_.


GROMACS
~~~~~~~

A standard installation of GROMACS is sufficient. GROMACS library
(``libgromacs.a`` or ``libgromacs.so``) and header files are required for compilation.

If GROMACS is not installed at standard location, define ``GMX_PATH`` environment variable as follows:

.. code-block:: bash

    export GMX_PATH=/opt/gromacs
    
Python3
~~~~~~~

To compile gmx_clusterByFeatures, Python3 developement files should be installed previously.

On Debian like distribution (Debian, Ubuntu, Linux Mint etc.), which uses apt as
package manager, python3-development files can be installed as follows:

.. code-block:: bash

  sudo apt-get install python3 python3-dev


On OS such as fedora/centos/RHEL, which uses yum as package manager, python3-
development files can be installed as follows:

.. code-block:: bash

  sudo yum install python3 python3-devel


Two other Python packages are required that can be installed as follows:

.. code-block:: bash

  sudo pip3 install sklearn matplotlib


Compilation and Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Clone the repository from github as directed above then follow these steps.

.. code-block:: bash

  cd gmx_clusterByFeatures # or gmx_clusterByFeatures-master (zip file download)
  sudo GMX_PATH=$GMX_PATH python3 setup.py install


Now, gmx_clusterByFeatures command will be accessible in terminal.
