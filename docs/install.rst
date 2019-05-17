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

**Not require to install GROMACS**

**Only available on Linux, MacOS-10.12 (Sierra), MacOS-10.13 (High Sierra) and MacOS-10.14 (Mojave)**

On Linux
~~~~~~~~

1. Python3 is available through package managers such as **yum** (Fedora, CentOs), **YaST** (OpenSuse) and **apt-get**
   (Ubuntu, Linux Mint). For example on ubuntu: run ``sudo apt-get install python3`` command to install Python3.

2. Install **gmx_clusterByFeatrues** by ``sudo pip3 install gmx-clusterByFeatures`` command.



On MacOS
~~~~~~~~

1. Python3 is available through |Homebrew| package manager. After installing Homebrew, run ``brew install python3`` command to install Python3.

2. Install **gmx_clusterByFeatrues** by ``pip3 install gmx-clusterByFeatures`` command.


.. note:: Presently, installation with pip on MacOS is restricted to **10.12 (Sierra)**, **10.13 (High Sierra)** 
          and 10.14 (Mojave) versions. For other MacOS versions, install gmx_clusterByFeatrues from source as 
          described further below.


Updating gmx_clusterByFeatrues
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To update the gmx_clusterByFeatrues package use following command:

.. code-block:: bash

    pip3 install --upgrade --no-deps gmx-clusterByFeatures


``--upgrade`` flag is used to update the package and ``--no-deps`` prevents
update of dependent packages like numpy, scipy, matplotlib etc.


****



Installation from source-code
-----------------------------

Requirements
~~~~~~~~~~~~~~

It depends on following two packages:
  * **GROMACS** : 2016 and above version
  * **Python** : 3.4 and above version


GROMACS
~~~~~~~

A standard installation of GROMACS is sufficient. GROMACS library
(``libgromacs.a`` or ``libgromacs.so``) and header files are required for compilation.

If GROMACS is not installed at standard location, define ``GMX_PATH`` environment variable as follows:

.. code-block:: bash

    export GMX_PATH=/path/to/installed/gromacs
    
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


Four other Python packages |numpy|, |scipy|, |scikit-learn|, and |matplotlib| are required 
that can be installed as follows:

.. code-block:: bash

  sudo pip3 install numpy scipy sklearn matplotlib


****

Downloading source-code
~~~~~~~~~~~~~~~~~~~~~~~~~

It can be downloaded using git as follows

.. code-block:: bash

  git clone -b master https://github.com/rjdkmr/gmx_clusterByFeatures

It can be also downloaded as `zip file <https://github.com/rjdkmr/gmx_clusterByFeatures/archive/master.zip>`_.



Compilation and Installation using python
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Clone the repository from github as directed above then follow these steps.

.. code-block:: bash

  cd gmx_clusterByFeatures # or gmx_clusterByFeatures-master (zip file download)
  export GMX_PATH=/path/to/installed/gromacs
  sudo GMX_PATH=$GMX_PATH python3 setup.py install


Now, gmx_clusterByFeatures command will be accessible in terminal.

Compilation and Installation using cmake for C++ IDEs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This method can be used for developement purpose using C++ IDE like QT creator and KDevelop etc.

**To install and use gmx_clusterByFeatrues from source location:**

.. code-block:: bash

  cd gmx_clusterByFeatures # or gmx_clusterByFeatures-master (zip file download)
  mkdir build
  cd build
  cmake -DGMX_PATH=/path/to/installed/gromacs -DINPLACE=ON ..
  make
  sudo make install  # Only needed for first time install 

In this installation, only ``gmx_clusterByFeatures`` executable file is installed at default
location (mostly ``/usr/local/bin``) while whole package remains at the source location.

This method is extremely useful for developement because ``make install`` is only required
for first time to install executable file. During subsequent developement, only command 
``make`` need to be repeated. In IDEs ``make`` command is executed by ``build``. In IDEs 
project build setting, cmake arguments ``-DGMX_PATH=/path/to/installed/gromacs -DINPLACE=ON``
needs to be added manually.
