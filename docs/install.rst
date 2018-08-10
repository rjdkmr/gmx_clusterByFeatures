Download and Installation
=========================

It can be downloaded using git as follows

.. code-block:: bash

  git clone -b master https://github.com/rjdkmr/gmx_clusterByFeatures

It can be also downloaded as `zip file <https://github.com/rjdkmr/gmx_clusterByFeatures/archive/master.zip>`_.

******

Requirements
------------
It depends on following two packages:
  * **GROMACS** : 2016 and above version
  * **Python** : 3.4 and above version

******

GROMACS
~~~~~~~

A standard installation of GROMACS is sufficient. GROMACS library
(``libgromacs.a`` or ``libgromacs.so``) and header files are required for compilation.

If GROMACS is not installed at standard location, ``-DGMX_PATH`` can be used to locate
GROMACS installation. e.g. ``-DGMX_PATH=/opt/gromacs``

******

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

******

Installation
------------

Clone the repository from github as directed above then follow these steps.

.. code-block:: bash

  cd gmx_clusterByFeatures # or gmx_clusterByFeatures-master (zip file download)
  mkdir build
  cd build
  cmake -DGMX_PATH="/opt/gromacs/gmx2016.5" ..
  make
  sudo make install


Now, gmx_clusterByFeatures command will be accessible in terminal.
