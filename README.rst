
.. image:: https://travis-ci.org/rjdkmr/gmx_clusterByFeatures.svg?branch=master
    :target: https://travis-ci.org/rjdkmr/gmx_clusterByFeatures

.. image:: https://readthedocs.org/projects/gmx-clusterbyfeatures/badge/?version=latest
    :target: https://gmx-clusterbyfeatures.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

gmx_clusterByFeatures
=====================
It can be used to cluster the conformations of a molecule in a molecular dynamics
trajectory using collection of features. The features could be any quantity as a
function of time such as Projections of eigenvector from PCA or dihedral-PCA,
distances, angles, channel radius etc.

**See details at:** `gmx_clusterByFeatures homepage <https://gmx-clusterbyfeatures.readthedocs.io>`_

.. note:: It is developed for **GROMACS MD trajectory.** However, it can be used with
  any other trajectory format after converting it to GROMACS format trajectory.

Installation on Linux and MacOS
-------------------------------

.. code-block:: bash

    sudo pip3 install gmx-clusterByFeatures

No dependency on GROMACS. Install and use it.

For more details, visit `download and installation <https://gmx-clusterbyfeatures.readthedocs.io/en/latest/install.html>`_ section. 

Usage
-----------

.. list-table:: List of sub-commands available in gmx_clusterByFeatures
    :widths: 1, 4
    :header-rows: 1
    :name: commands-table

    * - Command
      - Function

    * - `cluster <https://gmx-clusterbyfeatures.readthedocs.io/en/latest/commands/cluster.html>`_
      - Main module to perform clustering

    * - `featuresplot <https://gmx-clusterbyfeatures.readthedocs.io/en/latest/commands/featuresplot.html>`_
      - Feature vs Feature plot to check quality of clustering

    * - `distmat <https://gmx-clusterbyfeatures.readthedocs.io/en/latest/commands/distmat.html>`_
      - Distance-matrix related calculations

    * - `matplot <https://gmx-clusterbyfeatures.readthedocs.io/en/latest/commands/matplot.html>`_
      - To visulaize/plot matrix obtained from ``distmat``
      
    * - `hole <https://gmx-clusterbyfeatures.readthedocs.io/en/latest/commands/hole.html>`_
      - To calculate cavity/channel radius using HOLE program
      
    * - `holeplot <https://gmx-clusterbyfeatures.readthedocs.io/en/latest/commands/holeplot.html>`_
      - To calculate average and plot hole output radius file
    
    * - `holefeatures <https://gmx-clusterbyfeatures.readthedocs.io/en/latest/commands/holefeatures.html>`_
      - To write radius as a features for clustering
      
    * - `holeclustersplot <https://gmx-clusterbyfeatures.readthedocs.io/en/latest/commands/holeclustersplot.html>`_
      - To plot or write radius for clusters separately

For more details, visit `usage <https://gmx-clusterbyfeatures.readthedocs.io/en/latest/usage.html>`_ section. 
