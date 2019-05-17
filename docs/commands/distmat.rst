``distmat``
===========

Description
-----------

This tool can be used to calculate: 

**Average distance matrix:** It can be used to calculate average minimum-distance
matrix of residues between two atom-groups.

**MSF/RMSF in distance-matrix:** It can be used to calculate either variance 
(representing MSF) or standard-deviation (representing RMSF) of distance-matrices.

**Contact map:** It can be used to calculate contact-frequency map over the 
trajectory for the residues that are within a minimum distance given by 
``-ct`` option value.

**Fluctuation in second trajectory with reference to average of first trajectory:**
To calculate  fluctuations (MSF - variance  or RMSF - std. deviation in distance-matrix)
in a trajectory with respect to average distances from another trajectory, use ``-f
traj_for_average.xtc``  and ``-f2 traj_for_rmsf.xtc``. The averages will be calculated
from first trajectory ``traj_for_average.xtc``. Subsequently, variances and deviation will
be calculated for ``traj_for_variance.xtc`` with respect to previosly calculated averages.

**Trajectory and pdb for distance-matrix PCA:**
To speed up the calculation, it uses all available cores of the CPU using
multi-threading. Number of threads/cores could be change by "-nt" option.



Command summary 
----------------

.. code-block:: bash

  gmx_clusterByFeatures distmat [-f [<.xtc/.trr/...>]] [-s [<.tpr/.gro/...>]] [-n [<.ndx>]]
                                [-f2 [<.xtc/.trr/...>]] [-mean [<.dat>]] [-var [<.dat>]]
                                [-std [<.dat>]] [-cmap [<.dat>]] [-pca [<.xtc>]] [-b <time>]
                                [-e <time>] [-dt <time>] [-ct <real>] [-nt <int>] [-gx <int>]
                                [-gy <int>]


                                
.. list-table:: Options to specify input files to distmat
    :widths: 1, 1, 4
    :header-rows: 1
    :name: input-files-table-distmat
    :stub-columns: 1
    :align: left

    * - Option
      - Default
      - File type

    * - `-f [\<.xtc/.trr/...\>] <distmat.html#f-traj-xtc>`_
      - traj.xtc
      - Trajectory: xtc trr cpt gro g96 pdb tng

    * - `-s [\<.tpr/.gro/...\>] <distmat.html#s-topol-tpr>`_
      - topol.tpr
      - Structure+mass(db): tpr gro g96 pdb brk ent

    * - `-n [\<.ndx\>] <distmat.html#n-index-ndx>`_
      - index.ndx
      - Index file

    * - `-f2 [\<.xtc/.trr/...\>] <distmat.html#f2-traj-xtc>`_
      - traj.xtc
      - Trajectory: xtc trr cpt gro g96 pdb tng

.. list-table:: Options to specify output files to distmat
    :widths: 1, 1, 4
    :header-rows: 1
    :name: output-files-table-distmat
    :stub-columns: 1
    :align: left

    * - Option
      - Default
      - File type

    * - `-mean   [<.dat>] <distmat.html#mean-average-dat>`_
      - average.dat
      - Generic data file containing matrix

    * - `-var    [<.dat>] <distmat.html#var-variance-dat>`_
      - variance.dat
      - Generic data file containing matrix

    * - `-std    [<.dat>] <distmat.html#std-stdeviation-dat>`_
      - stdeviation.dat
      - Generic data file containing matrix

    * - `-cmap   [<.dat>] <distmat.html#cmap-contact-map-dat>`_
      - contact_map.dat
      - Generic data file containing matrix

    * - `-pca    [<.xtc>] <distmat.html#pca-pca-xtc>`_
      - pca.xtc
      - Trajectory format file containing distance-matrix of each frame 
      
      
.. list-table:: Other options to distmat
    :widths: 1, 1, 4
    :header-rows: 1
    :name: other-options-table-distmat
    :stub-columns: 1
    :align: left

    * - Option
      - Default
      - Description

    * - ``-b <real>``
      - 0
      - First frame (ps) to read from trajectory

    * - ``-e <real>``
      - 0
      - Last frame (ps) to read from trajectory

    * - ``-dt <real>``
      - 0
      - Only use frame when t MOD dt = first time (ps)

    * - `-ct <real>  <distmat.html#ct-0-4>`_
      - 0.4
      - cut-off distance (nm) for contact map

    * - `-nt <int> <distmat.html#nt-4>`_
      - All CPU cores
      - number of threads for multi-threading

    * - `-gx <int> <distmat.html#gx-5>`_
      - 5
      - Gap between residues along X-axis in distance-matrix for PCA

    * - `-gy <int> <distmat.html#gx-1>`_
      - 1
      - Gap between residues along Y-axis in distance-matrix for PCA

        
Options to specify input files
--------------------------------

``-f traj.xtc``
~~~~~~~~~~~~~~~~~~~~~~~~
Input trajectory file of ``xtc`` ``trr`` ``cpt`` ``gro`` ``g96`` ``pdb`` or
``tng`` format.


******

``-s topol.tpr``
~~~~~~~~~~~~~~~~~~~~~~~~
An input structure file of ``tpr`` ``gro`` ``g96`` or ``pdb`` format. It is **required**
if trajectory is given as input.

******

``-n index.ndx``
~~~~~~~~~~~~~~~~~~~~~~~~~
Two index groups from this file will be prompted for selection. Otherwise,
default index groups will be prompted for selection.

Minimum-distance matrix will be calculated between the two selected atom-groups.

******

``-f2 traj.xtc``
~~~~~~~~~~~~~~~~~~~~~~~~
Input trajectory file of ``xtc`` ``trr`` ``cpt`` ``gro`` ``g96`` ``pdb`` or
``tng`` format.

Second input trajectory. If this trajectory is provided, fluctuations in this trajectory
will be calculated with referece to average-distance matrix of first trajectory.

******

Options to specify output files
-------------------------------

``-mean average.dat``
~~~~~~~~~~~~~~~~~~~~~~~~~~
Output file containing average of minimum-distance matrix. 

******

``-var variance.dat``
~~~~~~~~~~~~~~~~~~~~~~
Output file containing variance of minimum-distance matrix over entire trajectory. 

******

``-std stdeviation.dat``
~~~~~~~~~~~~~~~~~~~~~~~~~~
Output file containing standard-deviation or RMSF of minimum-distance matrix over 
entire trajectory. 

******

``-cmap contact_map.dat``
~~~~~~~~~~~~~~~~~~~~~~~~~
Output file containing contact map over entire trajectory. The contact is determined
using the thershold distance given by ``-ct`` option;

******

``-pca pca.xtc``
~~~~~~~~~~~~~~~~~~~~~~~
Output file containing distance-matrices for each snapshot of the trajectory. This 
file can be used as input to ``gmx covar`` and ``gmx anaeig`` for distance matrix PCA.

A dummy pdb file is also dumped to use with ``gmx covar`` and ``gmx anaeig`` for 
distance matrix PCA.

.. warning:: These two outputs are not real trajectory and pdb file. These two files are
            dumped as a data-container to use with ``gmx covar`` and ``gmx anaeig``.
            For more details, see examples.
            

            
Other options
-------------

``-ct 0.4``
~~~~~~~~~~~~~~~~~~~~~~
cut-off distance (nm) for contact map. Minimum distance below this thershold will be 
considered to be in contact with each other.

******

``-nt 4``
~~~~~~~~~~~~~~~~~~~~~~
Number of parallel threads for distance-matrix computation. 

******

``-gx 5``
~~~~~~~~~~~~~~~~~~~~~~
Gap between residues in distance-matrix along **X-axis** dumped with option `-pca <distmat.html#pca-pca-xtc>`_
for further PCA. This gap reduces the distance-matrix size and subsequently speed-up
the PCA performance.

.. note:: This option **ONLY** affect output from ``-pca`` option.

******

``-gy 1``
~~~~~~~~~~~~~~~~~~~~~~
Gap between residues in distance-matrix along **Y-axis** dumped with option `-pca <distmat.html#pca-pca-xtc>`_
for further PCA. This gap reduces the distance-matrix size and subsequently speed-up
the PCA performance.

.. note:: This option **ONLY** affect output from ``-pca`` option.
