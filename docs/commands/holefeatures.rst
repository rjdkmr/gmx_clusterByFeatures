.. |colormaps| raw:: html

   <a href="https://matplotlib.org/tutorials/colors/colormaps.html#classes-of-colormaps" target="_blank">colormaps list</a>
   
   
   
``holeplot``
=============

Description
-----------

Write channel/cavity radius as features for clustering.

The output file can be used as input features for clustering of channel/cavity 
shape in `cluster <cluster.html>`_.


Command summary 
----------------

.. code-block:: bash

    gmx_clusterByFeatures holefeatures [-h] [-i radius.dat] [-o output.xvg]
                                       [-pca 5] [-xmin XMIN] [-xmax XMAX]
                                       [-endrad ENDRAD] [-ax Z] [-gap 1]
                                       [-b 0] [-e -1] [-do 90]
                                       

Options 
---------

``-i radius.dat``, ``--input radius.dat``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Name of input radius file. Radius file should be obtained from ``hole`` as an 
output file.

******

``-o output.xvg``, ``--output output.xvg``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Name of output file containing radius as function of time at each axis points.
This file can be used as features file for clustering. This file can be
also used to plot radius vs time with external plotting program.

The file name should end with xvg extension, which is recognized by 
"cluster" command.

******

``-pca 5``, ``--pca-pcs 5``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Number of eigenvectors to be considered for the features.
In place for taking radius as features, this option enable PCA of radii
and the resultant projections on eigenvectors can be used as features.

******

                        
``-xmin XMIN``, ``--axis-min XMIN``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Minimum value of axis point after which radius value will be considered for plot.

If not supplied, minimum axis value will be extracted from input radius file.

******

``-xmax XMAX``, ``--axis-max XMAX``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Maximum value of axis point after which radius value will be discarded from plot.

If not supplied, maximum axis value will be extracted from input radius file.

******

``-endrad ENDRAD``, ``--end-radius ENDRAD``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
End/Opening radius.
If radius is larger than this value, radius will not considered 
for average calculation and features output. This option value might be equal or
less than ``-endrad`` value supplied with ``hole`` sub-command.

******

``-ax Z``, ``--axis Z``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Principal axis parallel to the channel or cavity.

******

``-gap 1``, ``--gap 1``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Gap between axis-points in Angstroms
It should be either equal to or larger than ``-sample`` value supplied 
with ``hole`` sub-command.

******

``-b 0``, ``--begin 0``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
First frame in time to read from the input file

******

``-e -1``, ``--end -1``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Last frame in time to read from the input file.
By default ( ``-e -1``), all frames till the end will be read.

******

``-do 90``, ``--data-occupancy 90``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Percentage of radius-data occupancy for axis-points.
If an axis-point has radius-data less than this percentage of frames, 
the axis-point will not be considered for average calculation and 
features output.

This is critical for axis-points, which are at the opening of channel/cavity. 
In several frames, radius-value could be missing and therefore, ``dataOccupancy`` 
threshold could be used to discard those axis points with lots of missing 
radius values over the trajectories.

