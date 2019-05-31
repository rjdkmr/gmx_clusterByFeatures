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
                                       [-xmin XMIN] [-xmax XMAX]
                                       [-endrad ENDRAD] [-ax Z] [-gap 1]
                                       [-b BEGIN] [-e 1] [-do 90]

Options 
---------

``-i radius.dat``, ``--input radius.dat``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Name of input radius file. Radius file shoudl be obtained from ``hole`` as an 
output file.

******

``-o output.xvg``, ``--output output.xvg``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Name of output file containing radius as function of time at each axis points.
This file can be used as features file for clustersing. This file can be
also used to plot radius vs time with external plotting program.

The file name should end with xvg extension, which is recognized by 
"cluster" command.

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

``-b BEGIN``, ``--begin BEGIN``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
First frame to read from the input file

******

``-e 1``, ``--end 1``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Last frame to read from the input file.
By default ( ``-e -1``), all frames till the end will be read.

******

``-do 90``, ``--data-occupancy 90``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Precentage of radius-data occupancy for axis-points.
If an axis-point has radius-data less than this percentage of frames, 
the axis-point will not be considered for average calculation and 
features output.

This is critical for axis-points, which are at the opening of cahnnel/cavity. 
In several frames, radius-value could be missing and therefore, ``dataOccupancy`` 
thershold could be used to discard those axis points with lots of missing 
radius values over the trajectories.

