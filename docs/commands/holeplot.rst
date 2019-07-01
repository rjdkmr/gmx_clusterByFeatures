.. |colormaps| raw:: html

   <a href="https://matplotlib.org/tutorials/colors/colormaps.html#classes-of-colormaps" target="_blank">colormaps list</a>
   
   
   
``holeplot``
=============

Description
-----------

It can be used to generate plots for outputs generated from `hole <hole.html>`_.
It generates plot of average radius with standard deviation as a 
function of axis-points. It also shows the distribution of residues 
that outlines the channel/cavity. 


Command summary 
----------------

.. code-block:: bash

    gmx_clusterByFeatures holeplot  [-h] [-i radius.dat] [-resplot]
                                    [-o output.png] [-vplot]
                                    [-csv output.csv] [-xmin XMIN]
                                    [-xmax XMAX] [-endrad ENDRAD] [-ax Z]
                                    [-gap 1] [-b 0] [-e -1] [-do 90]
                                    [-rfreq 50] [-ymin YMIN] [-ymax YMAX]
                                    [-fs 18] [-rlsize 10] [-wd 6] [-ht 6]
                                    [-dpi 300]

                                  
Options 
---------

``-i radius.dat``, ``--input radius.dat``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Name of input radius file. Radius file shoudl be obtained from ``hole`` as an 
output file.

******

``-resplot``, ``--residues-plot``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Plot distributions of outlining residues to cavity/channel. By default, these 
distributions are not plotted. This option enables the plotting of distributions.

                        
``-o output.png``, ``--output output.png``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Name of the output plot file. The extension will be used to determine the output
format.
                        
Following output formats (system dependent) might be available:
    * ps : Postscript
    * eps : Encapsulated Postscript
    * pdf : Portable Document Format
    * pgf : PGF code for LaTeX
    * png : Portable Network Graphics
    * raw : Raw RGBA bitmap
    * rgba : Raw RGBA bitmap
    * svg : Scalable Vector Graphics
    * svgz : Scalable Vector Graphics
    * jpg : Joint Photographic Experts Group
    * jpeg : Joint Photographic Experts Group
    * tif : Tagged Image File Format
    * tiff : Tagged Image File Format

.. note:: To list the output formats, use ``gmx_clusterByFeatures holeplot -h``.

******

``-csv output.csv``, ``--out-csv output.csv``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Output csv file. The radius as a function of axis-points in csv formatted file. This
file can be read in external data-plotting program.

******

``-vplot``, ``--violinplot``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In place of normal line-plot, it plots radius distribution as violins.
It is useful because this plot gives distribution of radius values over 
entire trajectory for each axis-points

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

******

``-rfreq 50``, ``--residue-frequency 50``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Frequency percentage of residue occurrence during the simulations at a
given axis points. If frequency is less than this threshold, it will 
not considered for plotting. 

******

``-ymin YMIN``, ``--y-axis-min YMIN``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Minimum value at Y-axis. If not supplied minimum value from data will be used. 
It can be useful to minimum and maximum values of Y-axis when several plots 
are compared together.

******

``-ymax YMAX``, ``--y-axis-max YMAX``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Maximum value at Y-axis. If not supplied maximum value from data will be used.
It can be useful to minimum and maximum values of Y-axis when several plots 
are compared together.

******
                        
``-rlsize 10``, ``--rlabel-size 10``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Fontsize of residue label along Y-axis

******

``-fs 14``, ``--font-size 14``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Font-size of all texts in plot

******


``-wd 8``, ``--width 8``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Width of plot in inch

******

``-ht 8``, ``--height 8``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Height of plot in inch

******

``-dpi 300``, ``--dpi 300``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Resolution of plot
