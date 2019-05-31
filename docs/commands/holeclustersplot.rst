.. |colormaps| raw:: html

   <a href="https://matplotlib.org/tutorials/colors/colormaps.html#classes-of-colormaps" target="_blank">colormaps list</a>
   
   
   
``holeclustersplot``
=============

Description
-----------

It can be used to plot radius of cavity/channel for clusters seperately.
It reads radius file from `hole <hole.html>`_ and cluster-id file from
`cluster <cluster.html>`_, and extract radius of each cluster separately 
and plot them in one plot. This plot could be extremely
useful to compare radius along the channel/cavity in all clusters.


Command summary 
----------------

.. code-block:: bash

    gmx_clusterByFeatures holeclustersplot [-h] [-i radius.dat]                                                                                                                                                                           
                                           [-clid clid.xvg] [-o output.png]
                                           [-csv output.csv] [-xmin XMIN]
                                           [-xmax XMAX] [-endrad ENDRAD]
                                           [-ax Z] [-gap 1] [-b BEGIN]
                                           [-e -1] [-do 90] [-stdbar]
                                           [-dl 0] [-rmargin 0.15]
                                           [-lcols 1] [-fs 18] [-wd 6]
                                           [-ht 6] [-dpi 300]
                                  

Options 
---------

``-i radius.dat``, ``--input radius.dat``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Name of input radius file. Radius file shoudl be obtained from ``hole`` as an 
output file.

******

``-clid clid.xvg``, ``--clid clid.xvg``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Input file containing cluster-id as a function of time.
The number of frames in this file should be same as in input radius file.

******
                        
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

.. note:: To list the output formats, use ``gmx_clusterByFeatures holeclustersplot -h``.

******

``-csv output.csv``, ``--out-csv output.csv``
Output csv file.
The radius as a function of axis-points in csv formatted file. This
file can be read in external data-plotting program.

                        
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

******

``-stdbar``, ``--stddev-bar``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To show standard deviation as error-bar
If it is supplied, standard deviation will be shown as an error-bar in the plot.

******

``-dl 0``, ``--discard-lasts 0``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Number of smallest clusters to discard from the plotting.
It can be useful to filter out few smallest clusters because these may 
contain small number of frames.

******

``-rmargin 0.15``, ``--right-margin 0.15``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Margin at right side of the plots.
If legends overflow into the plot area, margin can be increased to fit the legend.

******

``-lcols 1``, ``--legend-cols 1``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Number of legend columns
If legend overflow into the plot area, legends can be made of more than 
one column to accomodate all legends.


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
