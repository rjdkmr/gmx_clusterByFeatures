.. |colormaps| raw:: html

   <a href="https://matplotlib.org/tutorials/colors/colormaps.html#classes-of-colormaps" target="_blank">colormaps list</a>
   
   
   
``matplot``
===========

Description
-----------

``distmat`` produces several output files containing 
`matrix data <distmat.html#output-files-table>`_. ``matplot`` 
can be used to visualize these data as a 2D map plot.

Command summary 
----------------

.. code-block:: bash

    gmx_clusterByFeatures matplot [-h] [-i distmat.dat] [-o output.png]
                                  [-xs 1] [-ys 1] [-xl Residue]
                                  [-yl Residue] [-cbl nm] [-a auto]
                                  [-cmap binary] [-vmin VMIN] [-vmax VMAX]
                                  [-fs 14] [-cbor vertical] [-wd 8] [-ht 8]
                                  [-dpi 300]
                                  
Options 
---------

``-i distmat.dat``, ``--input distmat.dat``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Input file containing matrix-data. This file is obtained as a output from ``distmat``.

******

``-o output.png``, ``--output output.png``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Name of the output matrix-plot file. The extension will be used to determine the output
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

.. note:: To list the output formats, use ``gmx_clusterByFeatures matplot -h``.

******

``-xs 1``, ``--x-start 1``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
First residue number along X-axis. Input file does not contain information about 
residue number. Therefore, this option can be used to set the number of residues 
along X-axis.

******

``-ys 1``, ``--y-start 1``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
First residue number along Y-axis. Input file does not contain information about 
residue number. Therefore, this option can be used to set the number of residues 
along Y-axis.

******

``-xl Residue``, ``--x-label Residue``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
X-axis label

******

``-yl Residue``, ``--y-label Residue``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Y-axis label

******

``-cbl (nm)``, ``--colorbar-label (nm)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Label for color bar

******

``-a auto``, ``--image-aspect auto``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Controls the aspect ratio of the axes.
The aspect is of particular relevance for images since it may distort 
the image, i.e. pixel will not be square.

Following two options are available:
    * ``-a equal`` : Ensures an aspect ratio of 1. Pixels will be square.
    * ``-a auto``  : The axes is kept fixed and the aspect is adjusted so
                     that the data fit in the axes. In general, this will 
                     result in non-square pixels.


******

``-cmap binary``, ``--colormap binary``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Name of colormap by which matrix image will be colored.
To preview the available colormaps, visit |colormaps|.


Following colormaps might be available:
::

    Accent           Blues            BrBG             BuGn             
    BuPu             CMRmap           Dark2            GnBu             Greens           
    Greys            OrRd             Oranges          PRGn             Paired           
    Pastel1          Pastel2          PiYG             PuBu             PuBuGn           
    PuOr             PuRd             Purples          RdBu             RdGy             
    RdPu             RdYlBu           RdYlGn           Reds             Set1             
    Set2             Set3             Spectral         Wistia           YlGn             
    YlGnBu           YlOrBr           YlOrRd           afmhot           autumn           
    binary           bone             brg              bwr              cividis          
    cool             coolwarm         copper           cubehelix        flag             
    gist_earth       gist_gray        gist_heat        gist_ncar        gist_rainbow     
    gist_stern       gist_yarg        gnuplot          gnuplot2         gray             
    hot              hsv              inferno          jet              magma            
    nipy_spectral    ocean            pink             plasma           prism            
    rainbow          seismic          spring           summer           tab10            
    tab20            tab20b           tab20c           terrain          viridis          
    winter           

**Reverse of the available colormaps** are also available with same name suffixed by 
"_r". For example, reverse of binary colormap is binary_r, reverse of gist_earth
is gist_earth_r etc.

.. note:: To list all available colormaps, use ``gmx_clusterByFeatures matplot -h``.

******

``-vmin VMIN``, ``--min-value VMIN``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Minimum value to begin color-mapping.
If not provided, minimum value of whole matrix will be considered.

******

``-vmax VMAX``, ``--max-value VMAX``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Maximum value to end color-mapping.
If not provided, maximum value of whole matrix will be considered.

******

``-fs 14``, ``--font-size 14``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Font-size of all texts in plot

******

``-cbor vertical``, ``--colorbar-orientation vertical``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Orientation of color bar
Following keywords are available:
* ``-cbor vertical`` - vertical colorbar at right side
* ``-cbor horizontal`` - horizontal colorbar at top

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
