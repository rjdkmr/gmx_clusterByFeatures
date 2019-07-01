.. |colormaps| raw:: html

   <a href="https://matplotlib.org/tutorials/colors/colormaps.html#classes-of-colormaps" target="_blank">colormaps list</a>
   
   
   
``featuresplot``
=================

Description
-----------

Features vs Features plot

This can be used to generate plots for features vs features data.
These type of plots are useful to check quality of clustering.

``gmx_clusterByFeatures cluster`` with ``-plot`` option also produces 
features vs features plot. However, the obtained plot is fixed 
and cannot be changed. Therefore, this sub-command can be used 
to obtained plots for desired features with several different 
options to customize the plot.

Command summary 
----------------

.. code-block:: bash

    gmx_clusterByFeatures featuresplot  [-h] [-i radius.dat]
                                        [-feat features.xvg]
                                        [-clid clid.xvg] [-o output.png]
                                        [-b 0] [-e -1] [-tmargin 0.1]
                                        [-lcols 5] [-fs 18] [-wd 8] [-ht 10]
                                        [-dpi 300]

                                  
Options 
---------

``-i input.txt``, ``--input input.txt``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Name of input text file. It should contains two features and their
respective labels in each row. All these values should be separated 
by comma. Each row in file should be in following format:

.. code-block:: bash
    
    [feature no. at X-axis], [feature no. at Y-axis],[X-Label],[Y-Label]


For example, following input will result in four plots:

.. code-block:: bash

    1,2,PC-1,PC-2
    2,3,PC-2,PC-3
    1,3,PC-1,PC-3
    1,4,PC-1,PC-4

    
******

``-feat features.xvg``, ``--features features.xvg``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Input features file.
This file should be same as supplied to ``gmx_clusterByFeatures cluster`` with ``-feat`` option.

******

``-clid clid.xvg``, ``--cluster-id clid.xvg``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Input file containing cluster-id as a function of time.
The number of frames in this file should be same as in features file.

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

.. note:: To list the output formats, use ``gmx_clusterByFeatures holeplot -h``.

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

``-tmargin 0.1``, ``--top-margin 0.1``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Margin at top side of the plot.
If legends overflow into the plot area, margin can be increased to fit the legend.

******

``-lcols 5``, ``--legend-cols 5``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Number of legend columns.
If legend overflow the plot area, legends can be made of more than 
one rows by limiting number of columns to accommodate all legends.

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
