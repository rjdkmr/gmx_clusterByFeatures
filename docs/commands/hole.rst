.. |hole| raw:: html

   <a href="http://www.holeprogram.org" target="_blank">HOLE</a>
   
.. |hole-rad-link| raw:: html

   <a href="https://github.com/osmart/hole2/tree/master/rad" target="_blank">link</a>
   
   
``hole``
===========

Description
-----------

It can be used to calculate radius of protein channel/cavity for GROMACS MD
trajectory. It uses |hole| program to calculate radius of cavity/channel
and dumps the output to a text file as a function of tiime. It also extract
channel's outlining residues and dumps to same output file. This output file
can be further read to perform final statistcal operations and plotting.

Please cite the original publication of hole:
  O.S. Smart, J.M. Goodfellow and B.A. Wallace (1993).
  The Pore Dimensions of Gramicidin A. Biophysical Journal 65:2455-2460.
  
Command summary 
----------------

.. code-block:: bash

    gmx_clusterByFeatures hole  [-f [<.xtc/.trr/...>]] [-s [<.tpr/.gro/...>]]
                                [-n [<.ndx>]] [-o [<.dat>]] [-pdb [<.pdb>]] [-b <time>]
                                [-e <time>] [-dt <time>] [-tu <enum>] [-[no]fit]
                                [-endrad <real>] [-sample <real>] [-cvect <vector>]
                                [-cpoint <vector>] [-catmid <int>] [-rad <enum>]

                               
.. list-table:: Options to specify input files to hole
    :widths: 1, 1, 4
    :header-rows: 1
    :name: input-files-table-hole
    :stub-columns: 1
    :align: left

    * - Option
      - Default
      - File type

    * - `-f [\<.xtc/.trr/...\>] <hole.html#f-traj-xtc>`_
      - traj.xtc
      - Trajectory: xtc trr cpt gro g96 pdb tng

    * - `-s [\<.tpr/.gro/...\>] <hole.html#s-topol-tpr>`_
      - topol.tpr
      - Structure+mass(db): tpr gro g96 pdb brk ent

    * - `-n [\<.ndx\>] <hole.html#n-index-ndx>`_
      - index.ndx
      - Index file

      
.. list-table:: Options to specify output files to hole
    :widths: 1, 1, 4
    :header-rows: 1
    :name: output-files-table-hole
    :stub-columns: 1
    :align: left

    * - Option
      - Default
      - File type

    * - `-o [<.dat>] <hole.html#o-radius-dat>`_
      - radius.dat
      - Generic data file containing matrix

    * - `-pdb [<.pdb>] <hole.html#pdb-sphpdb-pdb>`_
      - sphpdb.pdb
      - PDB file for scpheres along the radius

.. list-table:: Other options to hole
    :widths: 10, 12, 40
    :header-rows: 1
    :name: other-options-table-hole
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
      
    * - ``-tu <keyword>``
      - 0
      - Unit for time values: fs, ps, ns, us, ms, s

    * - `-[no]fit <hole.html#fit-nofit>`_
      - Enable
      - Enable fitting and superimposition of the atoms groups.

    * - `-endrad <real> <hole.html#endrad-5>`_
      - 5
      - Radius value (A) after which calculation is stopped.

    * - `-sample <real> <hole.html#sample-0-5>`_
      - 0.5
      - The distance between the planes for the sphere centers.

    * - `-cvect  <vector> <hole.html#cvect-0-0-1>`_
      - 0 0 1
      - Vector along the channel
    
    * - `-cpoint <vector> <hole.html#cpoint-999-999-999>`_
      - 999 999 999
      - Coordinate within a channel as seed for channel/cavity.

    * - `-catmid <int> <hole.html#catmid-1>`_
      - -1
      - Serial number of atom, whoose coordinate acts as seed for channel/cavity.

    * - `-rad <keyword> <hole.html#rad-bondi>`_
      - bondi
      - Radius type for atoms.
      
      
Options to specify input files
--------------------------------

``-f traj.xtc``
~~~~~~~~~~~~~~~~~~~~~~~~
Input trajectory file of ``xtc`` ``trr`` ``cpt`` ``gro`` ``g96`` ``pdb`` or
``tng`` format.


******

``-s topol.tpr``
~~~~~~~~~~~~~~~~~~~~~~~~
An input structure file of ``tpr`` ``gro`` ``g96`` or ``pdb`` format.

******

``-n index.ndx``
~~~~~~~~~~~~~~~~~~~~~~~~~
Two index groups from this file will be prompted for selection. Otherwise,
default index groups will be prompted for selection.

Minimum-distance matrix will be calculated between the two selected atom-groups.

******

Options to specify output files
-------------------------------

``-o radius.dat``
~~~~~~~~~~~~~~~~~~~~~~~~~~
Output file containing channel/cavity axis, radius and outlining residues 
as a function of time. 

******

``-pdb sphpdb.pdb``
~~~~~~~~~~~~~~~~~~~~~~
Output PDB file containing coordinates of sphere inside the channel/Cavity.
Radius of these spehere is channel/cavity radius.
This file can be used to visualize whether obtained radii are from inside
the intended channel/cavity.

******

Other options
-------------

``-fit/-nofit``
~~~~~~~~~~~~~~~~~~~~~~
Weather to superimpose structures on reference structure using least-square fitting.


``-endrad 5``
~~~~~~~~~~~~~~~~~~~~~~
Radius (A) above which the |hole| program regards a result as an indication that the
end of the pore has been reached

``-sample 0.5``
~~~~~~~~~~~~~~~~~~~~~~
The distance (A) between the planes for the sphere centers

``-cvect 0 0 1``
~~~~~~~~~~~~~~~~~~~~~~
This specified a vector which lies in the direction of the channel/cavity.

``-cpoint 999 999 999``
~~~~~~~~~~~~~~~~~~~~~~~~
A point which lies within the channel. If not given, center of mass 
will be used. This point will be used a seed to start calcualtion for
channel/cavity radius.

.. note:: Due to this option, superimposition of structures on reference 
          structure is neccessary.
          
.. note:: Conformations changtes during the simulations, therefore, this 
          coordinate may not be inside the cavity. To dynamically select seed coordinate, 
          use ``-catmid`` option.
          

``-catmid -1``
~~~~~~~~~~~~~~~~~~~~~~~~
Serial number of atom, which lies within the channel and acts
as a seed for channel/cavity. If not given, center of mass will be
used. It can be used to assign seed-coordinate dynamically.

           
``-rad bondi``
~~~~~~~~~~~~~~~~~~~~~~
Radius of atoms considered during channel/cavity calculation.

Accepted categories of radii are:
  * ``bondi`` - For all-atom force-fields, this category can be used.
  * ``amberuni`` - For united-atom force-fields, this category can be used.
  * ``downscaled``
  * ``hardcore``
  * ``simple`` - For united-atom force-fields, this category can be used.
  * ``xplor``
    
These radii are taken from original HOLE implementation. For values of these radii,
please follow this |hole-rad-link|.
