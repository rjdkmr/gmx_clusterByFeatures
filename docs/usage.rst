How to use gmx_clusterByFeatures?
=================================

It contains several sub-commands for different purposes. 

**Other tools are presently in development.**

.. list-table:: List of sub-commands available in gmx_clusterByFeatures
    :widths: 1, 4
    :header-rows: 1
    :name: commands-table

    * - Command
      - Function

    * - `cluster <commands/cluster.html>`_
      - Main module to perform clustering

    * - `featuresplot <commands/featuresplot.html>`_
      - Feature vs Feature plot to check quality of clustering
      
    * - `distmat <commands/distmat.html>`_
      - Distance-matrix related calculations

    * - `matplot <commands/matplot.html>`_
      - To visulaize/plot matrix obtained from ``distmat``
      
    * - `hole <commands/hole.html>`_
      - To calculate cavity/channel radius using HOLE program
      
    * - `holeplot <commands/holeplot.html>`_
      - To calculate average and plot hole output radius file
    
    * - `holefeatures <commands/holefeatures.html>`_
      - To write radius as a features for clustering
      
    * - `holefeatures <commands/holefeatures.html>`_
      - To plot or write radius for clusters separately

sub-commands
------------
.. toctree::
   :maxdepth: 1

   cluster <commands/cluster.rst>
   featuresplot <commands/featuresplot.rst>
   distmat <commands/distmat.rst>
   matplot <commands/matplot.rst>
   hole <commands/hole.rst>
   holeplot <commands/holeplot.rst>
   holefeatures <commands/holefeatures.rst>
   holeclustersplot <commands/holeclustersplot.rst>
