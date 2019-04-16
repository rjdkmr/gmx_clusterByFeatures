Clustering conformations using distances between atoms
======================================================

In this example, conformations of G-Quadruplex DNA is clustered according to distances between
three atom-pairs. These atom-pairs form hydrogen bonds during the simulations. However, either 
only one atom-pair among them could form hydrogen bond at a time or niethr can form hydrogen 
bond. The formation of hydrogen bonds are extremely fluctuating between these atom-pairs, and
therefore, clustering will filter conformations based on these hydrogen bonds and distances
between atom-pairs.

Calculation of distances
-------------------------

At first distances between atom-pairs are calculated using ``gmx pairdist`` tool as follows.

.. code-block:: bash

    gmx pairdist -s input.tpr -f input_traj.xtc -ref "resid 1 and atomname N7"  -sel "resid 17 and atomname H62" -o r1N7-r17H62
    gmx pairdist -s input.tpr -f input_traj.xtc -ref "resid 1 and atomname H62" -sel "resid 17 and atomname N3"  -o r1H62-r17N3
    gmx pairdist -s input.tpr -f input_traj.xtc -ref "resid 1 and atomname H61" -sel "resid 17 and atomname N1"  -o r1H61-r17N1
    
In next step, all above files are merged to a single file to use as a `feature input file <../commands/cluster.html#feat-features-xvg>`_.

.. code-block:: bash

    cat r1N7-r17H62.xvg  > distances.xvg
    printf "\n& \n\n"   >> distances.xvg
    cat r1H62-r17N3.xvg >> distances.xvg
    printf "\n& \n\n"   >> distances.xvg
    cat r1H61-r17N1.xvg >> distances.xvg
    printf "\n& \n\n"   >> distances.xvg
    
In next step, clustering is performed.

.. code-block:: bash

    echo 0 1 7 | gmx_clusterByFeatures cluster -s input.tpr -f input_traj.xtc -n input.ndx -feat distances.xvg \
                                               -method kmeans -nfeature 3 -cmetric ssr-sst -ncluster 10 -fit2central \
                                               -sort features -ssrchange 2 -cpdb clustered-trajs/central.pdb \
                                               -fout clustered-trajs/cluster.xtc -plot features_cluster.png \

K-means clustering was used with maximum number of 10 clusters (``-ncluster 10``). It means, clustering will be performed 10 times,
and in each iteration, starting from two, one more cluster was generated. Subsequently, 5 clusters were accepted as final clusters
using change in SSR/SST ratio (``-cmetric ssr-sst`` and ``-ssrchange 2``)

.. note:: Check carefully order of index groups selected in the above command.
          
          **a.** First index group - output in central structures and clustered trajectories
          
          **b.** Second index group - RMSD group, here it is whole G-Quadruplex DNA.
          
          **c.** Third group - Used for superposition by least-square fitting, here it is four tettrads of G-Quadruplex DNA.
          
**Central structures of each cluster:**

.. code-block:: bash

    Cluster-ID  Central Frame   Total Frames 
    1           27707           15640
    2           19260           8435
    3           30851           6338
    4           24630           5332
    5           39369           4894



**RMSD (nm) between central structures:**

.. code-block:: bash

    c1      c2      c3      c4      c5	
    0.000   0.512   0.274   0.266   0.401	
    0.512   0.000   0.483   0.443   0.397	
    0.274   0.483   0.000   0.259   0.484	
    0.266   0.443   0.259   0.000   0.385	
    0.401   0.397   0.484   0.385   0.000	
    
    
**Output files generated:**

a. ``-g cluster.log`` : log output containing information about the clusters.
b. ``-clid clid.xvg`` : Cluster-id as a function of time.
c. ``-fout clustered-trajs/cluster.xtc`` : 5 clustered trajectories were extracted with name cluster_c{ID}.xtc
d. ``-cpdb clustered-trajs/central.pdb`` : 5 central structures PDB files were extracted with name central_c{ID}.pdb
e. ``-plot features_cluster.png`` : Plots of feature-vs-feature with different colors as clusters and central structure.
   This plot can be used for visual inspection of clustering.
   
**Overall Results:**
    
    * Cluster-1: conformations with Hydrogen bonds between A1.N7 and A17.H62 atoms
    * Cluster-2: conformations where distance between all these atom-pairs are between 0.5 to 1.5 nm
    * Cluster-3: conformations with Hydrogen bonds between A1.H61 and A17.N1 atoms
    * Cluster-4: conformations with Hydrogen bonds between A1.H62 and A17.N3 atoms
    * Cluster-5: conformations where distance between all these atom-pairs are between 1.5 to 2.5 nm
    
These results demonstrate that clustering is able to filter out the conformations based on these distances.
