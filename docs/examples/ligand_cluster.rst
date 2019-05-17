Clustering ligand conformations using cartesian PCA
=====================================================

In this example, conformation of ligands were clustered with respect to receptor.

At first, PCA was performed using atom-coordinates (cPCA), and subsequently, projections on eigenvectors were used as the features

Atom-coordinates PCA
---------------------

1. **Covariance, eigenvector and eigenvalue caculcations**

.. code-block:: bash
    
    echo 13 14 | gmx covar -s input-files/input.tpr -f input-files/trajectory.xtc -n input-files/input.ndx
    
Here, ``13`` is index group of receptor atoms, which were used for superposition by least-square fitting.
``14`` is index group of ligand without any hydrogen atoms. Above command generated ``eigenvec.trr`` and 
``eigenval.xvg`` files. ``eigenvec.trr`` is neccessary in next command as input.

2. **Projections on eigenvectors**

.. code-block:: bash

    echo 13 14 | gmx anaeig -s input-files/input.tpr -f input-files/trajectory.xtc -n input-files/input.ndx -proj -first 1 -last 20
    
In the above command, ``-v eigenvec.trr`` was used by default and eigenvectors were read from this file.
A new output file ``proj.xvg`` is generated containing projections on first 20 eigenvectors.
This file is used as an input file in ``gmx_clusterByFeatures``.

Clustering
-----------

.. code-block:: bash

    echo 0 14 13 | gmx_clusterByFeatures cluster -s input-files/input.tpr -f input-files/trajectory.xtc -n input-files/input.ndx \
                                                 -feat proj.xvg -method kmeans -nfeature 20 -cmetric ssr-sst -ncluster 15 \
                                                 -fit2central -sort features -cpdb clustered-trajs/central.pdb \
                                                 -fout clustered-trajs/cluster.xtc -plot pca_cluster.png\

K-means clustering was used with maximum number of 15 clusters (``-ncluster 15``). It means, clustering were performed 15 times,
and in each iteration, starting from two, one more cluster was generated. Subsequently, 9 clusters were accepted as final clusters
using change in SSR/SST ratio (``-cmetric ssr-sst`` and ``-ssrchange 2``)

.. note:: Check carefully order of index groups selected in the above command.
          
          **a.** First index group - output in central structures and clustered trajectories
          
          **b.** Second index group - clustering group, here it is ligand without hydrogen atoms
          
          **c.** Third group - Used for superposition by least-square fitting.
          

Outputs
--------

**Central structures of each cluster:**

.. code-block:: bash

    Cluster-ID      Central Frame   Total Frames 
    1               45447           19639
    2               51211           15441
    3               36523           10488
    4               63595           9101
    5               70685           6909
    6               41378           6157
    7               3166            5891
    8               21937           4756
    9               7755            2166



**RMSD (nm) between central structures:**

.. code-block:: bash

    c1      c2      c3      c4      c5      c6      c7      c8      c9
    0.000   0.292   0.701   0.444   0.484   0.498   1.076   0.411   0.883
    0.292   0.000   0.684   0.428   0.418   0.552   1.063   0.439   0.844
    0.701   0.684   0.000   0.834   0.574   0.360   0.860   0.588   0.812
    0.444   0.428   0.834   0.000   0.571   0.705   0.940   0.733   0.763
    0.484   0.418   0.574   0.571   0.000   0.351   0.947   0.670   0.961
    0.498   0.552   0.360   0.705   0.351   0.000   0.959   0.548   0.967
    1.076   1.063   0.860   0.940   0.947   0.959   0.000   1.165   0.614
    0.411   0.439   0.588   0.733   0.670   0.548   1.165   0.000   0.890
    0.883   0.844   0.812   0.763   0.961   0.967   0.614   0.890   0.000
    
    
**Output files generated:**

a. ``-g cluster.log`` : log output containing information about the clusters.
b. ``-clid clid.xvg`` : Cluster-id as a function of time.
c. ``-fout clustered-trajs/cluster.xtc`` : 9 clustered trajectories were extracted with name cluster_c{ID}.xtc
d. ``-cpdb clustered-trajs/central.pdb`` : 9 central structures PDB files were extracted with name central_c{ID}.pdb
e. ``-plot pca_cluster.png`` : Plots of feature-vs-feature with different colors as clusters and central structure.
   This plot can be used for visual inspection of clustering. 

