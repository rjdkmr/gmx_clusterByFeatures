Clustering conformations using distance-matrix PCA
======================================================

In this example, conformational clustering of a flexible protein will be performed using the distance-matrix PCA (dmPCA).
This protein is extremely fleixble and, superposition of conformations are not accurate that is required during the 
conventional PCA. Therefore, to avoid the superposition step, distance-matrix can be used in place of atom-coordinates
for PCA.

Calculation of distance-matrix
------------------------------

At first, distance-matrix over the trajectory can be calculated using `distmat <../commands/distmat.html>`_ command.

.. code-block:: bash

    echo 3 3 | gmx_clusterByFeatures distmat -f input_traj.xtc -s input.tpr -n input.ndx -pca -gx 5
    
Above command produces two outputs:
    * ``pca.xtc``: This file is a container for distance-matrices over the entire trajectory in xtc file format. This 
      is **not a real** trajectory file.
    * ``pca_dummy.pdb``: This is a *dummy* pdb file containing same number of entries as obtained in above xtc file.
    
.. note:: ``-gx 5`` is used to reduce the size of distance-matrix. It means that there is a gap of 4 residues along X-axis
          in distance-matrix. For example, if a protein contains 100 residues, distance-matrix size is 100x100. If ``-gx 5``
          is used, new size is 20x100. 
          
.. note:: ``-gx`` and ``-gy`` options **ONLY** affect output produced with ``-pca`` option of ``distmat``.


Distance-matrix PCA
--------------------

The ``distmat`` produces files ``pca.xtc`` and ``pca_dummy.pdb`` in the above command. These two files are compatible to use with
GROMACS PCA tools. Following steps are used to perform dmPCA.

1. **Covariance, eigenvector and eigenvalue caculcations**

.. code-block:: bash

    echo 0 0 | gmx covar -f pca.xtc -s pca_dummy.pdb -nofit -nomwa -nopbc
    
Above command generated ``eigenvec.trr`` and  ``eigenval.xvg`` files. ``eigenvec.trr`` is neccessary in next command as input.

    
2. **Projections on eigenvectors**

.. code-block:: bash

    echo 0 0 | gmx anaeig -f pca.xtc -s pca_dummy.pdb -first 1 -last 10 -proj
    
In this command, ``-v eigenvec.trr`` was used by default and eigenvectors were read from this file.
A new output file ``proj.xvg`` is generated containing projections on first 10 eigenvectors.
This file is used as an input file in ``gmx_clusterByFeatures``.

Clustering
----------

.. code-block:: bash

    echo 0 3 3 | gmx_clusterByFeatures cluster -s input.tpr -f input_traj.xtc -n input.ndx -feat proj.xvg -method kmeans \
                                               -nfeature 5 -cmetric ssr-sst -ncluster 20 -fit2central -sort rmsdist \
                                               -ssrchange 2 -cpdb clustered-trajs/central.pdb -fout clustered-trajs/cluster.xtc \
                                               -plot pca_cluster.png -rmsd rmsdist/raw.xvg
                                               
K-means clustering was used with maximum number of 20 clusters (``-ncluster 20``). It means, clustering were performed 20 times,
and in each iteration, starting from two, one more cluster was generated. Subsequently, 8 clusters were accepted as final clusters
using change in SSR/SST ratio (``-cmetric ssr-sst`` and ``-ssrchange 2``). RMSD in distance-matrix is used for sorting 
(``-sort rmsdist``) frames in clustered trajectory to avoid the superposition of structure.

.. note:: Check carefully order of index groups selected in the above command.
          
          **a.** First index group - output in central structures and clustered trajectories
          
          **b.** Second index group - rmsd group, here it is C-alpha atoms of protein. **ONLY** used in 
                 in calculation of RMSD matrix, which is dumped in log file with ``-g`` option.
          
          **c.** Third group - Used for superposition by least-square fitting. **ONLY** used in 
                 clustered trajectories to align with central structure.
          
.. note:: In the clustered trajectories, conformations are sorted on the basis of ``rmsdist`` (RMSD in distance-matrix).
          With both ``-sort rmsdist`` and ``-rmsd`` options, rmsdist is calculated for each clustered trajectory.
          
          

Outputs
--------

**Central structures of each cluster:**

.. code-block:: bash

    Cluster-ID  Central Frame   Total Frames 
    1           20876           6715
    2           4902            5803
    3           22646           4958
    4           7717            4721
    5           8287            3137
    6           13989           2791
    7           24749           2090
    8           13740           1801
    

**Output files generated:**

a. ``-g cluster.log`` : log output containing information about the clusters.
b. ``-clid clid.xvg`` : Cluster-id as a function of time.
c. ``-fout clustered-trajs/cluster.xtc`` : 8 clustered trajectories were extracted with name cluster_c{ID}.xtc
d. ``-cpdb clustered-trajs/central.pdb`` : 8 central structures PDB files were extracted with name central_c{ID}.pdb
e. ``-plot pca_cluster.png`` : Plots of feature-vs-feature with different colors as clusters and central structure.
   This plot can be used for visual inspection of clustering. 
