.. |kmeans| raw:: html

   <a href="https://en.wikipedia.org/wiki/K-means_clustering" target="_blank">K-means clustering</a>

.. |DBSCAN| raw:: html

   <a href="https://en.wikipedia.org/wiki/DBSCAN" target="_blank">DBSCAN - Density-based spatial clustering of applications with noise</a>

.. |gmixture| raw:: html

   <a href="https://en.wikipedia.org/wiki/Mixture_model" target="_blank">Gaussian mixture model clustering</a>

.. |elbow| raw:: html

   <a href="https://en.wikipedia.org/wiki/Elbow_method_(clustering)" target="_blank">Elbow method</a>

.. |DBI| raw:: html

  <a href="https://en.wikipedia.org/wiki/Davies%E2%80%93Bouldin_index" target="_blank">DBI : Davies–Bouldin index</a>


gmx_clusterByFeatures
=====================
it can be used to cluster the conformations of a molecule in a molecular dynamics
trajectory using collection of features. The features could be any quantity as a
function of time such as Projections of egienvector from PCA or dihedral-PCA,
distances, angles, channel radius etc.

**See details at:** `gmx_clusterByFeatures homepage <https://gmx-clusterbyfeatures.readthedocs.io>`_

Clustering methods
------------------
Presently three methods are implemented:
  * |kmeans|
  * |DBSCAN|
  * |gmixture|


Clustering metrics
------------------
To determine the number of clustering, following metrics are implemented:
  * RMSD : Root Mean Square deviation between central structures of clusters.
  * SSR/SST ratio ( |elbow| ) : Relative change in SSR/SST ratio in percentage.
  * pFS : Psuedo F-statatics determined from SSR/SST ratio.
  * |DBI|
