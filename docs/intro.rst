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
During the Molecular Dynamics Simulations, molecule conformations changes considerably
and identifying the conformations is very important to study the biomolecular dynamics.
Conformational clustering can be performed to identify different conformations
sampled during the simulations.

Most widely approach for conformational clustering is to calculate Root Mean Square
Deviations between all conformations and cluster them according to these deviations.
However, for large MD trajectories, this RMSD matrix could be huge and takes very
long time to calculate. Therefore, an alternative method such as features based
clustering can be used to identify the cluster of conformations.

**gmx_clusterByFeatures** can be used to cluster the conformations of a molecule
in a molecular dynamics trajectory using collection of features. The features
could be any quantity as a function of time such as Projections of egienvector
from PCA or dihedral-PCA, distances, angles, channel radius etc.

.. note:: It is developed for **GROMACS MD trajectory**. However, it can be used with
  any other trajectory format after converting it to GROMACS format trajectory.

When Projections of egienvector from PCA or dihedral-PCA is used as features,
it yields clusters depending on the largest conformational changes during the
simulations. Depending on the Clustering metrics, a cluster may contain small
conformational fluctuations around the respective central structure.

When other features such as distances, angles, channel radius etc are used as the
features, the obtained clusters of conformations depends on these features. It can
be used to study the specific conformations given the features while ignoring all
other conformational fluctuations.

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
