#ifndef DO_CLUSTER_H
#define DO_CLUSTER_H

#include "python_helper.h"

class PythonCode {

public:
    static PythonInterpreter* python; //Python class that is wrapper to python interpreter
    static void InitPythonAndLoadFunc(); //Initialize Python
    static std::string importModules();  //Import some modules
    static std::string DoClusteringClass(); // Read DoClustering class defintion

    /*
     * Initialize DoClustering Class as doCluster object in python
     */
    static void initializeClustering(const char* filename, int nPC=2, const char* algo="kmeans", float dbscan_eps=0.5, int dbscan_min_samples=20);

    /*
     * Perform clustering for a given number of cluster
     */
    static void performClustering(int n_clusters);

    /*
     * Get cluster-id for a given number of cluster from Python to C++
     */
    static std::vector< int > getClusterLabels(int n_clusters);

    /*
     * Get SSR/SST ratio and Psuedo F-statistics from Python to C++
     */
    static void getSsrSstStats(int n_clusters, double *ratio, double *pFS);

    /*
     * Plot the features with clusters.
     */
    static void plotFeaturesClusters( int n_clusters, const char *plotfile, std::vector< long > central_id, int fsize=14, float width=12, float height=20);
};

PythonInterpreter* InitPythonAndLoadFunc();

#endif // DO_CLUSTER_H

