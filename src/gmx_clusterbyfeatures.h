#ifndef GMX_CLUSTERBYFEATURES_H
#define GMX_CLUSTERBYFEATURES_H

#include <vector>
#include <map>

#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xtcio.h"

#include "logstream.h"

const char *clusterMetrics[] = { NULL, "prior", "rmsd", "ssr-sst", "pFS", "DBI", NULL };
enum { ePriorClusterMetric = 1, eCRmsdClusterMetric, eSsrSstClusterMetric, ePfsClusterMetric, eDbiClusterMetric };

// Trajectory stuffs, clumped in one data structure for easy handling
struct TrajectoryStuffs {
    const char *filename;
    t_atoms atoms;
    t_trxstatus *status;
    rvec *x;
    real time;
    real dTime = 1;
    int ePBC;
    matrix box;
    int natoms;
    gmx_output_env_t *oenv;
};


//////////////////////////////////////////////////////////////////////////////////////////////
//// ClusteringStuffs Class declaration Start ////////

class ClusteringStuffs {

public:
    // These two will be same for all cluster-stuffs
    static std::vector< std::vector<real> > features; // features
    static std::vector< real > timeInInput; // Real time
    static bool bSortByFeatures;    // If sort by feature is on
    static std::map< int, real > ssrSstRatio, pFS, dbi; //Davies–Bouldin index

    // other variables that are different in case of different cluster-stuffs
    std::vector< int > clidAlongTime; // Clusterid along number of frames
    std::map< int, std::vector< long > > clusterDict;  // Clid -> [idx1, idx2, idx3 ...] indexes corresponds to features, clidAlongTime, timeInClIdInput
    int totalClustNum = 0;
    std::map< int, int > centralStructDict;  // // Clid -> idx, where idx is index of central structure.
    std::vector< int > clusterIndex; // Store cluster-index, order will be used for output files, It does not contain cluster with less than minimum frames
    rvec **centralCoords = NULL;  //Store coordinates of central-structures
    std::vector< std::vector< real > > centrlRmsdMatrix; // RMSD between central structures

    std::map< int, std::vector< real > > avgDistanceToAllDict; // Avg. distance to all other points in trjaectory

    /*
     * Read input features file and store DATA in:
     *     ClusteringStuffs::features
     *     ClusteringStuffs::timeInInput
     *     ClusteringStuffs.clusterIndex
     */
    static int read_features_input(const char *fnDataIn,
                                int minFeatures,
                                LogStream *lstream);

    /*
     * Calculate or extract Cluster-metrics: SSR/SST ratio, Psuedo F-statistics and Davies–Bouldin index.
     * Stored as dictionary in:
     *     ClusteringStuffs::ssrSstRatio
     *     ClusteringStuffs::pFS
     *     ClusteringStuffs::dbi
     */
    static int performClusterMetrics(int eClusterMetrics, int n_clusters, real ssrSstChangeCutoff, LogStream *lstream);

    /*
     * After determining the clusters, this is used to construct:
     *   ClusteringStuffs.clusterDict
     *   ClusteringStuffs.clidAlongTime
     *
     */
    int constructClusterDict(int numMinFrameCluster,
                             LogStream *lstream);

    /*
     * It reads the cluster-id input file. If file also contains features,
     * it also reads these features.
     * Following data are stored:
     *     ClusteringStuffs.clusterDict
     *     ClusteringStuffs.clidAlongTime
     *     ClusteringStuffs.clusterIndex
     *     ClusteringStuffs::timeInInput
     *     ClusteringStuffs::features (When features are present in file.)
     */
    int read_cluster_input(const char *fnDataIn,
                           gmx_bool *bFeatures,
                           int numMinFrameCluster,
                           LogStream *lstream);

    /*
     * Calculates central structures for the clusters.
     * It construct dictionary of {cluster-id -> frame-index} as
     * ClusteringStuffs.centralStructDict
     */
    int calculate_central_struct(LogStream *lstream);

    /*
     * Calculates and returns index of central-structure of an input cluster
     *
     */
    long get_index_central_struct(int clusterID, real **features, long *clusterFrameIndex,
                                  long clusterFrameIndexSize,
                                  int nFeatures);

    /*
     * Extract coordinates and write (optional) PDB file of central structures.
     * The extracted coordinates are stored in ClusteringStuffs.centralCoords
     *
     */
    int write_central_pdbfiles(std::vector < std::string > pdbNames,
                                  int *atomIndex, int atomIndexSize,
                                  TrajectoryStuffs inpTrajStuff);

    /*
     * Calculates RMSD between central structures and stored in
     * ClusteringStuffs.centrlRmsdMatrix
     *
     */
    int rmsd_bw_central_structure(int *fitAtomIndex, int fitAtomIndexSize,
                                  int *rmsdAtomIndex, int rmsdAtomIndexSize,
                                  TrajectoryStuffs inpTrajStuff,
                                  LogStream *lstream);

    /*
     * Convenience method to get central-structures frame-index
     */
    std::vector< long > get_central_ids();

    /*
     * Check whether any of the RMSDs between central structures
     * is less than the input threshold.
     */
    int any_central_rmsd_below_thershold(real thres);

    /*
     * Calculates Davies–Bouldin index and stored in ClusteringStuffs::dbi
     */
    void calculateDaviesBouldinIndex();


};

//// ClusteringStuffs Class declaration END ////////
//////////////////////////////////////////////////////////////////////////////////////////////


// General functions declarations
void CopyRightMsg();

void write_clustered_trajs(const char *fname, ClusteringStuffs *clustStuff,
                           int *atomIndex, int atomIndexSize,
                           TrajectoryStuffs inpTrajStuff,
                           gmx_bool bAlignTrajToCentral,
                           int *fitAtomIndex, int fitAtomIndexSize,
                           LogStream *lstream);

std::vector< std::vector< real > > calculate_rmsd(ClusteringStuffs *clustStuff,
                                                  int *fitAtomIndex, int fitAtomIndexSize,
                                                  int *rmsdAtomIndex, int rmsdAtomIndexSize,
                                                  TrajectoryStuffs inpTrajStuff);

void write_rmsd( std::vector< std::vector< real > > rmsd,
                 const char* fnOutRMSD,
                 std::vector< int > clusterIndex,
                 TrajectoryStuffs inpTrajStuff);


void sort_cluster_frame(std::vector< std::vector< real > > sorter,
                        ClusteringStuffs *clustStuff,
                        std::vector< std::vector< real > > *rmsd );

int gmx_clusterByFeatures(int argc,char *argv[]);



#endif // GMX_CLUSTERBYFEATURES_H

