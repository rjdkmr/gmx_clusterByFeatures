/*
 * This file is part of gmx_clusterByFeatures
 *
 * Author: Rajendra Kumar
 * Copyright (C) 2018  Rajendra Kumar
 *
 * gmx_clusterByFeatures is a free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * gmx_clusterByFeatures is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with gmx_clusterByFeatures.  If not, see <http://www.gnu.org/licenses/>.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 */

#include <cstdio>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>
#include <map>
#include <utility>
#include <numeric>
#include <algorithm>


#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/cmdlineinit.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/index.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/fileio/oenv.h"

#include "gmx_clusterbyfeatures.h"
#include "logstream.h"
#include "do_cluster.h"


template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
    	if (! item.empty() )
    		*(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

template <typename T>
std::vector<size_t> argsort(std::vector<T> const& values) {
    std::vector<size_t> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<size_t>(0));

    std::sort(
        begin(indices), end(indices),
        [&](size_t a, size_t b) { return values[a] < values[b]; }
    );
    return indices;
}


std::vector< int > getSortedKeys(std::map< int, std::vector< long > > dict) {
    std::vector< int > keys;
    for (std::map< int, std::vector< long > >::iterator it= dict.begin(); it!= dict.end(); ++it) {
        keys.push_back(it->first);
    }
    std::sort(keys.begin(), keys.end());
    return keys;
}

template <class T>
T **convert2DVectorToArray(std::vector< std::vector< T > > *Vector, int N){
    T **array;
    snew(array, N);
    for(int i=0; (i < N); i++)  {
        array[i] = &Vector->at(i)[0];
    }
    return array;
}

real calculate_distance(real *v1, real *v2, int n, bool square)   {
    int i=0;
    real distance2 = 0;
    real distance = 0;

    for(i=0; i<n; i++)  {
        distance2 += ( (v1[i] - v2[i]) * (v1[i] - v2[i]) );
    }

    distance = std::sqrt(distance2);

    if(square)
        return distance2;
    else
        return distance;
}

/*
 * get the vector of filenames according to the cluster index and a prefix.
 */
std::vector < std::string > get_outFile_names(const char *fname, std::vector< int > clusterIndex, std::string ext, std::string prefix) {
    std::vector < std::string > outNames, temp;
    std::string filename(fname), basename;
    temp = split(filename, '.');
    basename = temp[temp.size()-2];

    for (size_t i=0; i<clusterIndex.size(); i++)  {
        if (! prefix.empty() )
            outNames.push_back(basename + "_" + prefix +"_c" + std::to_string(clusterIndex[i]) + ext);
        else
            outNames.push_back(basename + "_c" + std::to_string(clusterIndex[i]) + ext);

    }

    /*
    for (int i=0; i<clustNum; i++)  {
        std::cout<<outNames[i]<<"\n";
    }*/

    return outNames;
}


void set_dTime(TrajectoryStuffs *inpTrajStuff)   {
    real oldTime, dTime;

    // Determine dt from original trajectory
    rewind_trj(inpTrajStuff->status);
    oldTime = inpTrajStuff->time;
    do {
        dTime = inpTrajStuff->time - oldTime;
        if (nframes_read(inpTrajStuff->status)>3)
            break;
        oldTime = inpTrajStuff->time;
        // std::cout<<"\n\n Input Trajectory dt = "<<inpTrajStuff->time<<"\n\n";
    } while( read_next_x(inpTrajStuff->oenv, inpTrajStuff->status, &inpTrajStuff->time, inpTrajStuff->x, inpTrajStuff->box) );
    rewind_trj(inpTrajStuff->status);

    inpTrajStuff->dTime = dTime;

    std::cout<<"\n\n Input Trajectory dt = "<<dTime*output_env_get_time_factor(inpTrajStuff->oenv)<<" "<<output_env_get_time_unit(inpTrajStuff->oenv)<<"\n\n";
}

rvec * copy_rvec_coord(rvec *inpRvec, int natoms) {

    rvec *outRvec;
    snew(outRvec, natoms);
    for (int i=0; i < natoms; i++) {
        outRvec[i][XX] = inpRvec[i][XX];
        outRvec[i][YY] = inpRvec[i][YY];
        outRvec[i][ZZ] = inpRvec[i][ZZ];
    }

    return outRvec;
}


// ############################### ClusteringStuffs ##########################


// Initialization of static variables
std::vector< std::vector<real> > ClusteringStuffs::features;
std::vector< real > ClusteringStuffs::timeInInput;
bool ClusteringStuffs::bSortByFeatures = false;
std::map< int, real > ClusteringStuffs::ssrSstRatio;
std::map< int, real > ClusteringStuffs::pFS;
std::map< int, real > ClusteringStuffs::dbi;


int ClusteringStuffs::write_central_pdbfiles(std::vector < std::string > pdbNames,
                                              int *outAtomIndex, int outAtomIndexSize,
                                              TrajectoryStuffs inpTrajStuff)  {
    std::string title;
    std::vector< int > sortedClusterIds = getSortedKeys(this->clusterDict);
    int ftp = fn2ftp(inpTrajStuff.filename);
    t_fileio *fio = trx_get_fileio(inpTrajStuff.status);
    gmx_bool bRet;
    rvec **centralCoords;
    bool bWriteFile = false;

    // If a null pointer is provided, only return central coords, no output file
    if (pdbNames.empty())  {
        std::cout<<"\n\nExtracting coordinates of the central structure...\n";
    }
    else {
        bWriteFile = true;
        std::cout<<"\n\nWriting central structure to pdb-files...\n";
    }

    snew(centralCoords, sortedClusterIds.size());

    if (ftp != efXTC)   {
        do {
            for (size_t c=0; c < sortedClusterIds.size(); c++) {
                if (inpTrajStuff.time == this->timeInInput[this->centralStructDict[sortedClusterIds[c]]] )   {
                    if(bWriteFile) {
                        title = "Cluster - " + std::to_string(sortedClusterIds[c]) + "; Time = " + std::to_string(inpTrajStuff.time);
                        //std::cout<<inpTrajStuff.time<<" "<<ClusteringStuffs::timeInClIdInput[this->centralStructDict[sortedClusterIds[c]]]<<"\n";
                        write_sto_conf_indexed(pdbNames[c].c_str(), title.c_str(), &inpTrajStuff.atoms, inpTrajStuff.x, \
                                               NULL, inpTrajStuff.ePBC, inpTrajStuff.box, outAtomIndexSize, outAtomIndex);
                    }
                    centralCoords[c] = copy_rvec_coord(inpTrajStuff.x, inpTrajStuff.natoms);
                }
            }

        } while(read_next_x(inpTrajStuff.oenv, inpTrajStuff.status, &inpTrajStuff.time, inpTrajStuff.x, inpTrajStuff.box));
    }
    else {
        for (size_t c=0; c < sortedClusterIds.size(); c++) {
            bRet = xtc_seek_time(fio, this->timeInInput[this->centralStructDict[sortedClusterIds[c]]], inpTrajStuff.natoms, FALSE);
            if (bRet == -1) {
                gmx_fatal(FARGS, "Frame for this time is not found in trajectory");
            }
            read_next_x(inpTrajStuff.oenv, inpTrajStuff.status, &inpTrajStuff.time, inpTrajStuff.x, inpTrajStuff.box);

            if(bWriteFile) {
                title = "Cluster - " + std::to_string(sortedClusterIds[c]) + "; Time = " + std::to_string(inpTrajStuff.time);
                write_sto_conf_indexed(pdbNames[c].c_str(), title.c_str(), &inpTrajStuff.atoms, inpTrajStuff.x, \
                                   NULL, inpTrajStuff.ePBC, inpTrajStuff.box, outAtomIndexSize, outAtomIndex);
            }
            centralCoords[c] = copy_rvec_coord(inpTrajStuff.x, inpTrajStuff.natoms);
        }
    }


    this->centralCoords = centralCoords;

    return TRUE;
}


int ClusteringStuffs::rmsd_bw_central_structure( int *fitAtomIndex, int fitAtomIndexSize,
                                                 int *rmsdAtomIndex, int rmsdAtomIndexSize,
                                                 TrajectoryStuffs inpTrajStuff,
                                                 LogStream *lstream)    {

    real *w_rls, *w_rms,  temp;
    int i = 0;
    rvec x_shift_ref, x_shift_i;
    // int cluster;
    // int clusterLength
    std::vector< int > sortedClusterIds = getSortedKeys(this->clusterDict);

    std::vector< std::vector< real > > centrlRmsdMatrix( sortedClusterIds.size(),  std::vector< real >(sortedClusterIds.size(), 0));

    std::cout<<"\n\nCalculating RMSD between central structures...\n";

    // Weight factor initialization for fitting and RMSD calculation
    snew(w_rms, inpTrajStuff.atoms.nr);
    snew(w_rls, inpTrajStuff.atoms.nr);

    if (fitAtomIndexSize < 3)
        gmx_fatal(FARGS, "Need >= 3 points to fit!\n" );

    // Assign weight-factor for fitting
    for (i = 0; i < fitAtomIndexSize; i++) {
        if (inpTrajStuff.atoms.atom[fitAtomIndex[i]].m != 0)
            w_rls[fitAtomIndex[i]] = inpTrajStuff.atoms.atom[fitAtomIndex[i]].m;
        else
            w_rls[fitAtomIndex[i]] = 1;
    }

    // Assign weight-factor for RMSD calculation
    for (i = 0; i < rmsdAtomIndexSize; i++)
            w_rms[rmsdAtomIndex[i]] = 1;


    //Loop over cluster start here
    for (size_t cref=0; cref < sortedClusterIds.size(); cref++) {

        // cluster = sortedClusterIds[cref];
        // clusterLength = this->clusterDict[cluster].size();


        // Reset to origin and store the translation factor of reference coordinate
        copy_rvec(this->centralCoords[cref][0], x_shift_ref);
        reset_x(fitAtomIndexSize, fitAtomIndex, inpTrajStuff.atoms.nr, NULL, this->centralCoords[cref], w_rls);
        rvec_dec(x_shift_ref, this->centralCoords[cref][0]);

        for (size_t ci=0; ci < cref; ci++) {
            copy_rvec(this->centralCoords[ci][0], x_shift_i);
            reset_x(fitAtomIndexSize, fitAtomIndex, inpTrajStuff.atoms.nr, NULL, this->centralCoords[ci], w_rls);
            rvec_dec(x_shift_i, this->centralCoords[ci][0]);


            do_fit(inpTrajStuff.natoms, w_rls, this->centralCoords[cref], this->centralCoords[ci]);
            temp = calc_similar_ind(FALSE, rmsdAtomIndexSize, rmsdAtomIndex, w_rms, this->centralCoords[ci], this->centralCoords[cref]);
            centrlRmsdMatrix.at(cref).at(ci) = temp;
            centrlRmsdMatrix.at(ci).at(cref) = temp;

            // Translate the central structure to original position
            for (i = 0; i < inpTrajStuff.natoms; i++)
                rvec_inc(this->centralCoords[ci][i], x_shift_i);

        }

        // Translate the central structure to original position
        for (i = 0; i < inpTrajStuff.natoms; i++)
            rvec_inc(this->centralCoords[cref][i], x_shift_ref);
    }
    //Loop over cluster end here

    sfree(w_rms);
    sfree(w_rls);

    lstream->setprecision(3);
    *lstream<<"\n\n=====================================\n";
    *lstream<<" Central structurs - RMSD matrix \n";
    *lstream<<"=====================================\n";
    for (size_t cref=0; cref < sortedClusterIds.size(); cref++)
        *lstream<<"c"<<sortedClusterIds[cref]<<"\t";
    *lstream<<"\n";
    for (size_t cref=0; cref < sortedClusterIds.size(); cref++) {
        for (size_t ci=0; ci < sortedClusterIds.size(); ci++) {
            *lstream<<centrlRmsdMatrix[cref][ci]<<"\t";
        }
        *lstream<<"\n";
    }
    *lstream<<"\n\n=====================================\n";
    lstream->resetprecision();

    this->centrlRmsdMatrix = centrlRmsdMatrix;

    return TRUE;
}

long ClusteringStuffs::get_index_central_struct(int clusterID, real **features, long *clusterFrameIndex,
                                                unsigned long clusterFrameIndexSize, int nFeatures)   {
    real dist = 0, dist_prev=9999999, avg_dist=0, prev_avg_dist=0;
    long idx = 0, prev_idx=0;
    std::vector < real > distAll(clusterFrameIndexSize, 0.0);
    std::vector < real > avgPoint(nFeatures, 0.0);
    int same_idx_counter=0;
    unsigned long start = 1, for_limit = 50;


    // Get the average point
    for(int i=0; i < nFeatures; i++ )  {
        for(unsigned long j=0; j < clusterFrameIndexSize; j++ )  {
            avgPoint[i] += features[clusterFrameIndex[j]][i];
        }
        avgPoint[i] = avgPoint[i]/clusterFrameIndexSize;
    }

    // determine a point which is closest to average point
    for(unsigned long i=0; i < clusterFrameIndexSize; i++ )  {
        dist = calculate_distance(features[clusterFrameIndex[i]], avgPoint.data(), nFeatures, true);
        if( dist < dist_prev) {
            idx = clusterFrameIndex[i];
            dist_prev = dist;
        }
    }
    prev_idx = idx;


    // Now determine the central point with iterative process
    while(1) {
        std::vector<size_t> sortedIndex;
        avg_dist = 0;

        // Calculate distance to all other points from current central point
        for(unsigned long i=0; i < clusterFrameIndexSize; i++ )  {
            dist = calculate_distance(features[clusterFrameIndex[i]], features[idx], nFeatures, false);
            distAll[i] = dist;
            avg_dist += dist;
        }

        // Get average distance for current central point
        avg_dist = avg_dist/(clusterFrameIndexSize-1);
        prev_avg_dist = avg_dist;

        // Sort the all distances to current central point to get its nearest neighbour points
        sortedIndex = argsort(distAll);


        // Now check which nearest neighbour point has less average distance value as compared with current central point
        for(unsigned long i=start; i < sortedIndex.size(); i++ ) {
            avg_dist = 0;
            for(unsigned long j=0; j < clusterFrameIndexSize; j++ )  {
                dist = calculate_distance(features[clusterFrameIndex[j]], features[clusterFrameIndex[sortedIndex[i]]], nFeatures, false);
                avg_dist += dist;
            }
            avg_dist = avg_dist/(clusterFrameIndexSize-1);

            // If a point is found, assigned this point as current central point
            if( avg_dist < prev_avg_dist) {
                idx = clusterFrameIndex[sortedIndex[i]];
                prev_avg_dist = avg_dist;
                break;
            }

            // Only check first 50 nighbours in first iteration
            // if central point not found, check next 50 and so on in subsequent iterations
            if ( (i >= for_limit) || (i > clusterFrameIndexSize) )
                break;
        }

        // std::cout<<same_idx_counter<<" "<<start<<" "<<idx<<" "<<prev_idx<<"\n";

        // Check if new central point is similar to that of previous iteration
        if (prev_idx == idx) {               // If similar, then check for new central point among next 50 nighbours in above step
            same_idx_counter += 1;
            start = for_limit * same_idx_counter;
        }
        else {                               // If not similar, then reset counters
            prev_idx = idx;
            same_idx_counter = 0;
            start = 1;
        }

        if (start > clusterFrameIndexSize)
            break;

        // In case of 20 iterations when no new central points found, stop the process
        // and current central point is considered as final central point
        if (same_idx_counter >= 20)
            break;

    }

    // Store distances for later use
    if (this->bSortByFeatures) {
        for(unsigned long j=0; j < clusterFrameIndexSize; j++ )  {
            distAll[j] = 0.0;
            dist = calculate_distance(features[clusterFrameIndex[j]], features[idx], nFeatures, false);
            distAll[j] = dist;
        }
        this->avgDistanceToAllDict.emplace(clusterID, distAll);
    }
    return idx;
}


int ClusteringStuffs::calculate_central_struct(LogStream *lstream){
    std::vector< long > clusterFrameIndex;
    std::vector< int > sortedKeys = getSortedKeys(this->clusterDict);
    real **featuresArray=NULL;
    long minIdx;
    int nFeatures = this->features.at(0).size();

    // Convert to array for speed-up
    featuresArray = convert2DVectorToArray(&this->features, this->features.size());

    printf("\n");
    for (size_t i=0; i<sortedKeys.size(); i++)  {
        printf("\rCalculating central structure for cluster-%d ...",sortedKeys[i]) ;
        fflush(stdout);
        clusterFrameIndex = this->clusterDict.at(sortedKeys[i]);
        minIdx = this->get_index_central_struct(sortedKeys[i], featuresArray, &clusterFrameIndex[0], clusterFrameIndex.size(), nFeatures);
        //std::cout<<it->first<<"\t\t"<<minIdx<<"\t\t"<<clusterFrameIndex.size()<<"\n";
        this->centralStructDict.emplace(sortedKeys[i], minIdx);
    }
    printf("\n");

    *lstream<<"\n===========================================";
    *lstream<<"\nCluster-ID\tCentral Frame\tTotal Frames \n";
    for (size_t i=0; i<sortedKeys.size(); i++)
        *lstream<<sortedKeys[i]<<"\t\t"<<this->centralStructDict[sortedKeys[i]]<<"\t\t"<<this->clusterDict.at(sortedKeys[i]).size()<<"\n";
    *lstream<<"===========================================\n\n";


    sfree(featuresArray);  //Free array

    return TRUE;
}

std::vector< long > ClusteringStuffs::get_central_ids(){
    std::vector< int > sortedKeys = getSortedKeys(this->clusterDict);
    std::vector< long > central_id;
    //Loop over cluster start here
    for (size_t i=0; i<sortedKeys.size(); i++)  {
        central_id.push_back( this->centralStructDict.at(sortedKeys[i]) );
    }
    return central_id;
}

void ClusteringStuffs::calculateDaviesBouldinIndex() {
    std::vector< int > sortedKeys = getSortedKeys(this->clusterDict);
    int n_cluster = sortedKeys.size();

    if (n_cluster == 1){
        ClusteringStuffs::dbi.emplace(n_cluster, 0);
        return;
    }

    std::vector< long > clusterFrameIndex;
    int nFeatures = this->features.at(0).size();
    std::vector< std::vector < real > > centroids(n_cluster, std::vector< real > (nFeatures, 0.0));
    std::vector< real > variances(n_cluster, 0.0);
    real sums = 0, di_sum = 0, dbi = 0, rij = 0;

    // Calculate Ai and Si
    for(size_t c=0; c < sortedKeys.size(); c++)    {
        clusterFrameIndex = this->clusterDict.at(sortedKeys[c]);

        // Get the centeroids (Ai)
        for(int i=0; i < nFeatures; i++ )  {
            for(unsigned long j=0; j < clusterFrameIndex.size(); j++ )  {
                centroids[c][i] += features[clusterFrameIndex[j]][i];
            }
            centroids[c][i] = centroids[c][i]/clusterFrameIndex.size();
        }

        // Calculate variances (Si) around the centroids in each cluster
        sums = 0;
        for(unsigned long j=0; j < clusterFrameIndex.size(); j++ )  {
            sums += calculate_distance(features[clusterFrameIndex[j]].data(), centroids[c].data(), nFeatures, false);
        }
        variances[c] = sums/clusterFrameIndex.size();
    }

    // calculate Di
    for(int i=0; i < n_cluster; i++ )  {
        real max_rij = 0;
        for(long j=0; j < n_cluster; j++ )  {
            if (i==j)
                continue;

            // Calculate Rij =  (Si + Sj)/Mij, where Mi is distance between two centroids
            rij = (variances[i] + variances[j]) / calculate_distance(centroids[i].data(), centroids[j].data(), nFeatures, false);


            // Determine max Rij
            if(rij > max_rij)
                max_rij = rij;
        }

        // Di = Max{Rij}, here added to calculate average later
        di_sum += max_rij;
    }

    // dbi = Di/C
    dbi = di_sum / n_cluster;
    ClusteringStuffs::dbi.emplace(n_cluster, dbi);

    // std::cout<<"\n ###### ncluster, DBI = "<<n_cluster<<" "<<dbi<<"\n";

    return;
}

int ClusteringStuffs::constructClusterDict(int numMinFrameCluster,
                                           LogStream *lstream) {

    std::vector<long> tempIndexVector;   // store a value for making first key-value pair
    int tmpClid;
    long index;

    for(unsigned long i=0; i < this->clidAlongTime.size(); i++) {
        tmpClid = this->clidAlongTime[i];

        // Ignore cluster with -1, particulalry genrated in DBSCAN
        if (tmpClid == -1) {
            continue;
        }

        index = i;

        // Make and expand clid->index dictionary
        tempIndexVector.clear();
        tempIndexVector.shrink_to_fit();
        if ( this->clusterDict.empty() ) {
            tempIndexVector.push_back(index);
            this->clusterDict.emplace( tmpClid, tempIndexVector );
        }
        else    {
            if (this->clusterDict.count(tmpClid) == 0)   {
                tempIndexVector.push_back(index);
                this->clusterDict.emplace(tmpClid, tempIndexVector);
            }
            else  {
                tempIndexVector = this->clusterDict.at(tmpClid);
                this->clusterDict.erase(tmpClid);
                tempIndexVector.push_back(index);
                this->clusterDict.emplace(tmpClid, tempIndexVector);
            }
        }
    }


    // Store cluster-index, order will be used for output files
    *lstream<<"\n===========================\nCluster-ID\tTotalFrames\n";
    for (std::map< int, std::vector<long> >::iterator it=this->clusterDict.begin(); it!=this->clusterDict.end(); ++it)  {

        if (it->second.size() < (unsigned int)numMinFrameCluster) {
            continue;
        }

        *lstream<<it->first<<"\t\t"<<it->second.size()<<"\n";
        this->clusterIndex.push_back(it->first);
    }
    *lstream<<"===========================\n\n";

    this->totalClustNum = this->clusterDict.size();

    return TRUE;
}


int ClusteringStuffs::read_cluster_input(const char *fnDataIn,
                                         gmx_bool *bFeatures,
                                         int numMinFrameCluster,
                                         gmx_output_env_t *oenv,
                                         LogStream *lstream) {

    std::ifstream fpDataIn;
    std::string line;
    std::vector< std::string > temp;
    real time, tempFeature;
    std::vector<long> tempIndexVector;   // store a value for making first key-value pair
    std::vector<real> tempFeatureVector;   // store a value for making first key-value pair
    int tmpClid;

    long index = 0;

    // Open data file, Do it in C++ way for easy reading and parsing
    fpDataIn.open(fnDataIn, std::ifstream::in);
    while(1)	{

        std::getline(fpDataIn, line); // Read each line of the file

        // Stop reading when reached end of the file
        if (fpDataIn.eof())
            break;

        // Split each line with a whitespace
        temp = split(line, ' ');

        // If blank line skip
        if ( temp.empty() )
            continue;

        // Skip line start with # and @, mostly in xvg file
        if ( (temp[0][0] == '#') || (temp[0][0] == '@') )
            continue;

        // Check if features are present
        if (temp.size() > 2)
            *bFeatures = TRUE;

        // Get current time and clid
        // Also convert time to ps (default unit in trajectory file)
        time = (real)std::stod(temp[0]) * output_env_get_time_invfactor(oenv); 
        tmpClid = std::stoi(temp[1]);

        // Store current time and clid in the variable
        timeInInput.push_back(time);
        this->clidAlongTime.push_back(tmpClid);

        // Make and expand clid->index dictionary
        tempIndexVector.clear();
        tempIndexVector.shrink_to_fit();
        if ( this->clusterDict.empty() ) {
            tempIndexVector.push_back(index);
            this->clusterDict.emplace( tmpClid, tempIndexVector );
        }
        else    {
            if (this->clusterDict.count(tmpClid) == 0)   {
                tempIndexVector.push_back(index);
                this->clusterDict.emplace(tmpClid, tempIndexVector);
            }
            else  {
                tempIndexVector = this->clusterDict.at(tmpClid);
                this->clusterDict.erase(tmpClid);
                tempIndexVector.push_back(index);
                this->clusterDict.emplace(tmpClid, tempIndexVector);
            }
        }

        // Read Feature values here if present
        if (*bFeatures) {
            tempFeatureVector.clear();
            tempFeatureVector.shrink_to_fit();
            for (size_t i=2; i < temp.size(); i++)  {
                tempFeature = (real)std::stod(temp[i]);
                tempFeatureVector.push_back( tempFeature );
            }
            this->features.push_back(tempFeatureVector);
        }
        index += 1;
    }

    // Store cluster-index, order will be used for output files
    *lstream<<"\n===========================\nCluster-ID\tTotalFrames\n";
    for (std::map< int, std::vector<long> >::iterator it=this->clusterDict.begin(); it!=this->clusterDict.end(); ++it)  {

        if (it->second.size() < (size_t)numMinFrameCluster)
            continue;

        *lstream<<it->first<<"\t\t"<<it->second.size()<<"\n";
        this->clusterIndex.push_back(it->first);
    }
    *lstream<<"===========================\n\n";

    this->totalClustNum = this->clusterDict.size();

    return TRUE;
}

int ClusteringStuffs::read_features_input(const char *fnDataIn,
                                       int minFeatures,
                                       gmx_output_env_t *oenv,
                                       LogStream *lstream) {

    std::ifstream fpDataIn;
    bool bPushTime = true;
    std::string line;
    std::vector< std::string > temp;
    real time, tempFeature;
    std::vector<real> tempFeatureVector;   // store a value for making first key-value pair
    std::vector< std::vector<real> > features_local;


    int nFeatures = 0;

    // Open data file, Do it in C++ way for easy reading and parsing
    fpDataIn.open(fnDataIn, std::ifstream::in);
    while(1)	{

        std::getline(fpDataIn, line); // Read each line of the file

        // Stop reading when reached end of the file
        if (fpDataIn.eof())
            break;

        // Split each line with a whitespace
        temp = split(line, ' ');

        // If blank line skip
        if ( temp.empty() )
            continue;

        // Skip line start with # and @, mostly in xvg file
        if ( (temp[0][0] == '#') || (temp[0][0] == '@') )
            continue;


        if (temp[0][0] == '&') {
            bPushTime = false;
            features_local.push_back(tempFeatureVector);
            nFeatures += 1;

            if(minFeatures == nFeatures)
                break;

            // Check if feature array size is same
            if(!bPushTime) {
                if(features_local[0].size() != tempFeatureVector.size()) {
                    gmx_fatal(FARGS,"Size of features array does not match between feature-1 and feature-%d...\n", features_local.size()+1);
                }
            }
            
            // clear and remove memory for tempFeatureVector
            tempFeatureVector.clear();
            tempFeatureVector.shrink_to_fit();

            continue;
        }

        // Push current feature
        tempFeature = (real)std::stod(temp[1]);
        tempFeatureVector.push_back(tempFeature);

        // Store current time
        if(bPushTime) {
            // Also convert time to ps (default unit in trajectory file)
            time = (real)std::stod(temp[0]) * output_env_get_time_invfactor(oenv);
            ClusteringStuffs::timeInInput.push_back(time); 
        }
        
    }

    // change (n_features, time) shape to (time, n_features)
    std::vector< std::vector<real> > features_trans(features_local[0].size(), std::vector<real>(features_local.size()));
    for(size_t i=0;i<features_local.size(); i++) {
        for (size_t j=0;j<features_local[i].size(); j++){
                features_trans[j][i] = features_local[i][j];
        }
    }
    ClusteringStuffs::features = features_trans;

    return TRUE;
}


int ClusteringStuffs::any_central_rmsd_below_thershold(real thres){
    for(size_t i=0; i<this->centrlRmsdMatrix.at(0).size(); i++) {
        for(size_t j=0; j<i; j++) {
            if(centrlRmsdMatrix[i][j] < thres)
                return TRUE;
        }
    }

    return FALSE;
}

int ClusteringStuffs::performClusterMetrics(int eClusterMetrics, int n_clusters, real ssrSstChangeCutoff, LogStream *lstream){
    int finalClustersNumber = 1;
    bool bGotFinalClusterNumber = false;
    real prevSsrSstRatio = ClusteringStuffs::ssrSstRatio.at(1), changeInSsrSstRatio = 0;

    lstream->setprecision(3);
    *lstream<< "\n\n##################### Cluster Metrics Summary #####################\n";
    *lstream<<"Clust. No.\tssr/sst (%)\tDelta(ssr/sst)\tpsuedo F-stat\tDBI\n";
    for(int i = 2; i <= n_clusters; i++){
        changeInSsrSstRatio = ClusteringStuffs::ssrSstRatio.at(i) -prevSsrSstRatio;

        *lstream<<i<<"\t\t";
        *lstream<<ClusteringStuffs::ssrSstRatio.at(i)<<"\t\t";
        *lstream<<changeInSsrSstRatio<<"\t\t";
        *lstream<<ClusteringStuffs::pFS.at(i)<<"\t";
        *lstream<<ClusteringStuffs::dbi.at(i)<<"\n";

        prevSsrSstRatio = ClusteringStuffs::ssrSstRatio.at(i);

        if(finalClustersNumber == 1) {
            finalClustersNumber = i;
            continue;
        }


        if (eClusterMetrics == eSsrSstClusterMetric) {
            if ((changeInSsrSstRatio < ssrSstChangeCutoff) && (!bGotFinalClusterNumber))    {
                finalClustersNumber = i-1;
                bGotFinalClusterNumber = true;
            }
        }

        if (eClusterMetrics == ePfsClusterMetric) {
            //std::cout<<ClusteringStuffs::pFS.at(finalClustersNumber)<<" "<<ClusteringStuffs::pFS.at(i)<<" "<<bGotFinalClusterNumber<<" "<<finalClustersNumber<<"\n";

            if(  ClusteringStuffs::pFS.at(i) > ClusteringStuffs::pFS.at(finalClustersNumber) ) {
                if(!bGotFinalClusterNumber)
                    finalClustersNumber = i;
            }
            else    {
                bGotFinalClusterNumber = true;
            }
        }

        if (eClusterMetrics == eDbiClusterMetric) {
            // std::cout<<ClusteringStuffs::dbi.at(finalClustersNumber)<<" "<<ClusteringStuffs::dbi.at(i)<<" "<<bGotFinalClusterNumber<<" "<<finalClustersNumber<<"\n";

            if(  ClusteringStuffs::dbi.at(i) < ClusteringStuffs::dbi.at(finalClustersNumber)) {
                if(!bGotFinalClusterNumber) {
                    finalClustersNumber = i;
                }
            }
            else {
                bGotFinalClusterNumber = true;
            }

        }
    }
    *lstream<< "##################### ######################### ################### \n";
    lstream->resetprecision();

    return finalClustersNumber;
}


// ############################### ClusteringStuffs END ##########################


void CopyRightMsg() {

    std::string msg = R"~(
               :-)  gmx_clusterByFeatures (-:

             Author: Rajendra Kumar

       Copyright (C) 2018  Rajendra Kumar


gmx_clusterByFeatures is a free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

gmx_clusterByFeatures is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with gmx_clusterByFeatures.  If not, see <http://www.gnu.org/licenses/>.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    )~";

        std::cerr<<msg<<"\n";
}


void write_clustered_trajs(const char *fname, ClusteringStuffs *clustStuff,
                           int *atomIndex, int atomIndexSize,
                           TrajectoryStuffs inpTrajStuff,
                           gmx_bool bAlignTrajToCentral,
                           int *fitAtomIndex, int fitAtomIndexSize,
                           LogStream *lstream)    {

    std::vector< std::string > fnOutTrajs = get_outFile_names(fname, clustStuff->clusterIndex, ".xtc", "");
    t_trxstatus *outStatus;
    int ftp = fn2ftp(inpTrajStuff.filename);
    t_fileio *fio = trx_get_fileio(inpTrajStuff.status);
    real outTime, currentTime;
    int cluster;
    long clusterLength;
    int bRet;

    real *w_rls = nullptr;
    rvec x_shift = {0.0, 0.0, 0.0};
    int i = 0;

    std::cout<<"\n\nWriting trajectory for each cluster...\n";

    //Loop over cluster start here
    for (size_t c=0; c < clustStuff->clusterIndex.size(); c++) {
        outTime = 0;
        cluster = clustStuff->clusterIndex[c];
        clusterLength = clustStuff->clusterDict[cluster].size();

        // Initial fitting stuff
        if( bAlignTrajToCentral) {
            snew(w_rls, inpTrajStuff.atoms.nr);

            if (fitAtomIndexSize < 3)
                gmx_fatal(FARGS, "Need >= 3 points to fit!\n" );

            // Assign weight-factor for fitting
            for (i = 0; i < fitAtomIndexSize; i++) {
                if (inpTrajStuff.atoms.atom[fitAtomIndex[i]].m != 0)
                    w_rls[fitAtomIndex[i]] = inpTrajStuff.atoms.atom[fitAtomIndex[i]].m;
                else
                    w_rls[fitAtomIndex[i]] = 1;
            }

            // Reset to origin and store the translation factor of reference coordinate
            copy_rvec(clustStuff->centralCoords[c][fitAtomIndex[0]], x_shift);
            reset_x(fitAtomIndexSize, fitAtomIndex, inpTrajStuff.atoms.nr, NULL, clustStuff->centralCoords[c], w_rls);
            rvec_dec(x_shift, clustStuff->centralCoords[c][fitAtomIndex[0]]);

        }

        // Open output trajectory here
        outStatus = open_trx(fnOutTrajs[c].c_str(),"w");

        // Loop for frame of each cluster start here
        for(int n = 0; n < clusterLength; n++)  {
            currentTime = ClusteringStuffs::timeInInput[clustStuff->clusterDict[cluster][n]];
            outTime = n*inpTrajStuff.dTime;

             if (ftp == efXTC) {
                bRet = xtc_seek_time(fio, currentTime, inpTrajStuff.natoms, FALSE);
                if (bRet == -1) {
                    gmx_fatal(FARGS, "Frame for this time is not found in trajectory");
                }
                read_next_x(inpTrajStuff.oenv, inpTrajStuff.status, &inpTrajStuff.time, inpTrajStuff.x, inpTrajStuff.box);
            }
            else {
                 if (currentTime < inpTrajStuff.time)
                     rewind_trj(inpTrajStuff.status);


                while (inpTrajStuff.time !=  currentTime)   {
                    bRet = read_next_x(inpTrajStuff.oenv, inpTrajStuff.status, &inpTrajStuff.time, inpTrajStuff.x, inpTrajStuff.box);
                    if (!bRet)
                        break;
                }
            }

            if( bAlignTrajToCentral)    {
                reset_x(fitAtomIndexSize, fitAtomIndex, inpTrajStuff.atoms.nr, NULL, inpTrajStuff.x, w_rls);
                do_fit(inpTrajStuff.natoms, w_rls, clustStuff->centralCoords[c], inpTrajStuff.x);
                for (i = 0; i < inpTrajStuff.natoms; i++)
                    rvec_inc(inpTrajStuff.x[i], x_shift);
            }

            write_trx(outStatus, atomIndexSize, atomIndex, &inpTrajStuff.atoms, 0, outTime, inpTrajStuff.box, inpTrajStuff.x, NULL, NULL);

        } // Loop for frame of each cluster end here

        if( bAlignTrajToCentral)    {
            for (i = 0; i < inpTrajStuff.natoms; i++)
                rvec_inc(clustStuff->centralCoords[c][i], x_shift);

            sfree(w_rls);
        }

        close_trx(outStatus);

    } //Loop over cluster END here
}

std::vector< std::vector< real > > calculate_rmsd(ClusteringStuffs *clustStuff,
                                                  int *fitAtomIndex, int fitAtomIndexSize,
                                                  int *rmsdAtomIndex, int rmsdAtomIndexSize,
                                                  TrajectoryStuffs inpTrajStuff)    {

    real *w_rls, *w_rms, currentTime, temp;
    int i = 0;
    int ftp = fn2ftp(inpTrajStuff.filename);
    t_fileio *fio = trx_get_fileio(inpTrajStuff.status);
    std::vector< std::vector< real > > clusterRMSD;
    rvec x_shift;
    long clusterLength;
    int cluster, bRet;

    std::cout<<"\n\nCalculating RMSD from central structure for each cluster...\n";

    // Weight factor initialization for fitting and RMSD calculation
    snew(w_rms, inpTrajStuff.atoms.nr);
    snew(w_rls, inpTrajStuff.atoms.nr);

    if (fitAtomIndexSize < 3)
        gmx_fatal(FARGS, "Need >= 3 points to fit!\n" );

    // Assign weight-factor for fitting
    for (i = 0; i < fitAtomIndexSize; i++) {
        if (inpTrajStuff.atoms.atom[fitAtomIndex[i]].m != 0)
            w_rls[fitAtomIndex[i]] = inpTrajStuff.atoms.atom[fitAtomIndex[i]].m;
        else
            w_rls[fitAtomIndex[i]] = 1;
    }

    // Assign weight-factor for RMSD calculation
    for (i = 0; i < rmsdAtomIndexSize; i++) {
        if (inpTrajStuff.atoms.atom[rmsdAtomIndex[i]].m != 0)
            w_rms[rmsdAtomIndex[i]] = inpTrajStuff.atoms.atom[rmsdAtomIndex[i]].m;
        else
            w_rms[rmsdAtomIndex[i]] = 1;
    }

    //Loop over cluster start here
    for (size_t c=0; c < clustStuff->clusterIndex.size(); c++) {
        std::vector< real > rmsd;

        cluster = clustStuff->clusterIndex[c];
        clusterLength = clustStuff->clusterDict[cluster].size();


        // Reset to origin and store the translation factor of reference coordinate
        copy_rvec(clustStuff->centralCoords[c][0], x_shift);
        reset_x(fitAtomIndexSize, fitAtomIndex, inpTrajStuff.atoms.nr, NULL, clustStuff->centralCoords[c], w_rls);
        rvec_dec(x_shift, clustStuff->centralCoords[c][0]);

        // Loop for frame of each cluster start here
        for(int n = 0; n < clusterLength; n++)  {
            currentTime = ClusteringStuffs::timeInInput[clustStuff->clusterDict[cluster][n]];


            if (ftp == efXTC) {
                bRet = xtc_seek_time(fio, currentTime, inpTrajStuff.natoms, FALSE);
                if (bRet == -1) {
                    gmx_fatal(FARGS, "Frame for this time is not found in trajectory");
                }
                read_next_x(inpTrajStuff.oenv, inpTrajStuff.status, &inpTrajStuff.time, inpTrajStuff.x, inpTrajStuff.box);
            }
            else {
                // If current time is larger than time in trajectory rewind back
                if (currentTime < inpTrajStuff.time)
                    rewind_trj(inpTrajStuff.status);

                while (inpTrajStuff.time !=  currentTime)   {
                    bRet = read_next_x(inpTrajStuff.oenv, inpTrajStuff.status, &inpTrajStuff.time, inpTrajStuff.x, inpTrajStuff.box);
                    if (!bRet)
                        break;
                }
            }
            // Loop for frame of each cluster end here

            // Fitting to reference structure and RMSD calculation
            reset_x(fitAtomIndexSize, fitAtomIndex, inpTrajStuff.atoms.nr, NULL, inpTrajStuff.x, w_rls);
            do_fit(inpTrajStuff.natoms, w_rls, clustStuff->centralCoords[c], inpTrajStuff.x);
            temp = calc_similar_ind(FALSE, rmsdAtomIndexSize, rmsdAtomIndex, w_rms, inpTrajStuff.x, clustStuff->centralCoords[c]);
            rmsd.push_back( temp );
        }

        clusterRMSD.push_back(rmsd);

        // Translate the central structure to original position
        for (i = 0; i < inpTrajStuff.natoms; i++)
            rvec_inc(clustStuff->centralCoords[c][i], x_shift);

    } //Loop over cluster end here

    sfree(w_rms);
    sfree(w_rls);

    return clusterRMSD;
}


void write_rmsd( std::vector< std::vector< real > > rmsd,
                 const char* fnOutRMSD,
                 std::vector< int > clusterIndex,
                 TrajectoryStuffs inpTrajStuff)    {

    if (fnOutRMSD == NULL) return;   // RETURN HERE IF NO RMSD FILE IS PROVIDED

    std::vector< std::string > fnOutRMSDs = get_outFile_names(fnOutRMSD, clusterIndex, ".xvg", "");
    FILE *fout;
    std::string title;

    for (size_t c=0; c < clusterIndex.size(); c++) {
        title = "RMSD: Cluster-" + std::to_string(clusterIndex[c]);
        fout = xvgropen(fnOutRMSDs[c].c_str(),title.c_str(), output_env_get_time_label(inpTrajStuff.oenv), "RMSD (nm)", inpTrajStuff.oenv);
        for (size_t n=0; n < rmsd[c].size(); n++) {
            //std::cout<<inpTrajStuff.dTime<<" "<<rmsd[c][n]<<std::endl;
            fprintf(fout,"%12.7f   %f \n",inpTrajStuff.dTime*n*output_env_get_time_factor(inpTrajStuff.oenv), rmsd[c][n]);
        }
        xvgrclose(fout);
    }
}


void sort_cluster_frame(std::vector< std::vector< real > > sorter,
                        ClusteringStuffs *clustStuff,
                        std::vector< std::vector< real > > *rmsd )   {

    std::vector< size_t > sortedIndex;
    std::vector< long > clusterFrameIndex;
    std::vector< long > sortedClusterFrameIndex;
    std::vector< real > tempRMSD;
    std::vector< real > sortedRMSD;
    int cluster;

    for(size_t c = 0; c < clustStuff->clusterIndex.size(); c++ )   {
        cluster = clustStuff->clusterIndex[c];
        clusterFrameIndex = clustStuff->clusterDict.at(cluster);
        sortedIndex = argsort(sorter[c]);

        // First reorder in clusterDict
        clustStuff->clusterDict.erase(cluster);      // Erase frame index vector
        for(size_t n=0; n < clusterFrameIndex.size(); n++)  {
            sortedClusterFrameIndex.push_back(clusterFrameIndex[sortedIndex[n]]);
        }
        clustStuff->clusterDict.emplace(cluster, sortedClusterFrameIndex);

        // Clear the sorted one
        sortedClusterFrameIndex.clear();
        sortedClusterFrameIndex.shrink_to_fit();

        // Reorder RMSD
        if ( ! rmsd->empty() ) {
            tempRMSD = rmsd->at(c);
            for(size_t n=0; n < clusterFrameIndex.size(); n++)  {
                sortedRMSD.push_back(tempRMSD[sortedIndex[n]]);
            }
            rmsd->at(c) = sortedRMSD;

            // Clear the sorted onels
            sortedRMSD.clear();
            sortedRMSD.shrink_to_fit();
        }
    }
}



int gmx_clusterByFeatures(int argc,char *argv[])    {

    const char *desc[] = {
        "\"gmx_clusterByFeatures cluster\" can be used to cluster the conformations using the input features.",
        "The features could be any data as a function of time such as Projections of egienvector",
        "from PCA or dihedral-PCA, distances, angles, channel radius etc.[PAR]",
        "[PAR] See more details at https://gmx-clusterbyfeatures.readthedocs.io [PAR]",
        "Clustering methods:",
        " * kmeans: K-means clustring (https://en.wikipedia.org/wiki/K-means_clustering)",
        " * dbscan: Density-based spatial clustering of applications with noise (https://en.wikipedia.org/wiki/DBSCAN)",
        " * gmixture: Gaussina mixture model clustering.",
        " [PAR]",
        "Cluster Metrics (only used with kmeans and gmixture):",
        " * prior: Number of clusters already known.",
        " * rmsd: RMSD between central structures.",
        " * ssr-sst: Relative change in SSR/SST ratio in percentage. Also known as elbow method (https://en.wikipedia.org/wiki/Elbow_method_(clustering))",
        " * pFS: Psuedo F-statatics determined from SSR/SST ratio. Highest value is considered.",
        " * DBI: Davies–Bouldin index (https://en.wikipedia.org/wiki/Davies%E2%80%93Bouldin_index). Lowest value is considered.[PAR]",
        "For summary of command line options. see more details here: https://gmx-clusterbyfeatures.readthedocs.io/en/latest/usage.html[PAR]",
        "For description of command line options, see more details here: https://gmx-clusterbyfeatures.readthedocs.io/en/latest/cmdline.html"
	};

    gmx_bool bAlignTrajToCentral=FALSE, bFit=TRUE;
    int numMinFrameCluster = 20, minFeatures=10;
    const char     *sortMethod[] = { NULL, "none", "rmsd", "features", "user", NULL };
    enum {eNoSort = 1, eSortByRMSD, eSortByFeatures, eSortByUser};
    int eSortMethod;

    int eClusterMetrics;
    real cmRmsdThershold = 0.1;
    real ssrSstChange = 2;

    const char     *clusterAlgo[] = { NULL, "kmeans", "dbscan", "gmixture", NULL };
    enum { eKmeans = 1, eDbscan, eGMixture };
    int eClusterMethod;
    real dbscan_eps = 0.5;
    int dbscan_min_samples = 20;
    int n_clusters=5;

    const char *plotfile = "pca_cluster.png";
    int fontsize = 14;
    real plotHeight = 20;
    real plotWidth = 12;

    t_pargs pa[] =
    {
        { "-method",         FALSE, etENUM, { clusterAlgo },        "Clustering methods. Accepted methods are:" },
        { "-nfeature",       FALSE, etINT,  {&minFeatures},         "Number of features to use for clustering" },
        { "-cmetric",        FALSE, etENUM, { clusterMetrics },     "Cluster metrics: Method to determine cluster number. Accepted methods are:" },
        { "-ncluster",       FALSE, etINT,  {&n_clusters},          "Number of clusters to generate for prior method. Maximum number of cluster for ctrmsd method." },
        { "-crmsthres",      FALSE, etREAL, {&cmRmsdThershold},     "RMSD (nm) threshold between central structures for RMSD cluster metric method." },
        { "-ssrchange",      FALSE, etREAL, {&ssrSstChange},        "Thershold relative change % in SSR/SST ratio for ssr-sst cluster metric method." },
        { "-db_eps",         FALSE, etREAL, {&dbscan_eps},          "The maximum distance between two samples for them to be considered as in the same neighborhood." },
        { "-db_min_samples", FALSE, etINT,  {&dbscan_min_samples},  "The number of samples (or total weight) in a neighborhood for a point to be considered as a core point. This includes the point itself." },
        { "-nminfr",         FALSE, etINT,  {&numMinFrameCluster},  "Number of nimimum frames in a cluster to output it as trajectory" },
        { "-fit",            FALSE, etBOOL, {&bFit},                "Enable fitting and superimposition of the atoms groups different from RMSD/clustering group before RMSD calculation." },
        { "-fit2central",    FALSE, etBOOL, {&bAlignTrajToCentral}, "Enable/Disable trajectory superimposition or fitting to central structure in the output trajectory" },
        { "-sort",           FALSE, etENUM, { sortMethod },         "Sort trajectory according to these values. Accepted methods are" },
        { "-plot",           FALSE, etSTR,  { &plotfile },          "To plot features with clusters in this file." },
        { "-fsize",          FALSE, etINT,  { &fontsize },          "Font size in plot." },
        { "-pltw",           FALSE, etREAL, { &plotWidth },         "Width (inch) of the plot." },
        { "-plth",           FALSE, etREAL, { &plotHeight },        "Height (inch) of the plot."}
    };

	t_filenm   fnm[] = {
            { efTRX, "-f",     NULL,           ffOPTRD },
            { efTPS, NULL,     NULL,           ffOPTRD },
            { efXVG, "-feat",  "feature",      ffOPTRD },
            { efNDX, NULL,     NULL,           ffOPTRD },
            { efXVG, "-clid",  "clid",         ffOPTWR },
            { efLOG, "-g",     "cluster.log",  ffOPTWR },
            { efTRX, "-fout",  "trajout.xtc",  ffOPTWR },
            { efPDB, "-cpdb",  "central.pdb",  ffOPTWR },
            { efXVG, "-rmsd",  "rmsd.xvg",     ffOPTWR }
	};

    #define NFILE asize(fnm)
    gmx_output_env_t *oenv;


    // Copyright message
    CopyRightMsg();

    // Parse command line argument and print all options
    if ( ! parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_TIME_UNIT,NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv) )	{
        return 0;
    }

    // Sorted and clustering algorithm method
    eSortMethod = nenum(sortMethod);
    eClusterMethod = nenum(clusterAlgo);
    eClusterMetrics = nenum(clusterMetrics);

    const char *fnClId, *fnFeatures;

    // Variables realting to clustering stuffs
    ClusteringStuffs *clustStuff = nullptr, *tempClustStuff;
    std::map < int, ClusteringStuffs* >  allClusterStuffs;
    std::vector< std::vector< real > > clusterRMSD; // RMSD of clusters with reference to central structure
    int finalClustersNumber = 1;
    PyCluster pycluster;

    //OUTPUT FILE STUFFS
    const char *fnOutPDB = NULL, *fnOutRMSD = NULL, *fnOutLog=NULL;
    std::vector< std::string > fnOutPDBs;

    //INDEX RELATED VARIABLES
    int outIndexSize, *outIndex; // Output atom index group
    int *fitAtomIndex, fitAtomIndexSize; // Fitting atom index group
    int *rmsdAtomIndex, rmsdAtomIndexSize; // Group for which clustering has been performed, also used for RMSD calculation
    char *grpnm;

    gmx_bool bFeatures = FALSE, bCentralPDB = FALSE, bDoCluster = FALSE, bTrajRMSD=FALSE;

    t_topology top;
    int ePBC;
    rvec *x;
    matrix box;

    //TRAJECTORY RELATED VARIABLES
    TrajectoryStuffs inpTrajStuff;
    const char *inpTrajName, *fnOutTraj;


    fnFeatures = opt2fn_null("-feat",NFILE,fnm);
    fnClId = opt2fn_null("-clid",NFILE,fnm);

    if((fnFeatures == NULL) && (fnClId == NULL)) {
        gmx_fatal(FARGS,"Input files with either features or cluster-ids and features are missing!!!\n");
    }

    if(fnFeatures != NULL){
        bDoCluster = TRUE;
        bFeatures = TRUE;
        fnClId = opt2fn("-clid",NFILE,fnm);
    }
    else {
        // Here it only extract the file name as file is input
        fnClId = opt2fn("-clid",NFILE,fnm);
    }

    // check whether sort by features is enabled and store it for later use
    ClusteringStuffs::bSortByFeatures = (eSortMethod == eSortByFeatures) ? true : false;

    //Initialize log output
    fnOutLog = opt2fn("-g", NFILE, fnm);
    LogStream lstream(fnOutLog);
    lstream<<"=======================\n";
    lstream<<"  Cluster Log output   \n";
    lstream<<"=======================\n";
    lstream<<"\nCommand:\n=======================\n";
    lstream<<gmx::CommandLineProgramContext(argc, argv).commandLine()<<"\n";
    lstream<<"=======================\n";

    // Read input cluster-id file, if it is the input
    if (!bDoCluster)    {
        clustStuff = new ClusteringStuffs();
        if (!(clustStuff->read_cluster_input(fnClId, &bFeatures, numMinFrameCluster, oenv, &lstream )))
            return EXIT_FAILURE;

            /* May be enable it in future
        // If given separately, read features file and store in ClusteringStuffs class static variable
        if( (!bFeatures) && (fnFeatures != NULL) )  {
            if (!(ClusteringStuffs::read_features_input(fnFeatures, minFeatures, oenv, &lstream)))    {
                gmx_fatal(FARGS,"Not able to read features file!!!\n");
              }
            bFeatures = TRUE;
          }*/
      }

    // Check if trajectory is given as input and save filename for later use
    inpTrajName = opt2fn_null("-f",NFILE,fnm);
    if(inpTrajName != NULL) {
      inpTrajStuff.filename = inpTrajName;
      inpTrajStuff.bTraj = true;
     }


    // Check if alignment to central structure is possible
    if ( (bAlignTrajToCentral) && (!bFeatures) ){
        lstream<<"\n\n======== WARNING ========\n";
        lstream<<"\nFeatures file is missing. Central structure cannot be calculated !!!\n";
        lstream<<"Switching-off -fit2central option.\n";
        bAlignTrajToCentral = FALSE;
    }


    // If features are in file, only calculate central structure otherwise exit here
    if((opt2fn_null("-cpdb", NFILE, fnm) != NULL) && (bFeatures) )
        bCentralPDB = TRUE;
    if ( !bFeatures && bCentralPDB ){
        lstream<<"\n\n======== WARNING ========\n";
        lstream<<"\nFeatures file is missing. Not able to write cnetral structures. \n";
        bCentralPDB = FALSE;
    }

    // If trajectory is not given, cnetral structures cannot be extracted. Also, alignment to trajectory cannot be done without trajectory
    if ( (bCentralPDB) && (!inpTrajStuff.bTraj) ){
        lstream<<"\n\n======== WARNING ========\n";
        lstream<<"\nTrajectory file is miising. Central structure needs to be extracted from trajectory!!!\n";
        lstream<<"Central structure will not be wriiten as pdb file. \n";
        bCentralPDB = FALSE;
        bAlignTrajToCentral = FALSE;
    }

    // If sorted method is given, check for RMSD calculation
    if ( (eSortMethod == eSortByRMSD) || (opt2fn_null("-rmsd", NFILE, fnm) != NULL)) {
        if (inpTrajStuff.bTraj)
            bTrajRMSD = TRUE;
        else    {
            lstream<<"\n\n======== WARNING ========\n";
            lstream<<"\nTrajectory file is miising. Central structure needs to be extracted from trajectory!!!\n";
            lstream<<"Therefore, RMSD w.r.t central structure cannot be calculated and trajectory cannot be sorted.\n";
            eSortMethod = eNoSort;
        }
    }

    // If trajectory is not given and RMSD cluster metric is given as input, change metric to ssr-sst
    if( (eClusterMetrics == eCRmsdClusterMetric) && (!inpTrajStuff.bTraj) ) {
        lstream<<"\n\n======== WARNING ========\n";
        lstream<<"\nTrajectory file is missing. Central structure needs to be extracted from trajectory!!!\n";
        lstream<<"Therefore, RMSD between central structure cannot be calculated.\n";
        lstream<<"Switching to ssr-sst metric with default value.\n";
        eClusterMetrics = eSsrSstClusterMetric;
    }

    // If trajectory is not given and output clustered trajectory is requested, exit here.
    if( (opt2fn_null("-fout", NFILE, fnm) != NULL) && (!inpTrajStuff.bTraj) ) {
        lstream<<"\n\n======== ERROR ========\n";
        lstream<<"\nInput trajectory file is missing while clustered trajectory is requested for output.\n";
        lstream<<"Exiting...\n";
        exit(1);
    }

    // Reading tpr, index and trajectory file
    if(inpTrajStuff.bTraj) {
        read_tps_conf(ftp2fn(efTPS,NFILE,fnm), &top, &ePBC, &x, NULL, box, FALSE);
        inpTrajStuff.atoms = top.atoms;
        inpTrajStuff.ePBC = ePBC;


        // Selection of output index group
        printf("\nChoose a group for the output:\n");
        get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&outIndexSize,&outIndex,&grpnm);

        if(fnFeatures != NULL) {
          printf("\nChoose a group for clustering/RMSD calculation:\n");
          get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&rmsdAtomIndexSize,&rmsdAtomIndex,&grpnm);

          // Selection for fitting group
          if (bFit)    {
              printf("\nChoose a group for fitting or superposition:\n");
              get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&fitAtomIndexSize,&fitAtomIndex,&grpnm);
            }
          else {
            fitAtomIndexSize = rmsdAtomIndexSize;
            fitAtomIndex = rmsdAtomIndex;
            }
          }

      // Read first frame of the input trajectory
      inpTrajStuff.natoms = read_first_x(oenv, &inpTrajStuff.status, inpTrajName, &inpTrajStuff.time, &inpTrajStuff.x, inpTrajStuff.box);
      inpTrajStuff.oenv = oenv;
      set_dTime(&inpTrajStuff);

      }


    if(bFeatures)  {
        // Initialize python and clustering code
        PyCluster pycluster = PyCluster();
        pycluster.InitPythonAndLoadFunc();
        pycluster.initializeClustering(fnFeatures, minFeatures, clusterAlgo[eClusterMethod], (float)dbscan_eps, dbscan_min_samples);
    }

    if(bDoCluster) {

        // Read features file and store in ClusteringStuffs class static variable
        if (!(ClusteringStuffs::read_features_input(fnFeatures, minFeatures, oenv, &lstream)))    {
            gmx_fatal(FARGS,"Not able to read features file!!!\n");
        }

        // Loop to find number of clusters, In case of prior method, run loop only once
        int curr_n_cluster;
        if ( (eClusterMetrics == eCRmsdClusterMetric) || (eClusterMetrics == ePriorClusterMetric)) {
            // Start with maximum number of cluster
            curr_n_cluster = n_clusters;
        }
        else {
            // Start with minimum number of clusters
            curr_n_cluster = 1;
        }

        while(1)    {
            double tempSsrSstRatio, tempPFS;
            lstream<<"\n###########################################\n";
            lstream<<"########## NUMBER OF CLUSTERS : "<<curr_n_cluster<<" ########\n";
            lstream<<"###########################################\n";

            // Initialize ClusteringStuffs
            clustStuff = new ClusteringStuffs();
            allClusterStuffs.emplace(curr_n_cluster, clustStuff);


            // Perform clustering
            pycluster.performClustering(curr_n_cluster);
            clustStuff->clidAlongTime = pycluster.getClusterLabels(curr_n_cluster);

            // Construct cluster dictionary and cluster-index
            clustStuff->constructClusterDict(numMinFrameCluster, &lstream);

            // Determine index of central structure
            if(!(clustStuff->calculate_central_struct(&lstream)))
                return EXIT_FAILURE;

            // Extract central structures and calculate RMSD between them
            if(inpTrajStuff.bTraj)  {
                clustStuff->write_central_pdbfiles(fnOutPDBs, outIndex, outIndexSize, inpTrajStuff);
                if (curr_n_cluster > 1) {
                    clustStuff->rmsd_bw_central_structure(fitAtomIndex, fitAtomIndexSize, \
                                                 rmsdAtomIndex, rmsdAtomIndexSize,\
                                                 inpTrajStuff, &lstream);
                  }
              }


            if ( (eClusterMetrics == ePriorClusterMetric) || (eClusterMethod == eDbscan) ){
                // If cluster-metric is not needed or DBSCAN method is used, iterate only once

                finalClustersNumber = curr_n_cluster;
                break;
            }
            else if (eClusterMetrics == eCRmsdClusterMetric)  {
                // If RMSD is used, start with maximum clusters number and reduce it untill criteria met

                // If cluster-count is one -- break here
                if (curr_n_cluster == 1) {
                    finalClustersNumber = curr_n_cluster;
                    break;
                }

                // use RMSD threshold here to break
                if (!clustStuff->any_central_rmsd_below_thershold(cmRmsdThershold)) {
                    finalClustersNumber = curr_n_cluster;
                    break;
                }

                // reduce cluster count for next iteration
                curr_n_cluster = curr_n_cluster - 1;
            }
            else {
                // If any of others cluster-metrics, ssr-sst ratio, Psuedo F-statistics and DB index is used,
                // Start with clusters number one and increase it to maximum clusters number.
                // For each increased number, compute all quantities
                pycluster.getSsrSstStats(curr_n_cluster, &tempSsrSstRatio,  &tempPFS);
                ClusteringStuffs::ssrSstRatio.emplace(curr_n_cluster, (real)tempSsrSstRatio);
                ClusteringStuffs::pFS.emplace(curr_n_cluster, (real)tempPFS);
                clustStuff->calculateDaviesBouldinIndex();

                if (curr_n_cluster == n_clusters) {
                    break;
                }
                curr_n_cluster = curr_n_cluster + 1;
            }
        }


        // If cluster-metrics, ssr-sst ratio, Psuedo F-statistics and DB index is used, here print summary
        // and select final clusters-number according to the input criteria.
        if ( ! ((eClusterMetrics == eCRmsdClusterMetric) || (eClusterMetrics == ePriorClusterMetric)) ) {
            finalClustersNumber = ClusteringStuffs::performClusterMetrics(eClusterMetrics, n_clusters, ssrSstChange, &lstream);
        }

        lstream<<"\n\n#####################################\n";
        lstream<<"Final number of cluster selected: "<<finalClustersNumber;
        lstream<<"\n#####################################\n";

        // Load the final one here
        clustStuff = allClusterStuffs.at(finalClustersNumber);

        // Remove all other ClustStuffs
        for (std::map< int, ClusteringStuffs* >::iterator it=allClusterStuffs.begin(); it!=allClusterStuffs.end(); ++it)  {
            // std::cout<<"Hellooo "<<it->first<<" "<<it->second<<"\n";
            if(it->first == finalClustersNumber)
                continue;
             tempClustStuff = allClusterStuffs.at(it->first);
             allClusterStuffs.emplace(it->first, nullptr);
             delete tempClustStuff;
        }
        allClusterStuffs.clear();

    }
    else {
        if (bFeatures)  {
            // Determine index of central structure
            if(!(clustStuff->calculate_central_struct(&lstream)))
                return EXIT_FAILURE;

            // Extract central structures and calculate RMSD between them
            if(inpTrajStuff.bTraj)  {
                clustStuff->write_central_pdbfiles(fnOutPDBs, outIndex, outIndexSize, inpTrajStuff);
                clustStuff->rmsd_bw_central_structure(fitAtomIndex, fitAtomIndexSize, \
                                                      rmsdAtomIndex, rmsdAtomIndexSize,\
                                                      inpTrajStuff, &lstream);
              }
        }
    }

    // Plot the features by cluster
    if((bFeatures) && (opt2parg_bSet("-plot", asize(pa), pa)) ){
        pycluster.plotFeaturesClusters(finalClustersNumber, plotfile, clustStuff->get_central_ids(), fontsize, plotHeight, plotWidth);
    }

    // Write cluster-id file
    if(bDoCluster) {
        FILE *fClId;
        fClId = xvgropen(fnClId, "Cluster-ID", "Time", "Cluster-ID", oenv);
        for(unsigned long i=0; i<clustStuff->timeInInput.size(); i++){
            fprintf(fClId, "%5.3f    %d\n", clustStuff->timeInInput.at(i), clustStuff->clidAlongTime.at(i));
        }
        xvgrclose(fClId);
    }

    // Write final central pdb files
    if (bCentralPDB) {
        fnOutPDB = opt2fn_null("-cpdb", NFILE, fnm);

        fnOutPDBs = get_outFile_names(fnOutPDB, clustStuff->clusterIndex, ".pdb", "");

        clustStuff->write_central_pdbfiles(fnOutPDBs, outIndex, outIndexSize, inpTrajStuff);

    }

    //std::cout<<sortedClusterIds.size()<<" "<<this->centralStructDict.size()<<"\n";


    // Calculate RMSD of trajectories
    if ((bTrajRMSD) || (eSortMethod == eSortByRMSD))    {
        clusterRMSD = calculate_rmsd(clustStuff, fitAtomIndex, fitAtomIndexSize, rmsdAtomIndex, rmsdAtomIndexSize, \
                                     inpTrajStuff);
    }



    if (eSortMethod == eSortByRMSD) {
        sort_cluster_frame(clusterRMSD, clustStuff, &clusterRMSD );
    }


    if ( (eSortMethod == eSortByFeatures) && (clustStuff->avgDistanceToAllDict.empty() == false) ) {
        std::vector< std::vector< real > > sorter;
        for(size_t i=0; i<clustStuff->clusterIndex.size(); i++) {
            sorter.push_back( clustStuff->avgDistanceToAllDict.at(clustStuff->clusterIndex.at(i)) );
        }
        sort_cluster_frame(sorter, clustStuff, &clusterRMSD );
    }

    // Write clustered trajectories, but only when input trajectory is present
    if(inpTrajStuff.bTraj)  {
        fnOutTraj = opt2fn_null("-fout", NFILE, fnm);
        if(fnOutTraj != NULL)
          write_clustered_trajs(fnOutTraj, clustStuff, outIndex, outIndexSize, \
                                inpTrajStuff, bAlignTrajToCentral, fitAtomIndex, fitAtomIndexSize, &lstream);
      }

    // Calculate and write RMSD
    if (bTrajRMSD)  {
        fnOutRMSD = opt2fn_null("-rmsd", NFILE, fnm);
        if (fnOutRMSD != NULL)
            write_rmsd( clusterRMSD, fnOutRMSD, clustStuff->clusterIndex, inpTrajStuff);
    }

    return 0;
}
