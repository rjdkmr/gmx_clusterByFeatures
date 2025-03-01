/*
 * This file is part of gmx_clusterByFeatures
 *
 * Author: Rajendra Kumar
 * Copyright (C) 2018-2025  Rajendra Kumar
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


#ifndef DO_CLUSTER_H
#define DO_CLUSTER_H

#include <pybind11/eval.h>

namespace py = pybind11;

class PyCluster {

public:
    // Evaluate in scope of main module
    static py::object scope;
    static void InitPythonAndLoadFunc(); //Initialize Python
    static std::string PyDoClusteringClusteringClassCode(); // Read DoClustering class defintion

    /*
     * Initialize DoClustering Class as doCluster object in python
     */
    static void initializeClustering(const char* filename, int nPC=2, const char* algo="kmeans", float dbscan_eps=0.5, int dbscan_min_samples=20, float silhouette_score_sample_size=20);

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
    static void getClusterMetrics(int n_clusters, double *ratio,  double *pFS, double *silhouette_score, double *davies_bouldin_score);

    /*
     * Plot the features with clusters.
     */
    static void plotFeaturesClusters( int n_clusters, const char *plotfile, std::vector< long > central_id, int fsize=14, float width=12, float height=20);
};

#endif // DO_CLUSTER_H

