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

#include <pybind11/eval.h>

#include <string>
#include <vector>
#include <iostream>
#include <sstream>


#include "do_cluster.h"

namespace py = pybind11;

py::object PyCluster::scope;

void PyCluster::InitPythonAndLoadFunc(){
    py::object scope = py::module::import("__main__").attr("__dict__");
    py::exec(PyCluster::PyDoClusteringClusteringClassCode(), scope);
    PyCluster::scope = scope;
}

std::string PyCluster::PyDoClusteringClusteringClassCode() {
    
    // Read hex sequences of code as a string 
    std::string code = { 
        #include "cluster.pyhex"
    };

    return code;
}



void PyCluster::initializeClustering(const char* filename, int nFeatures, const char* algo, float dbscan_eps, int dbscan_min_samples, float silhouette_score_sample_size) {
    std::stringstream code;
    code<<"doCluster = DoClustering( ";
    code<<"\'"<<filename<<"\', ";
    code<<"nFeatures= "<<nFeatures<<", ";
    code<<"algo=\'"<<algo<<"\', ";
    code<<"dbscan_eps="<<dbscan_eps<<", ";
    code<<"dbscan_min_samples="<<dbscan_min_samples<<", ";
    code<<"silhouette_score_sample_size="<<silhouette_score_sample_size;
    code<<")";
    py::exec(code.str(), PyCluster::scope);
}

void PyCluster::performClustering(int n_clusters){
    std::stringstream code;
    code<<"doCluster.calculate_clusters("<<n_clusters<<")";
    py::exec(code.str(), PyCluster::scope);
}

std::vector< int > PyCluster::getClusterLabels(int n_clusters){
    std::vector< int > labels;
    py::list pyList;

    // Run the code and return python list
    pyList = py::eval("doCluster.get_labels( " + std::to_string(n_clusters) + " ) \n", PyCluster::scope);

    if (pyList == NULL)   {
        std::cout<<"ERROR: Error in python execution. No python object returned from getClusterLabels(). \n";
        exit(1);
    }

    // Read Python list and convert back it to C++ vector
    for (size_t i =0; i<pyList.size(); i++){
        labels.push_back( pyList[i].cast<int>() );
    }

    return labels;
}

void PyCluster::getClusterMetrics(int n_clusters, double *ratio, double *pFS, double *silhouette_score, double *davies_bouldin_score) {
    py::tuple pyReturnValues;

    // Run the code and return python tuple
    pyReturnValues = py::eval("doCluster.get_cluster_metrics( " + std::to_string(n_clusters) + " ) \n", PyCluster::scope );

    if (pyReturnValues == NULL)   {
        std::cout<<"ERROR: Error in python execution. No python object returned from getSsrSstStats().\n";
        exit(1);
    }

    *ratio = pyReturnValues[0].cast<float>() ;
    *pFS = pyReturnValues[1].cast<float>() ;
    *silhouette_score = pyReturnValues[2].cast<float>() ;
    *davies_bouldin_score = pyReturnValues[3].cast<float>() ;

    //std::cout<<*ratio<<" "<<*pFS<<"\n";
}

void PyCluster::plotFeaturesClusters( int n_clusters, const char *plotfile,
                                     std::vector< long > central_id,
                                     int fsize, float width, float height) {
    std::stringstream code;

    // First assign central_id variable
    code<<"central_id = [ ";
    for (size_t i=0; i<central_id.size(); i++)
        code<<central_id[i]<<", ";
    code<<" ]";
    py::exec(code.str(), PyCluster::scope);

    code.str(""); // clear the stringstream
    code.clear();

    code<<"doCluster.plotFeaturesClusters("<<n_clusters<<", ";
    code<<" \'"<<plotfile<<"\', ";
    code<<"central_id=central_id, ";
    code<<"fsize="<<fsize<<", ";
    code<<"width="<<width<<", ";
    code<<"height="<<height;
    code<<")";
    py::exec(code.str(), PyCluster::scope);
}
