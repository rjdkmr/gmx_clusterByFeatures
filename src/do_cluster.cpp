/*
 * This file is part of gmx_clusterByFeatures
 *
 * Author: Rajendra Kumar
 * Copyright (C) 2018  Rajendra Kumar
 *
 * g_coordNdata is a free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * g_coordNdata is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with g_coordNdata.  If not, see <http://www.gnu.org/licenses/>.
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

    std::string code = R"~(import numpy as np
from sklearn import cluster as getCluster
from sklearn import mixture
import re, os, sys
import shlex, subprocess, shutil

class DoClustering:
    algo = 'kmeans'
    dbscan_eps = 0.5
    dbscan_min_samples=20

    features = None
    time = []
    labels = dict()
    sse = dict()

    #########################################################################################
    def __init__(self, filename, nFeatures=2, algo='kmeans', dbscan_eps=0.5, dbscan_min_samples=20):
        self.algo = algo
        self.dbscan_eps = dbscan_eps
        self.dbscan_min_samples = dbscan_min_samples

        # Read features file here
        fin = open(filename, 'r')
        nframes = None
        time = []
        coords = []
        features = []

        ppc = []
        ecount = 0
        for line in fin:
            line = line.lstrip().rstrip()
            if not line.strip():
                continue

            if re.search('&', line) is not None:
                ecount += 1
                features.append(np.asarray(ppc))
                ppc = []
                if ecount == nFeatures:
                    break
                continue

            if re.search('#|@', line) is None:
                temp = re.split('\s+', line)
                if ecount == 0:
                    self.time.append(float(temp[0]))
                ppc.append(float(temp[1]))

        fin.close()
        self.features = np.asarray(features).T

    #########################################################################################
    def calculate_clusters(self, n_clusters):
        if self.algo == 'kmeans':
            db = getCluster.KMeans(n_clusters=n_clusters)

        if self.algo == 'dbscan':
            db = getCluster.DBSCAN(eps=self.dbscan_eps, min_samples=self.dbscan_min_samples)

        if self.algo == 'gaussmix':
            db = mixture.GaussianMixture(n_components=n_clusters, covariance_type='full')

        db.fit(self.features)

        if hasattr(db, 'labels_'):
            labels = db.labels_.astype(np.int)
        else:
            labels = db.predict(features)

        trueIdx = np.nonzero(labels >= 0)
        labels[trueIdx] = labels[trueIdx] + 1

        self.labels[n_clusters] = list(self._sort_clusters(labels, n_clusters))
        self.sse[n_clusters] = db.inertia_

    #########################################################################################
    def _sort_clusters(self, labels, n_clusters):
        clusters = dict()
        clid = []

        # Store index of each cluster
        for i in range(len(labels)):
            if labels[i] == -1:
                continue
            if labels[i] not in clusters:
                clusters[labels[i]] = [ i ]
                clid.append(labels[i])
            else:
                clusters[labels[i]] = clusters[labels[i]] + [ i ]

        # Make a sorted list of clusters-id
        clid = list(sorted(clid))
        length = []
        for c in clid:
            length.append( len(clusters[c]) )

        # Change the cluster-ids using stored index above
        idx = np.argsort(length)[::-1]
        newIdx = 1
        for i in idx:
            labels[ clusters[ clid[i] ] ] = newIdx
            newIdx += 1

        return labels

    #########################################################################################
    def plotFeaturesClusters(self, n_clusters, plotfile, central_id=None, fsize=14, width=12, height=20):
        import matplotlib as mpl

        for gui in mpl.rcsetup.non_interactive_bk:
            try:
                mpl.use(gui, warn=True, force=True)
                from matplotlib import pyplot as plt
                break
            except:
                continue

        labels = self.labels[n_clusters]

        fig = plt.figure(figsize=(height, width))
        fig.subplots_adjust(top=0.95, bottom=0.1, wspace=0.3, hspace=0.5)
        mpl.rcParams['font.size'] = fsize
        handles, legend_labels = None, None
        xNewT = self.features.T

        length = len(xNewT)
        if length > 6:
            length = 6

        axCounter = 1
        for pc1 in range(length):
            for pc2 in range(pc1):
                ax = fig.add_subplot(8,2,axCounter)
                axCounter += 1
                for l in set(labels):
                    if l == -1:
                        ax.scatter(xNewT[pc1][labels == l], xNewT[pc2][labels == l], s=0.5, c='k')
                    else:
                        ax.scatter(xNewT[pc1][labels == l], xNewT[pc2][labels == l], s=0.5, label=str(l))

                    handles, legend_labels = ax.get_legend_handles_labels()

                if central_id is not None:
                    for t in central_id:
                        ax.scatter(xNewT[pc1][t], xNewT[pc2][t], s=12, c='k')

                ax.set_xlabel('feature-{0}'.format(pc1+1))
                ax.set_ylabel('feature-{0}'.format(pc2+1))

        fig.legend(handles, legend_labels, ncol=8, loc='upper center',scatterpoints=5,markerscale=6)
        plt.savefig(plotfile, dpi=300)

    #########################################################################################
    def get_labels(self, n_clusters):
        return self.labels[n_clusters]

    #########################################################################################
    def get_ssr_sst_stats(self, n_clusters):
        sst = self.sse[1]
        ssr = sst - self.sse[n_clusters]
        ratio = ssr/sst * 100
        if n_clusters != 1:
            pFS = (ssr/(n_clusters-1)) / (self.sse[n_clusters] /(self.features.shape[0]-n_clusters) )
        else:
            pFS = 0.0

        return (ratio, pFS)


)~";

    return code;
}



void PyCluster::initializeClustering(const char* filename, int nFeatures, const char* algo, float dbscan_eps, int dbscan_min_samples) {
    std::stringstream code;
    code<<"doCluster = DoClustering( ";
    code<<"\'"<<filename<<"\', ";
    code<<"nFeatures= "<<nFeatures<<", ";
    code<<"algo=\'"<<algo<<"\', ";
    code<<"dbscan_eps="<<dbscan_eps<<", ";
    code<<"dbscan_min_samples="<<dbscan_min_samples;
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

    if (pyList==NULL)   {
        std::cout<<"ERROR: Error in python execution. No python object returned from getClusterLabels(). \n";
        exit(1);
    }

    // Read Python list and convert back it to C++ vector
    for (size_t i =0; i<pyList.size(); i++){
        labels.push_back( pyList[i].cast<int>() );
    }

    return labels;
}

void PyCluster::getSsrSstStats(int n_clusters, double *ratio, double *pFS) {
    py::tuple pyReturnValues;

    // Run the code and return python tuple
    pyReturnValues = py::eval("doCluster.get_ssr_sst_stats( " + std::to_string(n_clusters) + " ) \n", PyCluster::scope );

    if (pyReturnValues==NULL)   {
        std::cout<<"ERROR: Error in python execution. No python object returned from getSsrSstStats().\n";
        exit(1);
    }

    *ratio = pyReturnValues[0].cast<float>() ;
    *pFS = pyReturnValues[1].cast<float>() ;

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
