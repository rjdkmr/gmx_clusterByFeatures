import numpy as np
from sklearn import cluster as getCluster
from sklearn import mixture
from sklearn import metrics
import re, os, sys
import shlex, subprocess, shutil
import matplotlib as mpl

mpl.backends.backend_registry.list_builtin(mpl.backends.BackendFilter.NON_INTERACTIVE)
from matplotlib import pyplot as plt

try:
    import hdbscan
    has_hdbscan = True
except ImportError:
    has_hdbscan = False

class DoClustering:
    algo = 'kmeans'
    dbscan_eps = 0.5
    dbscan_min_samples=20
    silhouette_score_sample_size = 50

    features = None
    nframes = None
    time = []
    labels = dict()
    sse = dict()
    silhouette_score = dict()
    davies_bouldin_score = dict()

    #########################################################################################
    def __init__(self, filename, nFeatures=2, algo='kmeans', dbscan_eps=0.5, dbscan_min_samples=20, silhouette_score_sample_size=20):
        
        if algo == 'hdbscan' and not has_hdbscan:
            print("HDBSCAN library not found. Please install it or use another clustering algorithm.")
            sys.exit(1)

        self.algo = algo
        self.dbscan_eps = dbscan_eps
        self.dbscan_min_samples = dbscan_min_samples
        self.silhouette_score_sample_size = None

        # Read features file here
        fin = open(filename, 'r')
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
        self.nframes = len(self.time)
        self.features = np.asarray(features).T

        if self.nframes > 10000:
            self.silhouette_score_sample_size = 10000
        else:
            self.silhouette_score_sample_size = silhouette_score_sample_size

    #########################################################################################
    def calculate_clusters(self, n_clusters):
        if self.algo == 'kmeans' and self.nframes <= 100000:
            db = getCluster.KMeans(n_clusters=n_clusters, n_init=5, random_state=np.random.RandomState(12345))

        if self.algo == 'kmeans' and self.nframes > 100000:
            db = getCluster.MiniBatchKMeans(n_clusters=n_clusters, n_init=5, random_state=np.random.RandomState(12345))

        if self.algo == 'dbscan':
            db = getCluster.DBSCAN(eps=self.dbscan_eps, min_samples=self.dbscan_min_samples)

        if self.algo == 'gmixture':
            db = mixture.GaussianMixture(n_components=n_clusters, covariance_type='full')

        if self.algo == 'hdbscan':
            db = hdbscan.HDBSCAN(cluster_selection_epsilon=self.dbscan_eps, min_samples=self.dbscan_min_samples )

        db.fit(self.features)

        if self.algo in ['dbscan', 'hdbscan']:
            n_clusters = len(set(db.labels_)) - (1 if -1 in db.labels_ else 0)

        if hasattr(db, 'labels_'):
            labels = db.labels_.astype(int)
        else:
            labels = db.predict(self.features)
            
        if n_clusters > 1:
            self.silhouette_score[n_clusters] = metrics.silhouette_score(self.features, labels, sample_size=self.silhouette_score_sample_size)
            self.davies_bouldin_score[n_clusters] = metrics.davies_bouldin_score(self.features, labels)
        else:
            self.silhouette_score[n_clusters] = 0
            self.davies_bouldin_score[n_clusters] = 0

        trueIdx = np.nonzero(labels >= 0)
        labels[trueIdx] = labels[trueIdx] + 1

        ordered_newLabels, old_to_new_ordered_map = self._sort_clusters(labels, n_clusters)
        self.labels[n_clusters] = list(ordered_newLabels)
        if hasattr(db, 'inertia_'):
            self.sse[n_clusters] = db.inertia_
        else:
            self.sse[n_clusters] = 1

        if self.algo == 'hdbscan':
            cmap_list = [(0, '#c2c0c1'), (0.25, '#46a6e4'), (0.75, '#c01755'), (1.0, '#000000')]
            cmap = mpl.colors.LinearSegmentedColormap.from_list('dummy', cmap_list, N=n_clusters)
            colors = cmap(np.linspace(0, 1, n_clusters))
            ordered_colors = [colors[old_to_new_ordered_map[i]-1] for i in range(1, n_clusters+1)]
            fig = plt.figure(figsize=(11, 8))
            ax = fig.add_subplot(1,1,1)
            db.condensed_tree_.plot(select_clusters=True, axis=ax, selection_palette=ordered_colors)
            fig.savefig('hdbscan_condensed_tree.png', dpi=300)

        return n_clusters

    #########################################################################################
    def _sort_clusters(self, labels, n_clusters):
        if n_clusters == 1:
            return labels, {1: 1}

        clusterIds = sorted(list(set(list(labels))))
        length = []
        for cid in clusterIds:
            if cid != -1:
                length.append(np.sum(labels == cid))

        # Change the cluster-ids using stored index above
        sorted_by_length_idx = np.argsort(length)[::-1]
        newIdx = 1
        newLabels = np.ones(labels.shape, dtype=int) * -1
        old_to_new_map = dict()
        for old_cid_idx in sorted_by_length_idx:
            newLabels[ np.nonzero(labels == old_cid_idx+1) ] = newIdx
            old_to_new_map[old_cid_idx+1] = newIdx
            newIdx += 1

        return newLabels, old_to_new_map

    #########################################################################################
    def plotFeaturesClusters(self, n_clusters, plotfile, central_id=None, fsize=14, width=12, height=20):
        labels = self.labels[n_clusters]

        fig = plt.figure(figsize=(width, height))
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
    def get_cluster_metrics(self, n_clusters):
        sst = self.sse[1]
        ssr = sst - self.sse[n_clusters]
        ratio = ssr/sst * 100
        
        if n_clusters != 1:
            pFS = (ssr/(n_clusters-1)) / (self.sse[n_clusters] /(self.features.shape[0]-n_clusters) )
        else:
            pFS = 0.0
            
        return (ratio, pFS, self.silhouette_score[n_clusters], self.davies_bouldin_score[n_clusters])

