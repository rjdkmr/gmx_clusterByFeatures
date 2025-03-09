import numpy as np
from sklearn import cluster as getCluster
from sklearn import mixture
from sklearn import metrics
import re, os, sys
import shlex, subprocess, shutil

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

        db.fit(self.features)

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

        self.labels[n_clusters] = list(self._sort_clusters(labels, n_clusters))
        if hasattr(db, 'inertia_'):
            self.sse[n_clusters] = db.inertia_
        else:
            self.sse[n_clusters] = 1

    #########################################################################################
    def _sort_clusters(self, labels, n_clusters):
        if n_clusters == 1:
            return labels

        clusterIds = sorted(list(set(list(labels))))
        length = []
        for cid in clusterIds:
            length.append(np.sum(labels == cid))

        # Change the cluster-ids using stored index above
        sorted_by_length_idx = np.argsort(length)[::-1]
        newIdx = 1
        newLabels = np.ones(labels.shape, dtype=int)
        for old_cid_idx in sorted_by_length_idx:
            newLabels[ np.nonzero(labels == old_cid_idx+1) ] = newIdx
            newIdx += 1

        return newLabels

    #########################################################################################
    def plotFeaturesClusters(self, n_clusters, plotfile, central_id=None, fsize=14, width=12, height=20):
        import matplotlib as mpl

        for gui in mpl.rcsetup.non_interactive_bk:
            try:
                mpl.use(gui, force=True)
                from matplotlib import pyplot as plt
                break
            except:
                continue

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

