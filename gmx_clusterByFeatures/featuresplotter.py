"""
This file is part of gmx_clusterByFeatures

Author: Rajendra Kumar

Copyright (C) 2014-2025  Rajendra Kumar

gmx_clusterByFeatures is a free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

"""

import csv
import numpy as np
import sys
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from collections import OrderedDict


class FeaturesPlotter:
    def __init__(self, inputfile, clidfile, featuresfile, nFeatures=None, clusterlogfile=None, begin=0, end=-1):
        self.inputfile = inputfile
        self.clidfile = clidfile
        self.featuresfile = featuresfile
        self.clusterlogfile = clusterlogfile
        self.begin = begin
        self.end = end
        self.nFeatures = nFeatures
        
        self.clids = None
        self.clidsvsframes = None
        self.clidstime = None
        self.features = None
        self.featurestime = None
        
        self.metaplotdata = None
        
    def plot_features(self, outfile, topmargin=0.1, legendcols=5, width=12, height=10, fontsize=12, dpi=300):
        
        if self.features is None:
            self.read_featuresfile()
        
        if self.clids is None:
            self.read_clidfile()
            
        self.checkFeaturesVsClidTime()
            
        if self.metaplotdata is None:
            self.read_input_file()
            
        for mpltdata in self.metaplotdata:
            if (   (mpltdata['x'] > len(self.features))
                or (mpltdata['y'] > len(self.features))):
                raise ValueError('The requested feature is not found or read from features file.')
            
        if len(self.metaplotdata) % 2 != 0:
            nrows = int(len(self.metaplotdata)/2) + 1
        else:
            nrows = int(len(self.metaplotdata)/2)
        
        # TODO: make input colormap
        cmap_list = [(0, '#c2c0c1'), (0.25, '#46a6e4'), (0.75, '#c01755'), (1.0, '#000000')]
        cmap = mpl.colors.LinearSegmentedColormap.from_list('dummy', cmap_list, N=len(self.clids))
        normalize = mpl.colors.Normalize(vmin=self.clids[0], vmax=self.clids[-1])
        colors = OrderedDict()
        for clid in self.clids:
            colors[clid] = cmap(normalize(clid))

        fig = plt.figure(figsize=(height, width))
        axes = []
        #fig.subplots_adjust(wspace=0.3, hspace=0.5)
        mpl.rcParams['font.size'] = fontsize
        handles, legend_labels = None, None
        
        
        for i in range(len(self.metaplotdata)):
            ax = fig.add_subplot(nrows,2,i+1)
            axes.append(ax)
            xid = self.metaplotdata[i]['x']
            yid = self.metaplotdata[i]['y']
            for clid in self.clids:
                x = self.features[xid-1][self.clidsvsframes[clid]]
                y = self.features[yid-1][self.clidsvsframes[clid]]
                ax.scatter(x, y, s=0.5, label=str(clid), color=colors[clid])
                ax.set_aspect('equal', 'datalim')
                
                if 'xlabel' in self.metaplotdata[i]:
                    ax.set_xlabel(self.metaplotdata[i]['xlabel'])
                else:
                    ax.set_xlabel('feature-{0}'.format(xid))

                if 'ylabel' in self.metaplotdata[i]:
                    ax.set_ylabel(self.metaplotdata[i]['ylabel'])
                else:
                    ax.set_ylabel('feature-{0}'.format(yid))
                    
        handles, legend_labels = axes[0].get_legend_handles_labels()

        fig.legend(handles, legend_labels, ncol=legendcols, loc='upper center',scatterpoints=4,markerscale=6)
        fig.set_tight_layout(tight={'rect':(None,None,None,1-topmargin)})
        plt.savefig(outfile, dpi=dpi)

        
    def read_input_file(self):
        metaplotdata = []
        fin = open(self.inputfile, 'r')
        
        for line in fin:
            line = line.lstrip().rstrip()
            if not line.strip():
                continue
            
            temp = re.split(',', line)
            
            if len(temp) < 2:
                raise ValueError('Two features are required to plot')
            
            tobeplot = OrderedDict()
            tobeplot['x'] = int(temp[0])
            tobeplot['y'] = int(temp[1])
            if len(temp) > 2:
                tobeplot['xlabel'] = r'{0}'.format(temp[2])
            if len(temp) > 3:
                tobeplot['ylabel'] = r'{0}'.format(temp[3])
                
            metaplotdata.append(tobeplot)
            
        self.metaplotdata = metaplotdata
        
        
    def read_featuresfile(self):
        fin = open(self.featuresfile, 'r')
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
                if ecount == 1:
                    self.featurestime = np.asarray(time)
                ppc = []
                if self.nFeatures is not None:
                    if ecount == self.nFeatures:
                        break
                continue

            if re.search('^#|^@', line) is None:
                temp = re.split('\s+', line)
                t = float(temp[0])
                
                # Do not append when end time reached
                if (self.end != -1):
                    if (t >= self.begin) and (t <= self.end):
                        if ecount == 0:
                            time.append(t)
                        ppc.append(float(temp[1]))
                else:
                    if (t >= self.begin):
                        if ecount == 0:
                            time.append(t)
                        ppc.append(float(temp[1]))
                    

        fin.close()
        self.features = np.asarray(features)
        self.nFeatures = ecount
        
    def read_clidfile(self):
        fin = open(self.clidfile, 'r')
        clids_data = OrderedDict()
        time = []

        index = 0
        for line in fin:
            line = line.rstrip().lstrip()
            if not line:
                continue

            if re.search('^#|^@', line) is not None:
                        continue

            temp = re.split('\s+', line)

            t = float(temp[0])
            time.append(t)
            clid = int(temp[1])

            if t >= self.begin:
                if clid not in clids_data:
                    clids_data[clid] = [ index ]
                else:
                    clids_data[clid] = clids_data[clid] + [ index ]

                index = index + 1

            if (self.end != -1):
                if (t >= self.end):
                    break


        fin.close()
        
        clids = sorted(clids_data.keys())
        
        self.clids = clids
        self.clidstime = np.asarray(time)
        self.clidsvsframes = clids_data
        
    def checkFeaturesVsClidTime(self):
        if self.clidstime.shape[0] != self.featurestime.shape[0]:
            raise AssertionError("time in features file does not match with time in cluster-id file.")
        
        if not np.allclose(self.clidstime, self.featurestime):
            raise AssertionError("time in features file does not match with time in cluster-id file.")
                
