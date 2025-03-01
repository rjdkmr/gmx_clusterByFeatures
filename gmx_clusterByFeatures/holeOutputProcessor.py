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
from numpy import ma
import sys
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from collections import OrderedDict
import math
from sklearn.decomposition import PCA

class HoleOutputProcessor:

    def __init__(self, filename, axis='Z', endrad=5, xmin=None, xmax=None, gap=1, begin=0, end=-1, dataOccupancy=90):
        ''' Class to read, process and visualzie hole output

        Parameters
        ----------
        filename : str
            Input file obtained from ``hole`` sub-command.
        axis : str
            Principal axis parallel to the channel or cavity.
        endrad : float
            If radius is larger than this value, this value will not considered
            for average calculation and features output. This value might be equal or
            less than ``-endrad`` value supplied with ``hole`` sub-command.
        xmin : float
            Minimum value of axis point after which radius value will be considered.
            If not supplied, minimum axis value will be extracted from input file.
        xmax : float
            Maximum value of axis point before which radius value will be considered.
            If not supplied, maximum axis value will be extracted from input file.
        gap : float
            Gap between axis points in Ångstroms. It should be either equal to or larger than
           ``-sample value supplied with hole sub-command.
        begin : float
            First frame in time to read from the input file
        end : float
            Last frame in time to read from the input file. If its ``end = -1``, All frames till
            the end will be read.
        dataOccupancy : float
            Percentage of radius-data occupancy for axis-points. If an axis-point has radius-data
            less than this percentage of frame, the axis-point will not be considered for average
            calculation and features output.

            This is critical for axis-points, which are at the
            opening of channel/cavity. In several frames, radius-value could be missing and therefore,
            dataOccupancy threshold could be used to discard those axis points with lots of missing
            radius values.

        '''

        self.filename = filename
        self.paxis = axis
        self.xmin = xmin
        self.xmax = xmax
        self.gap = gap
        self.begin = begin
        self.end = end
        self.dataOccupancy = dataOccupancy
        self.endrad = endrad

        self.axis2index = { 'X':0, 'Y':1 , 'Z':2}
        self.flag_radius = 99999

        self.max_axis_value = None
        self.min_axis_value = None
        self.output_range = None
        self.frame_number = 0

        self.axis_value = None
        self.radius = OrderedDict()
        self.rawResidues = OrderedDict()
        self.residuesByAxis = None
        self.residue_list = None
        self.fullDataAxis = None
        self.average = None
        self.error = None
        self.time = []


        if self.end != -1:
            if self.end < self.begin:
                raise AssertionError('Begin time cannot be larger than end time.')

        if self.dataOccupancy > 100 or self.dataOccupancy < 0:
            raise ValueError('Data occupancy thershold should be between 0 to 100.')

        if self.paxis not in self.axis2index.keys():
            raise ValueError('{0} is requested for axis. However, only X, Y and Z axis are accepted'.format(self.paxis))


    def plot_by_cluster(self, cluster_file, outfile, csvfile=None, stdbar=False, discard_lasts=0, ymin=None, ymax=None, width=6, height=4, fontsize=12, rightmargin=0.15, legendcols=1, dpi=300):
        '''Plot radius by Clusters
        It can be used to plot radius of cavity/channel for clusters seperately.
        It reads radius file from "hole" and cluster-id file  from "cluster",
        and extract radius of each cluster separately and plot them in one plot.
        This plot could be extremely useful to compare radius along the
        channel/cavity in all clusters.

        Parameters
        ----------
        cluster_file : str
            Input file containing cluster-id as a function of time. The number of frames in this
            file should be same as in radius file.

        outfile : str
            Output plot file. The radius as a function of axis-points is plotted in this file.

        csvfile : str
            Output csv file. The radius as a function of axis-points in csv formatted file. This
            file can be read in external data-plotting program.

        stdbar : bool
            If it is ``True``, standard deviation will be shown as an error-bar in the plot.

        discard_lasts : int
            Number of smallest clusters to discard from the plotting. It can be useful to filter
            out few smallest clusters because these may contain small number of frames.

        ymin : float
            Minimum value at Y-axis.
            If not supplied minimum value from data will be used. It can be useful to minimum and
            maximum values of Y-axis when several plots are compared together.

        ymax : float
            Maximum value at Y-axis.
            If not supplied maximum value from data will be used. It can be useful to minimum and
            maximum values of Y-axis when several plots are compared together.

        width : int
            Width of the plot

        height : int
            Height of the plot

        fontsize : int
            Font size in the plot

        rightmargin : float
            Margin at right side of the plots. If legends overflow into the plot area, margin can
            be increased to fit the legend.

        legendcols : int
            If legend overflow into the plot area, legends can be made of more than
            one column to accommodate all legends.

        '''

        # Read hole data if not read previously
        if self.axis_value is None:
            self.read_hole_data()

        # Calculate average and standard-deviation to mark the axis point containing missing data
        if self.average is None:
            self.calculate_average()

        clids_data = self.read_clid(cluster_file)
        clids = sorted(clids_data.keys())
        if discard_lasts > 0:
            clids = clids[:-1*discard_lasts]

        averages = OrderedDict()
        stddevs = OrderedDict()
        for clid in clids:
            means, sd = [], []
            for key in list(map(self._float2str, self.axis_value)):
                #if float(key) in self.fullDataAxis:
                means.append((self.radius[key])[clids_data[clid]].mean())
                sd.append((self.radius[key])[clids_data[clid]].std())
            averages[clid] = means
            stddevs[clid] = sd

        # Write data to a csv file
        if csvfile is not None:
            with open(csvfile, 'w', newline='') as csvfile:
                csvwriter = csv.writer(csvfile, dialect='excel')

                # Write header
                row = ['Axis']
                for clid in clids:
                    row.append(clid)
                csvwriter.writerow(row)

                # Write data
                for i in range(len(self.axis_value)):
                    row = [self.axis_value[i]]
                    for clid in clids:
                        row.append(averages[clid][i])
                    csvwriter.writerow(row)

        # TODO: make input colormap
        cmap_list = [(0, '#c2c0c1'), (0.25, '#46a6e4'), (0.75, '#c01755'), (1.0, '#000000')]
        cmap = mpl.colors.LinearSegmentedColormap.from_list('dummy', cmap_list, N=len(clids))
        normalize = mpl.colors.Normalize(vmin=clids[0], vmax=clids[-1])
        colors = OrderedDict()
        for clid in clids:
            colors[clid] = cmap(normalize(clid))

        mpl.rcParams['font.size'] = fontsize
        fig, ax1 = plt.subplots(1, sharex=True, gridspec_kw={'hspace': 0}, figsize=(width, height))
        #fig.subplots_adjust(top=0.9, bottom=0.1)


        # Plot for radius
        color=['#000000', '#BDBDBD']
        for clid in clids:
            if stdbar:
                eb = ax1.errorbar(self.axis_value, averages[clid], yerr=stddevs[clid], color=colors[clid], marker='o', markersize=3.0, label=str(clid))
                for line in eb.lines[2]:
                    line.set_alpha(0.3)
            else:
                ax1.plot(self.axis_value, averages[clid], label=str(clid), marker='o', ms=3, color=colors[clid])
        handles, legend_labels = ax1.get_legend_handles_labels()
        ax1.set_ylabel(r'Radius ($\AA$)')
        ax1.set_xlabel(r'Axis ($\AA$)')

        # Set ylimits if given
        ylims = ax1.get_ylim()
        if ymin is not None:
            ax1.set_ylim(ymin, ylims[1])
            ylims = [ymin, ylims[1]]
        if ymax is not None:
            ax1.set_ylim(ylims[0], ymax)
            ylims = [ylims[0], ymax]


        legend = fig.legend(handles, legend_labels, ncol=legendcols, loc='right',scatterpoints=5, markerscale=3)
        fig.set_tight_layout(tight={'rect':(None,None,1-rightmargin, None)})
        fig.savefig(outfile, dpi=dpi, bbox_extra_artists=(legend, ), bbox_inches='tight')

    def write_features(self, outfile, pca=None):
        ''' Write radius as features for clustering

        .. :note: If a frame does not contain radius value, zero is written.

        Parameters
        ----------
        outfile : str
            Name of output file. Output file containing radius as function of
            time at each axis points. This file can be used as features file
            for clustering. This file can be also used to plot radius vs time
            with external plotting program.

            The file name should end with xvg extension, which is
            recognized by ``cluster`` command.
            
        pca : int
            Number of eigenvectors to be considered for the features.
            In place for taking radius as features, this option enable PCA of radii
            and the resultant projections on eigenvectors can be used as features.
            


        '''
        toBePrinted = True
        msg =  "WARNING: radius value for the axis point was not calculated by hole. Zero will be written as feature value. "
        msg += "However, large number of missing values will lead to wrong clustering. Therefore, please try to minimize or "
        msg += "eliminate the missing values by changing axis-point range using xmin and xmax option and/or dataOccupancy option."
        
        # Read hole data if not read previously
        if self.axis_value is None:
            self.read_hole_data()

        # Calculate average and standard-deviation to mark the axis point containing missing data
        if self.average is None:
            self.calculate_average()


        if pca is not None:
            features = []
            for i in range(len(self.axis_value)):
                key = self._float2str(self.axis_value[i])
                if float(key) in self.axis_value:
                    x = self.radius[key]
                    if isinstance(x, ma.MaskedArray):
                        if toBePrinted:
                            print(msg)
                            toBePrinted = False
                        features.append(x.filled(0))
                    else:
                        features.append(x)
            features = np.asarray(features)
            print(features.shape)
            
            # Use Singular Value Decomposition
            pcaobj = PCA(pca)
            projs = pcaobj.fit_transform(features.T)
            projs = projs.T
            
            fout = open(outfile, 'w')
            for i in range(projs.shape[0]):
                for j in range(projs.shape[1]):
                    fout.write('{0:18.3f} {1:18.3f}\n'.format(self.time[j], projs[i][j]))
                fout.write('\n & \n')
            fout.close()
            
        else:
            fout = open(outfile, 'w')
            for key in list(map(self._float2str, self.axis_value)):
                # Axis point with any missing data are not written
                if float(key) in self.axis_value:
                    for i in range(len(self.radius[key])):
                        if self.radius[key][i] is ma.masked:
                            if toBePrinted:
                                print(msg)
                                toBePrinted = False
                            fout.write('{0:18.3f} {1:18.3f}\n'.format(self.time[i], 0))
                        else:
                            fout.write('{0:18.3f} {1:18.3f}\n'.format(self.time[i], self.radius[key][i]))
                    fout.write('\n & \n')
            fout.close()

    def plot_radius_residues(self, outfile, csvfile=None,  violinplot=False, residue_frequency=50, ymin=None, ymax=None, width=6, height=6, fontsize=18, rlabelsize=10, dpi=300):
        ''' Plot radius and residues as a function of axis points

        Parameters
        ----------
        outfile : str
            Output plot file

        csvfile : str
            Output csv file. The radius as a function of axis-points in csv formatted file. This
            file can be read in external data-plotting program.

        violinplot : bool
            In place of normal line-plot, it plots radius distribution as violins. It is useful
            because this plot gives distribution of radius values over entire trajectory for each
            axis-points

        residue_frequency : float
            Frequency (%) of residue occurrence during the simulations at a given axis points.
            If frequency is less than this threshold, it will not considered for plotting.

        ymin : float
            Minimum value at Y-axis.
            If not supplied minimum value from data will be used. It can be useful to minimum and
            maximum values of Y-axis when several plots are compared together.

        ymax : float
            Maximum value at Y-axis.
            If not supplied maximum value from data will be used. It can be useful to minimum and
            maximum values of Y-axis when several plots are compared together.


        width : int
            Width of the plot

        height : int
            Height of the plot

        fontsize : int
            Font size in the plot

        rlabelsize : float
            Fontsize of residue label along Y-axis.

        '''

        if self.average is None:
            self.calculate_average()

        if csvfile is not None:
            self.average2csvfile(csvfile)

        mpl.rcParams['font.size'] = fontsize
        fig = plt.figure(figsize=(width, height))
        fig, (ax1, ax2) = plt.subplots(2, sharex=True, gridspec_kw={'hspace': 0}, figsize=(width, height))        
        fig.subplots_adjust(hspace=0.0, wspace=0)

        # Plot for radius
        if violinplot:
            self._violinplot(ax1)
        else:
            color=['#000000', '#BDBDBD']
            ax1.errorbar(self.axis_value, self.average.values(), yerr=self.error.values(), ecolor=color[1], elinewidth=2.5, color=color[0], lw=1.5, marker='o', mfc=color[0], mew=0, ms=1)
            ax1.plot(self.axis_value, self.average.values(), '-r')
        ax1.set_ylabel(r'Radius ($\AA$)')
        ax1.grid(which='major', axis='x', linewidth=0.5, alpha=0.5)
        ax1.tick_params(axis='x', direction='in', labelbottom=False, top=True)

        # Set ylimits if given
        ylims = ax1.get_ylim()
        if ymin is not None:
            ax1.set_ylim(ymin, ylims[1])
            ylims = [ymin, ylims[1]]
        if ymax is not None:
            ax1.set_ylim(ylims[0], ymax)
            ylims = [ylims[0], ymax]


        # Plot for residues
        if self.residuesByAxis is None:
            self.processResiduesData()

        plot_coords, ylabels, median = [], [], []
        for res in self.residuesByAxis:
            if (len(self.residuesByAxis[res])/self.frame_number) < (residue_frequency/100):
                continue
            plot_coords.append(self.residuesByAxis[res])
            median.append(np.median(self.residuesByAxis[res]))
            ylabels.append(res)

        idx = np.argsort(median)
        for i in range(len(median)):
            violins = ax2.violinplot(plot_coords[idx[i]], positions=[i], vert=False, showmedians=True, showextrema=False)
            for pc in violins['bodies']:
                pc.set_edgecolor('black')
                pc.set_linewidth(0.5)
                pc.set_alpha(0.5)

        ax2.set_yticks(range(0,len(median),1))
        ax2.set_yticklabels(np.asarray(ylabels)[idx], fontsize=rlabelsize)
        ax2.grid(which='major', axis='both', linewidth=0.5, alpha=0.5)
        ax2.set_ylim(-1, len(median))
        ax2.set_xlabel(r'axis ($\AA$)')
        ax2.tick_params(axis='x', direction='inout')
        ax2.tick_params(axis='x', direction='in', top=True)

        fig.tight_layout()
        plt.savefig(outfile, dpi=dpi)

    def plot_radius(self, outfile, csvfile=None, violinplot=False, ymin=None, ymax=None, width=6, height=4, fontsize=18, dpi=300):
        ''' Plot radius and residues as a function of axis points

        Parameters
        ----------
        outfile : str
            Output plot file

        csvfile : str
            Output csv file. The radius as a function of axis-points in csv formatted file. This
            file can be read in external data-plotting program.

        violinplot : bool
            In place of normal line-plot, it plots radius distribution as violins. It is useful
            because this plot gives distribution of radius values over entire trajectory for each
            axis-points

        ymin : float
            Minimum value at Y-axis.
            If not supplied minimum value from data will be used. It can be useful to minimum and
            maximum values of Y-axis when several plots are compared together.

        ymax : float
            Maximum value at Y-axis.
            If not supplied maximum value from data will be used. It can be useful to minimum and
            maximum values of Y-axis when several plots are compared together.


        width : int
            Width of the plot

        height : int
            Height of the plot

        fontsize : int
            Font size in the plot

        '''

        if self.average is None:
            self.calculate_average()

        if csvfile is not None:
            self.average2csvfile(csvfile)

        mpl.rcParams['font.size'] = fontsize
        fig = plt.figure(figsize=(width, height))
        fig, (ax1, ax2) = plt.subplots(2, sharex=True, gridspec_kw={'hspace': 0}, figsize=(width, height))
        fig.subplots_adjust(hspace=0.0, wspace=0)


        # Plot for radius
        if violinplot:
            self._violinplot(ax1)
        else:
            color=['#000000', '#BDBDBD']
            ax1.errorbar(self.axis_value, self.average.values(), yerr=self.error.values(), ecolor=color[1], elinewidth=2.5, color=color[0], lw=1.5, marker='o', mfc=color[0], mew=0, ms=1)
            ax1.plot(self.axis_value, self.average.values(), '-r')
        ax1.set_ylabel(r'Radius ($\AA$)')
        ax1.grid(which='major', axis='x', linewidth=0.5, alpha=0.5)
        ax1.set_xlabel(r'Axis ($\AA$)')

        # Set ylimits if given
        ylims = ax1.get_ylim()
        if ymin is not None:
            ax1.set_ylim(ymin, ylims[1])
            ylims = [ymin, ylims[1]]
        if ymax is not None:
            ax1.set_ylim(ylims[0], ymax)
            ylims = [ylims[0], ymax]

        fig.tight_layout()
        plt.savefig(outfile, dpi=dpi)

    def calculate_average(self):
        '''Calculate average and standard deviation of radius along the axis
        '''

        # Read hole data if not read previously
        if self.axis_value is None:
            self.read_hole_data()

        # Remove other axis which does not found in all frames
        radius = OrderedDict()
        residues = OrderedDict()
        self.fullDataAxis = []
        keys_old = list(sorted(self.axis_value))
        keys_new = []
        for key in keys_old:
            key = self._float2str(key)
            radius_npy = np.asarray(self.radius[key])
            # Use real numpy array else change to masked numpy array
            if self.flag_radius not in self.radius[key]:
                keys_new.append(float(key))
                radius[key] = radius_npy
                residues[key] = self.rawResidues[key]
                self.fullDataAxis.append(float(key))
            else:
                data_found = (radius_npy != self.flag_radius)
                idx = np.nonzero(data_found)
                if data_found.sum()/self.frame_number >= self.dataOccupancy/100:
                    keys_new.append(float(key))
                    radius[key] = ma.masked_values(radius_npy, self.flag_radius)
                    residues[key] = np.asarray(self.rawResidues[key])

        self.axis_value = np.asarray(keys_new)
        self.radius = radius
        self.rawResidues = residues

        print('After removing axis points with any missing data:')
        print('------------------------------------------------------------')
        print('#{0:>12} {1:>18} {2:>18}'.format('Axis', 'Mean radius',  'Std. deviation'))
        print('------------------------------------------------------------')
        self.average = OrderedDict()
        self.error = OrderedDict()
        for key in list(map(self._float2str, self.axis_value)):
            #print(type(self.radius[key]))
            self.average[key] = self.radius[key].mean()
            self.error[key]  = self.radius[key].std()
            print('{0:12.3f} {1:16.3f} {2:>16.3f}'.format(float(key), self.average[key], self.error[key]))
        print('------------------------------------------------------------')

    def average2csvfile(self, csvfile):
        '''Write average of radius to a CSV output file

        Parameters
        ----------

        csvfile : str
            Name of output CSV file
        '''
        if self.average is None:
            self.calculate_average()

        # Write data to a csv file
        with open(csvfile, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile, dialect='excel')

            # Write header
            row = ['Axis', 'Radius', 'std-deviation']
            csvwriter.writerow(row)

            # Write data
            for key in list(map(self._float2str, self.axis_value)):
                row = [float(key), self.average[key], self.error[key]]
                csvwriter.writerow(row)


    def processResiduesData(self):
        self.residuesByAxis = OrderedDict()
        self.residue_list = []
        for key in list(map(self._float2str, self.axis_value)):
            for i in range(len(self.rawResidues[key])):
                temp = re.split(',', self.rawResidues[key][i])
                for tres in temp:
                    if tres.rstrip().lstrip():
                        if tres == 'dummy':
                            continue
                        if tres not in self.residue_list:
                            self.residue_list.append(tres)
                            self.residuesByAxis[tres] = []
                        self.residuesByAxis[tres].append(float(key))

    def _violinplot(self, ax):
        values = []
        positions = []

        for key in list(map(self._float2str, self.axis_value)):
            # Axis point with any missing data are not written
            if float(key) in self.axis_value:
                if hasattr(self.radius[key], 'compressed'):
                    values.append(self.radius[key].compressed())
                else:
                    values.append(self.radius[key])
                positions.append(float(key))

        violins = ax.violinplot(values, positions=positions, showmedians=True, showextrema=False)
        for pc in violins['bodies']:
            pc.set_edgecolor('black')
            pc.set_linewidth(0.5)
            pc.set_alpha(0.5)

    def _float2str(self, x):
        return str(np.round(x, 3))

    def _append_block_data(self, block_data, time):

        block_data_transpose = np.array(block_data).T
        input_axis_string = block_data_transpose[self.axis2index[self.paxis]]
        input_axis = np.asarray(list(map(float, block_data_transpose[self.axis2index[self.paxis]])))
        input_radius = list(map(float,block_data_transpose[3]))
        input_residues = block_data_transpose[4]
        print_new_output_range = False

        if np.round(abs(input_axis[1]-input_axis[0]), 2) > self.gap:
            raise ValueError('Gap {0} between slabs in hole output is larger than the input gap {1}. '
                             'Input gap should be either equal aur larger than gap between slabs in hole output'.format(abs(input_axis[1]-input_axis[0]), self.gap))

        if self.max_axis_value is None:
            self.max_axis_value = np.amax(input_axis)
        else:
            self.max_axis_value = max(np.amax(input_axis), self.max_axis_value)

        if self.min_axis_value is None:
            self.min_axis_value = np.amin(input_axis)
        else:
            self.min_axis_value = min(np.amax(input_axis), self.min_axis_value)

        # Check if axis value is outside of input range
        '''
        if self.xmin is not None:
            if self.xmin < self.min_axis_value:
                raise ValueError ("Input minimum axis value {0} is less than minimum axis value {1} at time {2}".format(self.xmin, self.min_axis_value, time))
        if self.xmax is not None:
            if self.xmax > self.max_axis_value:
                raise ValueError ("Input maximum axis value {0} is less than maximum axis value {1} at time {2}".format(self.xmax, self.max_axis_value, time))
        '''

        # Initialize output range for axis
        if self.output_range is None:
            self.output_range = []
            if self.xmin is not None:
                self.output_range.append(self.xmin)
            else:
                self.output_range.append(self.min_axis_value-(self.min_axis_value%self.gap))

            if self.xmax is not None:
                self.output_range.append(self.xmax)
            else:
                self.output_range.append(self.max_axis_value-(self.max_axis_value%self.gap)+self.gap)

            sys.stdout.write('\nMaximum axis range is: ({0:.3f} to {1:.3f})\n'.format(self.output_range[0], self.output_range[1]))
            sys.stdout.flush()

        # Initialize axis value
        if self.axis_value is None:
            self.axis_value = np.arange(self.output_range[0], self.output_range[1], self.gap)

        # Change axis-values list if input range is not given
        if (self.min_axis_value < self.output_range[0]) and (self.xmin is None):
            tmin = self.output_range[0]
            while(tmin >= self.min_axis_value):
                tmin = tmin - self.gap
            self.axis_value = np.hstack((np.arange(tmin, self.output_range[0], self.gap), self.axis_value))
            self.output_range[0] = tmin
            print_new_output_range = True

        # Change axis-values list if input range is not given
        if (self.max_axis_value > self.output_range[1]) and (self.xmax is None):
            tmax = self.output_range[1]
            while(tmax <= self.max_axis_value):
                tmax = tmax + self.gap
            self.axis_value = np.hstack((self.axis_value, np.arange(self.output_range[1], tmax, self.gap)))
            self.output_range[1] = tmax
   
        if print_new_output_range:
            sys.stdout.write('\nNew Maximum axis range at time {0} is: ({1:.3f} to {2:.3f})\n'.format(time, self.output_range[0], self.output_range[1]))
            sys.stdout.flush()


        # Initialize temorary dictionary of lists
        radius_list, residues_list = OrderedDict(), OrderedDict()
        for value in list(map(self._float2str, self.axis_value)):
            radius_list[value] = []
            residues_list[value] = ''

        # Read from block data and append in temorary dictionary of lists
        for i in np.argsort(input_axis):
            min_idx = np.abs((input_axis[i] - self.axis_value) - 0).argmin()
            if abs(input_axis[i] - self.axis_value[min_idx]) <= self.gap:
                x = self._float2str(self.axis_value[min_idx])

                # If a new axis point is found, all previous frame should be filled with this value
                if x not in self.radius:
                    self.radius[x] = [self.flag_radius] * self.frame_number
                if x not in self.rawResidues:
                    self.rawResidues[x] = ['dummy'] * self.frame_number

                # Append values in temporary list
                if input_radius[i] < self.endrad:
                    radius_list[x].append(input_radius[i])
                else:
                    radius_list[x].append(self.flag_radius)
                residues_list[x] += ',' + input_residues[i]

        # Append in final list after merging results
        for value in list(map(self._float2str, self.axis_value)):
            if radius_list[value]:
                # Remove dummy radius value if present in any bin
                if self.flag_radius in radius_list[value]:
                    t = np.asarray(radius_list[value])
                    yes_radius = (t != self.flag_radius)
                    if yes_radius.sum() != 0:
                        self.radius[value].append(np.amax(t[yes_radius]))
                    else:
                        self.radius[value].append(self.flag_radius)
                else:
                    self.radius[value].append(np.amax(radius_list[value]))
                self.rawResidues[value].append(residues_list[value])
            else:
                # Sometime input range is outside in first frame, therefore, here again initialized 
                if value not in self.radius:
                    self.radius[value] = [self.flag_radius] * self.frame_number
                if value not in self.rawResidues:
                    self.rawResidues[value] = ['dummy'] * self.frame_number
                self.radius[value].append(self.flag_radius)
                self.rawResidues[value].append('dummy')

        self.time.append(time)

    def _print_frame_number(self, time):
        if(    ((self.frame_number<100) and (self.frame_number%10==0))
           or  ((self.frame_number<1000) and (self.frame_number%100==0))
           or  ((self.frame_number<10000) and (self.frame_number%1000==0))
           or  ((self.frame_number<100000) and (self.frame_number%10000==0))
           or  ((self.frame_number<1000000) and (self.frame_number%100000==0))
           ):
            sys.stdout.write("\rReading frame %d at time %12.3f" %(self.frame_number, time))
            sys.stdout.flush()

    def read_hole_data(self):
        sys.stdout.write("\nReading file : %s\n" % self.filename)
        sys.stdout.flush()

        infile = open(self.filename,'r')
        block =[]
        frame_number = 0
        time = 0

        for line in infile:
            #Removing last new line charecter
            line = line.rstrip('\n')

            #Skipping blank/empty line
            if not line.strip():
                continue

            #Getting Time tag and time => Starting of new frame
            if(re.match('# Time',line)!=None):
                self._print_frame_number(time)

                if(len(block)>0):
                    if time >= self.begin:
                        self._append_block_data(block, time)
                        self.frame_number += 1

                    if (self.end != -1):
                        if (time >= self.end):
                            break

                block = []
                time = float( re.split('=', line)[1])
                continue

            #Skipping other lines starting with '#' tag'
            if(re.match('#',line)!=None):
                continue

            block.append(line.split())


        #For last frame
        if self.end == -1:
            self._append_block_data(block, time)
            self.frame_number += 1

        sys.stdout.write("\nFinishid reading.... Total number of frame read =  %d\n" % self.frame_number)
        sys.stdout.flush()
        infile.close()

    def read_clid(self, filename):
        fin = open(filename, 'r')
        data = OrderedDict()

        index = 0
        for line in fin:
            line = line.rstrip().lstrip()
            if not line:
                continue

            if re.search('^#|^@', line) is not None:
                        continue

            temp = re.split('\s+', line)

            time = float(temp[0])
            clid = int(temp[1])

            if time >= self.begin:
                if clid not in data:
                    data[clid] = [ index ]
                else:
                    data[clid] = data[clid] + [ index ]

                index = index + 1

            if (self.end != -1):
                if (time >= self.end):
                    break


        fin.close()

        return data


def get_error(time, x, sets, err_type='block', tool='gmx analyze'):
    """To estimate error using block averaging method
    .. warning::
            To calculate errors by using ``error = 'acf'`` or ``error = 'block'``,
            GROMACS tool ``g_analyze`` or ``gmx analyze`` should be present in ``$PATH``.
    Parameters
    ----------
    time : 1D list or array
        :attr:`DNA.time`
    x   : 2D list or array
        Shape of (nset, nframe); where *nset* is number of set and *nframe* is
        total number of frames. *nframe* should be equal to length of time
        list/array
    sets : int
        Number of sets (*nset*)
    err_type : str
        Error estimation by autocorrelation method ``err_type='acf'`` or
        block averaging method ``err_type='block'``
    tool : str
        GROMACS tool to calculate error. In older versions it is `g_analyze` while in
        newer versions (above 2016) it is `gmx analyze`.
    Returns
    -------
    error : 1D array
        Of length = number of sets (*nset*)
    """
    for i in range(sets):
        if (len(time) != len(x[i])):
            raise ValueError('\nError: number of frame in time {0} mismatched with {1} of x[{2}]!!\n' .format(
                len(time), len(x[i]), i))

    if not((err_type == 'block') or (err_type == 'acf')):
        print('\nWarning: Method {0} is not implemented. Switching to \'acf\'.\n' .format(
            err_type))
        err_type = 'acf'

    error = []
    char_set = string.ascii_lowercase
    name = ''.join(random.sample(string.ascii_lowercase, 10))

    filename = name + '.xvg'
    eefile = 'ee_' + name + '.xvg'
    acfile = 'acf_' + name + '.xvg'

    fout = open(filename, 'w')
    for i in range(len(time)):
        fout.write('{0}' .format(time[i]))
        for j in range(sets):
            fout.write('    {0}' .format(x[j][i]))
        fout.write("\n")
    fout.close()

    command = '{0} -f {1} -ee {2} -ac {3} -fitfn exp' .format(tool, filename, eefile, acfile)

    try:
        p = sub.Popen(command.split(), stdout=sub.PIPE, stderr=sub.PIPE, universal_newlines=True)
        out, outputerror = p.communicate()
    except Exception as e:
        os.remove(filename)
        raise e

    lines = out.split('\n')

    if (err_type == 'block'):
        for line in lines:
            if(re.match('Set', line)):
                temp = line.split()
                error.append(float(temp[3]))

    if (err_type == 'acf'):
        acf_time = []
        for line in lines:
            if re.match('COR: Correlation time', line) is not None:
                temp = line.split('=')
                acf_time.append(abs(float(temp[1].split()[0])))

        total_time = float(time[-1]) - float(time[0])
        dt = total_time / len(time)
        for i in range(sets):
            if(acf_time[i] >= dt):
                n_indp = total_time / acf_time[i]
                tmp_err = np.std(x[i]) / np.sqrt(n_indp)
            else:
                tmp_err = np.std(x[i]) / np.sqrt(len(time))
            error.append(tmp_err)

    os.remove(filename)
    os.remove(eefile)
    os.remove(acfile)
    if os.path.isfile('fitlog.log'):
        os.remove('fitlog.log')

    return np.array(error)
