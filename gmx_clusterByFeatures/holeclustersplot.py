"""
This file is part of gmx_clusterByFeatures

Author: Rajendra Kumar

Copyright (C) 2014-2025  Rajendra Kumar

gmx_clusterByFeatures is a free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

"""

import re
import sys
import argparse
import os
import numpy as np

from . import holeOutputProcessor

def get_output_formats(fmt=False):
    fig = holeOutputProcessor.plt.figure()
    output_formats_dict = fig.canvas.get_supported_filetypes()
    del fig
    
    if fmt:
        result = ''
        for key in output_formats_dict:
            result += key + ' : ' + output_formats_dict[key] + '\n'
        return result
    else:
        return output_formats_dict


description=\
"""DESCRIPTION
===========
It can be used to plot radius of cavity/channel for clusters separately.
It reads radius file from "hole" and cluster-id file  from "cluster",
and extract radius of each cluster separately and plot them in one plot.
This plot could be extremely useful to compare radius along the 
channel/cavity in all clusters.

"""

inputFileHelp=\
"""Name of input file.
Name of input radius file. Radius file should be obtained from "hole" as an 
output file.

"""

clusterFileHelp=\
"""Input file containing cluster-id as a function of time.
The number of frames in this file should be same as in input radius file.

"""

outputFileHelp=\
"""Name of output file.
Name of the output plot file. The extension will be used to determine
the output format.

Following output formats are available:
{0}

""".format(get_output_formats(fmt=True))

outCsvFileHelp=\
"""Output csv file.
The radius as a function of axis-points in csv formatted file. This
file can be read in external data-plotting program.

"""

xminHelp=\
"""Minimum value of axis-point
Minimum value of axis point after which radius value will be considered for plot.

If not supplied, minimum axis value will be extracted from input radius file.

"""

xmaxHelp=\
"""Maximum value of axis-point
Maximum value of axis point after which radius value will be discarded from plot.

If not supplied, maximum axis value will be extracted from input radius file.

"""

yminHelp=\
"""Minimum value at Y-axis.
If not supplied minimum value from data will be used. It can be useful to minimum and 
maximum values of Y-axis when several plots are compared together.

"""
            
ymaxHelp=\
"""Maximum value at Y-axis.
If not supplied maximum value from data will be used. It can be useful to minimum and 
maximum values of Y-axis when several plots are compared together.

"""


endradHelp=\
"""End/Opening radius
If radius is larger than this value, radius will not considered 
for average calculation and features output. This option value might be equal or
less than "-endrad" value supplied with "hole" sub-command.

"""

gapHelp=\
"""Gap between axis-points in Ångstroms
It should be either equal to or larger than "-sample" value supplied 
with "hole" sub-command.
"""

endHelp=\
"""Last frame in time to read from the input file.
If its ``end = -1``, All frames till the end will be read.

"""

dataOccupancyHelp=\
"""Percentage of radius-data occupancy for axis-points.
If an axis-point has radius-data
less than this percentage of frame, the axis-point will not be considered 
for average calculation and features output.

This is critical for axis-points, which are at the 
opening of channel/cavity. In several frames, radius-value could be missing
and therefore, dataOccupancy threshold could be used to discard those axis 
points with lots of missing radius values.

"""
stdbarHelp=\
"""To show standard deviation as error-bar
If it is supplied, standard deviation will be shown as an error-bar in the plot.

"""

discardLastsHelp=\
"""Number of smallest clusters to discard from the plotting.
It can be useful to filter out few smallest clusters because these may 
contain small number of frames.

"""

rightmarginHelp=\
"""Margin at right side of the plots.
If legends overflow into the plot area, margin can be increased to fit the legend.

"""

legendcolsHelp=\
"""Number of legend columns
If legend overflow into the plot area, legends can be made of more than 
one column to accommodate all legends.

"""

def main():
    parser, args = parseArguments()
    
    # Input file
    inputFile = None
    if args.inputFile is not None:
        inputFile = args.inputFile.name
        args.inputFile.close()
    else:
        showErrorAndExit(parser, "No Input File!!!\n")

    # Input file
    cluster_file = None
    if args.cluster_file is not None:
        cluster_file = args.cluster_file.name
        args.cluster_file.close()
    else:
        showErrorAndExit(parser, "No cluster id File!!!\n")
        
    # output file
    outputFile = None
    if args.outputFile is not None:
        outputFile = args.outputFile
    else:
        showErrorAndExit(parser, "No Output File!!!\n")
        
        
    # Determine output file-extension type
    output_formats = list(get_output_formats().keys())
    outputFileExtension = os.path.splitext(args.outputFile)[1][1:]
    if outputFileExtension not in output_formats:
        showErrorAndExit(parser, "File extension {0} is not recognized as an"
                         " acceptable extension.\n Use from following: {1}"
                         .format(outputFileExtension, output_formats))

    holeProcessor = holeOutputProcessor.HoleOutputProcessor(inputFile, axis=args.axis, endrad=args.endrad, 
                                                            xmin=args.xmin, xmax=args.xmax, gap=args.gap, 
                                                            begin=args.begin, end=args.end,
                                                            dataOccupancy=args.dataOccupancy)

    holeProcessor.plot_by_cluster(cluster_file, outputFile, csvfile=args.outCsvFile, stdbar=args.stdbar, 
                                  discard_lasts=args.discard_lasts, width=args.width, height=args.height,
                                  ymin=args.ymin, ymax=args.ymax, fontsize=args.fontsize, 
                                  rightmargin=args.rightmargin, legendcols=args.legendcols)



def parseArguments():

    output_formats = list(get_output_formats().keys())
    
    parser = argparse.ArgumentParser(prog='gmx_clusterByFeatures holeclustersplot',
                                     description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', action='store',
                        type=argparse.FileType('r'), metavar='radius.dat',
                        dest='inputFile', required=False, help=inputFileHelp)

    parser.add_argument('-clid', '--clid', action='store',
                        type=argparse.FileType('r'), metavar='clid.xvg',
                        dest='cluster_file', required=False, help=clusterFileHelp)
    
    parser.add_argument('-o', '--output', action='store',
                        type=str, metavar='output.png', 
                        dest='outputFile', help=outputFileHelp)
    
    parser.add_argument('-csv', '--out-csv', action='store',
                        type=str, metavar='output.csv', 
                        dest='outCsvFile', help=outCsvFileHelp)
    
    parser.add_argument('-xmin', '--axis-min', action='store',
                        type=float, dest='xmin', 
                        help=xminHelp)
    
    parser.add_argument('-xmax', '--axis-max', action='store',
                        type=float, dest='xmax', 
                        help=xmaxHelp)
    
    parser.add_argument('-endrad', '--end-radius', action='store',
                        type=float, dest='endrad', default=5, 
                        help=endradHelp)
    
    parser.add_argument('-ax', '--axis', action='store',
                        type=str, metavar='Z', default='Z', 
                        dest='axis', help='Principal axis parallel to the channel or cavity.')
    
    parser.add_argument('-gap', '--gap', action='store',
                        type=float, default=1, metavar=1,
                        dest='gap', help=gapHelp)
    
    parser.add_argument('-b', '--begin', action='store',
                        type=float, default=0, metavar=0,
                        dest='begin', help="First frame in time to read from the input file")

    parser.add_argument('-e', '--end', action='store',
                        type=float, metavar=-1, default=-1, 
                        dest='end', help=endHelp)
    
    parser.add_argument('-do', '--data-occupancy', action='store',
                        type=float, metavar=90, default=90, 
                        dest='dataOccupancy', help=dataOccupancyHelp)
    
    parser.add_argument('-stdbar', '--stddev-bar', action='store_true',
                        default=False, 
                        dest='stdbar', help=stdbarHelp)
    
    parser.add_argument('-dl', '--discard-lasts', action='store',
                        type=int, metavar=0, 
                        dest='discard_lasts', help=discardLastsHelp)
    
    parser.add_argument('-ymin', '--y-axis-min', action='store',
                        type=float, dest='ymin', 
                        help=yminHelp)
    
    parser.add_argument('-ymax', '--y-axis-max', action='store',
                        type=float, dest='ymax', 
                        help=ymaxHelp)
        
    parser.add_argument('-rmargin', '--right-margin', action='store',
                        type=float, metavar=0.15, default=0.15,
                        dest='rightmargin', help=rightmarginHelp)
    
    parser.add_argument('-lcols', '--legend-cols', action='store',
                        type=int, metavar=1, default=1,
                        dest='legendcols', help=legendcolsHelp)
    
    parser.add_argument('-fs', '--font-size', action='store',
                        type=int, metavar=18, default=18,
                        dest='fontsize', help="Font-size of all texts in plot")
    
    parser.add_argument('-wd', '--width', action='store',
                        type=float, metavar=6, default=6,
                        dest='width', help="width of plot in inch")
    
    parser.add_argument('-ht', '--height', action='store',
                        type=float, metavar=6, default=6,
                        dest='height', help="height of plot in inch")
    
    parser.add_argument('-dpi', '--dpi', action='store',
                        type=int, metavar=300, default=300,
                        dest='dpi', help="resolution of plot")

    idx = sys.argv.index("holeclustersplot") + 1
    args = parser.parse_args(args=sys.argv[idx:])

    return parser, args
    


def showErrorAndExit(parser, message):
    parser.print_help()
    print("\n===== ERROR =======")
    print(message)
    print("See Usage Above!!!")
    sys.exit(False)
    
if __name__=="__main__":
    main()

