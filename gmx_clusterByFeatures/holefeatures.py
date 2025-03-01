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

from .holeOutputProcessor import HoleOutputProcessor


description=\
"""DESCRIPTION
===========
Write channel/cavity radius as features for clustering

The output file can be used as features for clustering of channel/cavity.

"""

inputFileHelp=\
"""Name of input file.
Radius file obtained from hole as an output file.

"""

outputFileHelp=\
"""Name of output file.
Output file containing radius as function of time at each axis points.
This file can be used as features file for clustering. This file can be
also used to plot radius vs time with external plotting program.

The file name should end with xvg extension, which is recognized by 
"cluster" command.

"""

pcaHelp=\
"""Number of eigenvectors to be considered for the features.
In place for taking radius as features, this option enable PCA of radii
and the resultant projections on eigenvectors can be used as features.

"""

xminHelp=\
"""Minimum value of axis-point
Minimum value of axis point after which radius value will be considered.

If not supplied, minimum axis value will be extracted from input file.

"""

xmaxHelp=\
"""Maximum value of axis-point
Maximum value of axis point before which radius value will be considered.

If not supplied, maximum axis value will be extracted from input file.

"""

endradHelp=\
"""End/Opening radius
If radius is larger than this value, this value will not considered 
for average calculation and features output. This value might be equal or
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
If an axis-point has radius-data less than this percentage of frame, 
the axis-point will not be considered for average calculation and 
features output.

This is critical for axis-points, which are at the opening of 
channel/cavity. In several frames, radius-value could be missing 
and therefore, dataOccupancy threshold could be used to discard 
those axis points with lots of missing radius values.

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

    # output file
    outputFile = None
    if args.outputFile is not None:
        outputFile = args.outputFile
    else:
        showErrorAndExit(parser, "No Output File!!!\n")
        

    holeProcessor = HoleOutputProcessor(inputFile, axis=args.axis, endrad=args.endrad, xmin=args.xmin, 
                                        xmax=args.xmax, gap=args.gap, begin=args.begin, end=args.end,
                                        dataOccupancy=args.dataOccupancy)
    holeProcessor.write_features(outputFile, pca=args.pca)



def parseArguments():
    
    parser = argparse.ArgumentParser(prog='gmx_clusterByFeatures holefeatures',
                                     description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', action='store',
                        type=argparse.FileType('r'), metavar='radius.dat',
                        dest='inputFile', required=False, help=inputFileHelp)

    
    parser.add_argument('-o', '--output', action='store',
                        type=str, metavar='output.xvg', 
                        dest='outputFile', help=outputFileHelp)
    
    parser.add_argument('-pca', '--pca-pcs', action='store',
                        type=int, dest='pca', default=5,
                        metavar='5', help=pcaHelp)
    
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
                        type=float, default=0,  metavar=0,
                        dest='begin', help="First frame in time to read from the input file")

    parser.add_argument('-e', '--end', action='store',
                        type=float, metavar=-1, default=-1, 
                        dest='end', help=endHelp)
    
    parser.add_argument('-do', '--data-occupancy', action='store',
                        type=float, metavar=90, default=90, 
                        dest='dataOccupancy', help=dataOccupancyHelp)
    

    idx = sys.argv.index("holefeatures") + 1
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

