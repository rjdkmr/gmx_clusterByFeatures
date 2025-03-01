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

from . import featuresplotter

def get_output_formats(fmt=False):
    fig = featuresplotter.plt.figure()
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
Features vs Features plot

This can be used to generate plots for features vs features data.
These type of plots are useful to check quality of clustering.

"gmx_clusterByFeatures cluster" with "-plot" option also produces 
features vs features plot. However, the obtained plot is fixed 
and cannot be changed. Therefore, this sub-command can be used 
to obtained plots for desired features with several different 
options to customize the plot.

"""

inputFileHelp=\
"""Name of input file.
Name of input text file. It should contains two features and their
respective labels in each row. All these values should be separated 
by comma. Each row in file should be in following format:

[feature no. at X-axis], [feature no. at Y-axis],[X-Label],[Y-Label]

For example, following input will result in four plots:

1,2,PC-1,PC-2
2,3,PC-2,PC-3
1,3,PC-1,PC-3
1,4,PC-1,PC-4

"""

featuresfileHelp=\
"""Input features file.
This file should be same as supplied to "gmx_clusterByFeatures cluster" with "-feat" option.

"""

clidfileHelp=\
"""Input file containing cluster-id as a function of time.
The number of frames in this file should be same as in features file.
"""

outputFileHelp=\
"""Name of output file.
Name of the output plot file. The extension will be used to determine
the output format.

Following output formats are available:
{0}

""".format(get_output_formats(fmt=True))


endHelp=\
"""Last frame in time  to read from the input file.
By default ("-e -1"), all frames till the end will be read.

"""

topmarginHelp=\
"""Margin at top side of the plot.
If legends overflow into the plot area, margin can be increased to fit the legend.

"""

legendcolsHelp=\
"""Number of legend columns
If legend overflow the plot area, legends can be made of more than 
one rows by limiting number of columns to accommodate all legends.

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
    featuresfile = None
    if args.featuresfile is not None:
        featuresfile = args.featuresfile.name
        args.featuresfile.close()
    else:
        showErrorAndExit(parser, "No features File!!!\n")
        
    # Input file
    clidfile = None
    if args.clidfile is not None:
        clidfile = args.clidfile.name
        args.clidfile.close()
    else:
        showErrorAndExit(parser, "No cluster-id File!!!\n")

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
        
    featurePlot = featuresplotter.FeaturesPlotter(inputFile, clidfile, featuresfile, 
                                                  nFeatures=None, clusterlogfile=None, 
                                                  begin=args.begin, end=args.end)
    featurePlot.plot_features(outputFile, width=args.width, height=args.height,
                              topmargin=args.topmargin, legendcols=args.legendcols,
                              fontsize=args.fontsize, dpi=args.dpi)
        

def parseArguments():

    output_formats = list(get_output_formats().keys())
    
    parser = argparse.ArgumentParser(prog='gmx_clusterByFeatures featuresplot',
                                     description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', action='store',
                        type=argparse.FileType('r'), metavar='input.txt',
                        dest='inputFile', required=False, help=inputFileHelp)

    parser.add_argument('-feat', '--features', action='store',
                        type=argparse.FileType('r'), metavar='features.xvg',
                        dest='featuresfile', required=False, help=featuresfileHelp)
    
    parser.add_argument('-clid', '--cluster-id', action='store',
                        type=argparse.FileType('r'), metavar='clid.xvg',
                        dest='clidfile', required=False, help=clidfileHelp)
    
    parser.add_argument('-o', '--output', action='store',
                        type=str, metavar='output.png', 
                        dest='outputFile', help=outputFileHelp)
    
    parser.add_argument('-b', '--begin', action='store',
                        type=float, default=0, metavar=0,
                        dest='begin', help="First frame in time to read from the input file")

    parser.add_argument('-e', '--end', action='store',
                        type=float, metavar=-1, default=-1, 
                        dest='end', help=endHelp)

    parser.add_argument('-tmargin', '--top-margin', action='store',
                        type=float, metavar=0.1, default=0.1,
                        dest='topmargin', help=topmarginHelp)
    
    parser.add_argument('-lcols', '--legend-cols', action='store',
                        type=int, metavar=5, default=5,
                        dest='legendcols', help=legendcolsHelp)
    
    parser.add_argument('-fs', '--font-size', action='store',
                        type=int, metavar=18, default=18,
                        dest='fontsize', help="Font-size of all texts in plot")
    
    
    parser.add_argument('-wd', '--width', action='store',
                        type=float, metavar=8, default=8,
                        dest='width', help="width of plot in inch")
    
    parser.add_argument('-ht', '--height', action='store',
                        type=float, metavar=10, default=10,
                        dest='height', help="height of plot in inch")
    
    parser.add_argument('-dpi', '--dpi', action='store',
                        type=int, metavar=300, default=300,
                        dest='dpi', help="resolution of plot")

    idx = sys.argv.index("featuresplot") + 1
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

