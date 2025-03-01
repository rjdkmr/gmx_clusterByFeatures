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
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar


def get_output_formats(fmt=False):
    fig = plt.figure()
    output_formats_dict = fig.canvas.get_supported_filetypes()
    del fig
    
    if fmt:
        result = ''
        for key in output_formats_dict:
            result += key + ' : ' + output_formats_dict[key] + '\n'
        return result
    else:
        return output_formats_dict
    
def print_colormaps():
    cmaps = sorted(plt.colormaps())
    result = ''
    i = 0
    n = 1
    while (i<len(cmaps)):
        if cmaps[i][-2:] == '_r':
            i += 1
            continue
        if n%5 == 0:
            result += '\n'
            n = 0
        result += '{0:16} '.format(cmaps[i])
        n += 1
        i += 1
    
    return result


description="""DESCRIPTION
===========
Matrix plot

This can be used to generate plots for outputs generated from distmat.

"""

inputFileHelp=\
"""Name of input file.
Output matrix file obtained from distmat.

"""

outputFileHelp=\
"""Name of output file.
Name of the output matrix-plot file. The extension will be used to determine
the output format.

Following output formats are available:
{0}

""".format(get_output_formats(fmt=True))

minValueHelp=\
"""Minimum value to begin color-mapping
If not provided, minimum value of whole matrix will be considered.

"""

maxValueHelp=\
"""Maximum value to end color-mapping
If not provided, maximum value of whole matrix will be considered.

"""

aspectHelp=\
"""Controls the aspect ratio of the axes.
The aspect is of particular relevance for images since it may distort 
the image, i.e. pixel will not be square.

equal : Ensures an aspect ratio of 1. Pixels will be square.
auto  : The axes is kept fixed and the aspect is adjusted so
        that the data fit in the axes. In general, this will 
        result in non-square pixels.

"""

colormapHelp=\
"""Name of colormap.
Name of colormap by which matrix image will be colored.
To preview the available colormaps, please visit following link:
https://matplotlib.org/tutorials/colors/colormaps.html#classes-of-colormaps

Following colormaps are available:
{0}

All the above colormaps are also available in reverse with same name suffixed by 
"_r". For example, reverse of binary colormap is binary_r, reverse of gist_earth
is gist_earth_r etc.

""".format(print_colormaps())


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
        
    # Determine output file-extension type
    output_formats = list(get_output_formats().keys())
    outputFileExtension = os.path.splitext(args.outputFile)[1][1:]
    if outputFileExtension not in output_formats:
        showErrorAndExit(parser, "File extension {0} is not recognized as an"
                         " acceptable extension.\n Use from following: {1}"
                         .format(outputFileExtension, output_formats))


    mpl.rcParams['font.size'] = args.fontsize
    fig = plt.figure(figsize=(args.width, args.height))
    ax = fig.add_subplot(111)
    
    # Read file and plot
    data = read_mat_file(inputFile)
    image = ax.imshow(data, origin='lower', cmap=args.colormap, vmin=args.vmin, vmax=args.vmax, aspect=args.aspect)
    
    
    # Customize X-ticks and X-ticklabels
    # print(data.shape)
    xticks = ax.get_xticks()
    diff = int(xticks[1] - xticks[0])
    xticks = list(range(0, int(data.shape[1])+1, diff))
    xticklabels = list(range(args.xstart, args.xstart+int(data.shape[1])+1, diff))
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)

    # Customize Y-ticks and Y-ticklabels
    yticks = ax.get_yticks()
    diff = int(yticks[1] - yticks[0])
    yticks = list(range(0, int(data.shape[0])+1, diff))
    yticklabels = list(range(args.ystart, args.ystart+int(data.shape[0])+1, diff))
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    
    # X-label and ylabel
    ax.set_xlabel(args.xlabel)
    ax.set_ylabel(args.ylabel)
    
    # Set colorbar
    if args.orientation == 'horizontal':
        divider = make_axes_locatable(ax)
        # add an axes to the right of the main axes.
        cax = divider.append_axes("top", size="5%", pad="2%")
        cbar = plt.colorbar(image, cax=cax, orientation=args.orientation)
        cax.xaxis.set_ticks_position("top")
        cax.xaxis.set_label_position("top")
        cbar.ax.set_xlabel(args.cblabel)
    else:
        divider = make_axes_locatable(ax)
        # add an axes to the right of the main axes.
        cax = divider.append_axes("right", size="5%", pad="2%")
        cbar = plt.colorbar(image, cax=cax, orientation=args.orientation)
        cbar.ax.set_ylabel(args.cblabel)
    
    fig.tight_layout()
    plt.savefig(outputFile, dpi=args.dpi)


def read_mat_file(filename):
    fin = open(filename, 'r')

    data = []
    for line in fin:
        line = line.rstrip().lstrip()
        if not line:
            continue

        temp = re.split('\s+', line)
        data.append(list(map(float, temp)))

    data = np.asarray(data)

    return data.T


def parseArguments():

    output_formats = list(get_output_formats().keys())
    
    parser = argparse.ArgumentParser(prog='gmx_clusterByFeatures matplot',
                                     description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', action='store',
                        type=argparse.FileType('r'), metavar='distmat.dat',
                        dest='inputFile', required=False, help=inputFileHelp)

    
    parser.add_argument('-o', '--output', action='store',
                        type=str, metavar='output.png', 
                        dest='outputFile', help=outputFileHelp)
    
    parser.add_argument('-xs', '--x-start', action='store',
                        type=int, metavar=1, default=1,
                        dest='xstart', help="First residue number along X-axis")
    
    parser.add_argument('-ys', '--y-start', action='store',
                        type=int, metavar=1, default=1,
                        dest='ystart', help="First residue number along Y-axis")
    
    parser.add_argument('-xl', '--x-label', action='store',
                        type=str, metavar='Residue', default='Residue',
                        dest='xlabel', help="X-axis label")
    
    parser.add_argument('-yl', '--y-label', action='store',
                        type=str, metavar='Residue', default='Residue',
                        dest='ylabel', help="Y-axis label")
    
    parser.add_argument('-cbl', '--colorbar-label', action='store',
                        type=str, metavar='(nm)', default='(nm)',
                        dest='cblabel', help="label for color bar")
    
    parser.add_argument('-a', '--image-aspect', action='store',
                        type=str, metavar='auto', default='auto',
                        choices=['equal', 'auto'],
                        dest='aspect', help=aspectHelp)

    parser.add_argument('-cmap', '--colormap', action='store',
                        type=str, metavar='binary', default='binary',
                        choices=plt.colormaps(),
                        dest='colormap', help=colormapHelp)
    
    parser.add_argument('-vmin', '--min-value', action='store',
                        type=float, dest='vmin', 
                        help=minValueHelp)
    
    parser.add_argument('-vmax', '--max-value', action='store',
                        type=float, dest='vmax', 
                        help=maxValueHelp)

    parser.add_argument('-fs', '--font-size', action='store',
                        type=int, metavar=14, default=14,
                        dest='fontsize', help="Font-size of all texts in plot")

    parser.add_argument('-cbor', '--colorbar-orientation', action='store',
                        type=str, metavar='vertical', default='vertical',
                        choices=['vertical', 'horizontal'],
                        dest='orientation', help="Orientation of color bar")
    
    parser.add_argument('-wd', '--width', action='store',
                        type=float, metavar=8, default=8,
                        dest='width', help="width of plot in inch")
    
    parser.add_argument('-ht', '--height', action='store',
                        type=float, metavar=8, default=8,
                        dest='height', help="height of plot in inch")
    
    parser.add_argument('-dpi', '--dpi', action='store',
                        type=int, metavar=300, default=300,
                        dest='dpi', help="resolution of plot")

    idx = sys.argv.index("matplot") + 1
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
