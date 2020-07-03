# -*- coding: utf-8 -*-


# Global imports
from __future__ import unicode_literals
import os
import glob
import argparse as ap
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib
matplotlib.use("tkagg")
import seaborn as sns
from math import isnan

# Local imports
from PELEParseReports import *


# Script information
__author__ = ["Marti Municoy", "Sergi Rodà"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = ["Marti Municoy", "Sergi Rodà"]
__email__ = ["marti.municoy@bsc.es","sergi.rodallordes@bsc.es"]

# Possible font of the axis and the title
Dict_of_fonts={"title" : {'family': 'serif',
        'weight': 'bold',
        'size': 16,},
"axis" : {'family': 'serif',
        'weight': 'normal',
        'size': 14,}}


# Functions
def parseArgs():
    """Parse arguments from command-line

    RETURNS
    -------
    reports : string
              list of report files to look for data
    x_data : string
             data to parse and assign to the X axis
    y_data : string
             data to parse and assign to the Y axis
    z_data : string
             data to parse and assign to the colorbar or the 3rd axis
    z_data : string
             data to parse and assign to the colorbar of the 3D scatter plot
    z_max : float
            it sets the maximum range value of the colorbar
    z_min : float
            it sets the minimum range value of the colorbar
    output_path : string
                  output directory where the resulting plot will be saved
    title : string
                  title of the generated figure
    font: list of strings
		  key representing the font style in the font dictionary
    size: integer
          The font size of all the labels in the plot
    scatterplot: boolean
          Perform the 2D scatter plot of two PELE metrics and a 3rd metric as colorbar (hover function)
    densityplot: boolean
          Perform the density plot of one PELE metric
    pointplot: boolean
          Perform the 2D scatter plot of two PELE metrics and a 3rd metric as colorbar
    threeDscatterplot: boolean
          Perform the 3D scatter plot of three PELE metrics and a 4th metric as colorbar (hover function)
    """

    parser = ap.ArgumentParser(description='Script used to generate plots from the metrics \
        of the reports files from a PELE simulation')
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", required=True, metavar="FILE",
                          type=str, nargs='*', help="path to report files")
    optional.add_argument("-X", "--xaxis", metavar="INTEGER [METRIC]",
                          type=str, nargs='*', help="column numbers and " +
                          "metric to plot on the X axis", default=None)
    optional.add_argument("-Y", "--yaxis", metavar="INTEGER [METRIC]",
                          type=str, nargs='*', help="column numbers and " +
                          "metric to plot on the Y axis", default=None)
    optional.add_argument("-Z", "--zaxis", metavar="INTEGER [METRIC]",
                          type=str, nargs='*', help="column numbers and " +
                          "metric to represent in the colorbar", default=None)
    optional.add_argument("-Z2","--z2axis", metavar="INTEGER [METRIC]",
                          type=str, nargs='*', help="column numbers and " +
                          "metric to represent in the colorbar of the 3D scatter plot", default=None)
    optional.add_argument("-o", "--output", metavar="PATH", type=str,
                          help="output path to save figure", default=None)
    optional.add_argument("-r", "--Zmin", metavar="FLOAT", type=float,
                          help="minimum Z value for the colorbar (or Z2 value in the 3D scatter plot)",
                          default=None)
    optional.add_argument("-R", "--Zmax", metavar="FLOAT", type=float,
                          help="maximum Z value for the colorbar (or Z2 value in the 3D scatter plot)",
                          default=None)
    optional.add_argument("-T","--title", metavar="STRING", type=str,
			  help = "title of the figure", default="")
    optional.add_argument("-F","--font", metavar="STRING [STRING]", type=str,nargs='*',
                          help = "a list of the name of the font of the axis and the title", default="")
    optional.add_argument("-S","--size", metavar="INTEGER", type=int,
                          help = "the size of the font for all the plot", default=12)
    optional.add_argument("-SP","--scatterplot"
                          ,help = "Perform the archetypical PELEPlot", action = "store_true")
    optional.add_argument("-DP","--densityplot"
                          ,help = "Perform the densityPlot of the metric specified in X", action = "store_true")
    optional.add_argument("-PP","--pointplot"
                          ,help = "Perform the pointPlot of the metric specified in X and Y", action = "store_true")
    optional.add_argument("-TP","--threeDscatterplot"
                          ,help = "Perform the 3D scatter plot of the metric specified in X, Y, and Z [and Z 2]", action = "store_true")
    parser._action_groups.append(optional)
    args = parser.parse_args()

    reports = parseReports(args.input, parser)

    x_data = args.xaxis
    y_data = args.yaxis
    z_data = args.zaxis
    z2_data = args.z2axis

    if z_data is None:
        if args.Zmin is not None or args.Zmax is not None:
            print("No data to represent in the colorbar (Z axis). " +
                  "Custom ranges are ignored.")
        z_min = None
        z_max = None
    else:
        z_min = args.Zmin
        z_max = args.Zmax

    output_path = args.output
    title = args.title
    font=args.font
    if len(font)==0:
        font=["",""]
    size = args.size
    SP = args.scatterplot
    DP = args.densityplot
    PP = args.pointplot
    TP = args.threeDscatterplot

    return reports, x_data, y_data, z_data, z2_data, z_min, z_max, output_path, title, font, size, SP, DP, PP, TP


def addUnits(metric_name):
    """Add units according to the input metric

    PARAMETERS
    ----------
    metric_name : string
                  name of the metric to plot

    RETURNS
    -------
    label : string
            name of the metric to plot with the units that were added to it
    """

    if "energy" in metric_name.lower():
        label = metric_name + " ($kcal/mol$)"
    elif "energies" in metric_name.lower():
        label = metric_name + " ($kcal/mol$)"
    elif "distance" in metric_name.lower():
        label = metric_name + " ($\AA$)"
    elif "rmsd" in metric_name.lower():
        label = metric_name + " ($\AA$)"
    else:
        label = metric_name
    return label
            
def parseAxisData(axis_data):
    """It sets the columns and label of the data that wants to be plotted

    PARAMETERS
    ----------
    axis_data : string
                axis data to parse

    RETURNS
    -------
    parsed_data : tuple of a list and a string
                  the list specifies the report columns that want to be plotted
                  in the axis and the string sets the name of the axis
    """

    if axis_data is None:
        return ([None, ], None)
    else:
        try:
            rows = [int(axis_data[0]), ]
        except ValueError:
            print("Warning: axis data not recognized: {}".format(axis_data))
            return ([None, ], None)
        if len(axis_data) == 1:
            return (rows, None)
        elif len(axis_data) > 1:
            label_index = 1
            while axis_data[label_index] == "+":
                label_index += 2
                try:
                    rows.append(int(axis_data[2]))
                except (ValueError, IndexError):
                    print("Warning: axis data not recognized: " +
                          "{}".format(axis_data))
                    return ([None, ], None)
                if len(axis_data) == label_index:
                    return (rows, '?')
            label = addUnits(axis_data[label_index])
            return (rows, label)

        print("Warning: axis data not recognized: {}".format(axis_data))
        return ([None, ], None)


def scatterPlot(reports,
                x_rows = [None, ], y_rows = [None, ], z_rows = [None, ],
                x_name = None, y_name = None, z_name = None,
                output_path = None, z_max = None, z_min = None, title = "", font = ["",""], size=12):
    """Represent the scatter plot

    PARAMETERS
    ----------
    reports : string
              list of report files to look for data
    x_rows : list of integers
             integers which specify the report columns to represent in the X
             axis
    y_rows : list of integers
             integers which specify the report columns to represent in the Y
             axis
    z_rows : list of integers
             integers which specify the report columns to represent in the
             colorbar
    x_name : string
             label of the X axis
    y_name : string
             label of the Y axis
    z_name : string
             label of the colorbar
    output_path : string
                  output directory where the resulting plot will be saved
    z_max : float
            it sets the maximum range value of the colorbar
    z_min : float
            it sets the minimum range value of the colorbar
    title: string
	   it sets the title name of the plot
    size: integer
          The font size of all the labels in the plot
    """

    # The different variables are created and the size of the labels is set
    x_values, y_values, z_values, labels, annotations = [], [], [], [], []
    plt.rcParams.update({'font.size': size})

    # Set the rows and their labels to perform the scatter plot if they are not specified
    with open(reports[0], 'r') as report_file:
        line = report_file.readline()
        if None in x_rows:
            x_rows = [7, ]
            x_name = "RMSD ($\AA$)"
        if None in y_rows:
            y_rows = [5, ]
            y_name = "Energy ($kcal/mol$)"
        if x_name is None:
            x_name = str(line.split("    ")[x_rows[0] - 1])
        if y_name is None:
            y_name = str(line.split("    ")[y_rows[0] - 1])
        if (None not in z_rows) and (z_name is None):
            z_name = str(line.split("    ")[z_rows[0] - 1])
            z_name = addUnits(z_name)

    # Get the report files and save the directory where the report is saved and their number
    for report in reports:
        report_directory = os.path.dirname(report)
        report_number = os.path.basename(report).split('_')[-1].split('.')[0]
        
        # Open the report file and save the valeus that will be represented in the 2D scatter plot
        with open(report, 'r') as report_file:
            next(report_file)
            for i, line in enumerate(report_file):
                x_total = 0.
                y_total = 0.
                z_total = 0.

                for x_row in x_rows:
                    x_total += float(line.split()[x_row - 1])

                for y_row in y_rows:
                    y_total += float(line.split()[y_row - 1])

                if None not in z_rows:
                    for z_row in z_rows:
                        z_total += float(line.split()[z_row - 1])

                if isnan(x_total) or isnan(y_total) or isnan(z_total):
                    continue

                x_values.append(x_total)
                y_values.append(y_total)
                z_values.append(z_total)

                epoch = report_directory.split('/')[-1]
                if not epoch.isdigit():
                    epoch = '0'

                annotations.append("Epoch: " + epoch + "\n" +
                                   "Trajectory: " + report_number + "\n" +
                                   "Model: " + str(i + 1))

                labels.append(0)

    if z_max is None:
        z_max = max(z_values)

    if z_min is None:
        z_min = min(z_values)

    if z_min == z_max:
        cmap = plt.cm.autumn
    else:
        cmap = plt.cm.plasma

    norm = plt.Normalize(z_min, z_max)

    fig, ax = plt.subplots()

    if output_path is not None:
        s = 20
    else:
        s = None


    sc = plt.scatter(x_values, y_values, c=z_values, cmap=cmap, s=s,
                        norm=norm)

    ax.margins(0.05)
    ax.set_facecolor('white')
    if font[0]!="":
        plt.ylabel(y_name,Dict_of_fonts[font[0]])
        plt.xlabel(x_name,Dict_of_fonts[font[0]])
    else:
        plt.ylabel(y_name)
        plt.xlabel(x_name)
    if font[1]!="":
        plt.title(title,Dict_of_fonts[font[1]])
    else:
        plt.title(title)

    annot = ax.annotate("", xy=(0, 0), xytext=(20, 20),
                        textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    # Activate the colorbar only if the Z axis contains data to plot
    if None not in z_rows:
        cbar = plt.colorbar(sc, drawedges=False)
        if font[0]!="":
            cbar.set_label(z_name,Dict_of_fonts[font[0]])
        else:
            cbar.set_label(z_name)

    def update_annot(ind):
        """Update the information box of the selected point"""
        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        annot.set_text(annotations[int(ind["ind"][0])])
        annot.get_bbox_patch().set_facecolor(cmap(norm(
            z_values[ind["ind"][0]])))

    def hover(event):
        """Action to perform when hovering the mouse on a point"""
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    # Respond to mouse motion
    fig.canvas.mpl_connect("motion_notify_event", hover)

    # Save or display the plot depending on whether an output path was set or not
    if output_path is not None:
        plt.savefig(output_path)
    else:
        plt.show()

def densityPlot(reports, x_rows = [None, ], x_name = None, title = "", size = 12):
    """Represent the density plot

    PARAMETERS
    ----------
    reports : string
              list of report files to look for data
    x_rows : list of integers
             integers which specify the report columns to represent in the X
             axis
    x_name : string
             label of the X axis
    title: string
       it sets the title name of the plot
    size: integer
          The font size of all the labels in the plot
    """
    
    x_values = []
    plt.rcParams.update({'font.size': size})

    for report in reports:
        report_file = open(report,'r')
        i=0
        for line in report_file:
            if i!=0:
                x_values.append(float(line.split()[x_rows[0] - 1]))
            else:
                pass
            i+=1

    sns.distplot(x_values)
    plt.title(title)
    plt.xlabel(x_name)
    if x_name !=None:
        plt.ylabel("Density (1/{})".format(x_name.split("(")[1][0:-1]))
    plt.show()

    
def pointPlot(reports, x_rows = [None, ], x_name = None, y_rows = [None, ], y_name = None, title = "", size = 12):
    """Represent the scatter point plot

    PARAMETERS
    ----------
    reports : string
              list of report files to look for data
    x_rows : list of integers
             integers which specify the report columns to represent in the X
             axis
    x_name : string
             label of the X axis
    y_rows : list of integers
             integers which specify the report columns to represent in the Y
             axis
    y_name : string
             label of the Y axis
    title: string
       it sets the title name of the plot
    size: integer
          The font size of all the labels in the plot
    """    

    x_values,y_values = [],[]
    plt.rcParams.update({'font.size': size})

    for report in reports:
        report_file = open(report,'r')
        i=0
        for line in report_file:
            if i!=0:
                x_values.append(float(line.split()[x_rows[0] - 1]))
                y_values.append(float(line.split()[y_rows[0] - 1]))
            else:
                pass
            i+=1
            
    plt.scatter(x_values,y_values)
    plt.title(title)
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    plt.show()

def ThreeDPlot(reports, x_rows = [None, ], x_name = None, y_rows = [None, ], y_name = None,
    z_rows = [None, ], z_name = None, z2_rows = [None, ], z2_name = None,
    z_max = None, z_min = None, output_path = None, title = "", font = ["",""], size = 12):
    """Represent the 3D scatter plot of some PELE metrics

    PARAMETERS
    ----------
    reports : string
              list of report files to look for data
    x_rows : list of integers
             integers which specify the report columns to represent in the X
             axis
    x_name : string
             label of the X axis
    y_rows : list of integers
             integers which specify the report columns to represent in the Y
             axis
    y_name : string
             label of the Y axis
    z_rows : list of integers
             integers which specify the report columns to represent in the Z
             axis
    z_name : string
             label of the Z axis
    z2_rows : list of integers
             integers which specify the report columns to represent in the colorbar
    z2_name : string
             label of the colorbar
    output_path : string
                  output directory where the resulting plot will be saved
    z_max : float
            it sets the maximum range value of the colorbar
    z_min : float
            it sets the minimum range value of the colorbar
    title: string
       it sets the title name of the plot
    size: integer
          The font size of all the labels in the plot
    """

    # The different variables are created and the size of the labels is set
    x_values, y_values, z_values, z2_values, labels, annotations = [], [], [], [], [], []
    plt.rcParams.update({'font.size': size})

    # Set the rows and their labels to perform the scatter plot if they are not specified
    with open(reports[0], 'r') as report_file:
        line = report_file.readline()
        if None in x_rows:
            x_rows = [7, ]
            x_name = "RMSD ($\AA$)"
        if None in y_rows:
            y_rows = [5, ]
            y_name = "Energy ($kcal/mol$)"
        if None in z_rows:
            z_rows = [6, ]
            z_name = "SASA of the ligand"
        if x_name is None:
            x_name = str(line.split("    ")[x_rows[0] - 1])
        if y_name is None:
            y_name = str(line.split("    ")[y_rows[0] - 1])
        if z_name is None:
            z_name = str(line.split("    ")[z_rows[0] - 1])
        if (None not in z2_rows) and (z2_name is None):
            z2_name = str(line.split("    ")[z2_rows[0] - 1])
            z2_name = addUnits(z2_name)

    # Get the report files and save the directory where the report is saved and their number
    for report in reports:
        report_directory = os.path.dirname(report)
        report_number = os.path.basename(report).split('_')[-1].split('.')[0]

        # Open the report file and save the valeus that will be represented in the 3D scatter plot
        with open(report, 'r') as report_file:
            next(report_file)
            for i, line in enumerate(report_file):
                x_total = 0.
                y_total = 0.
                z_total = 0.
                z2_total = 0.

                for x_row in x_rows:
                    x_total += float(line.split()[x_row - 1])

                for y_row in y_rows:
                    y_total += float(line.split()[y_row - 1])

                for z_row in z_rows:
                    z_total += float(line.split()[z_row - 1])

                if None not in z2_rows:
                    for z2_row in z2_rows:
                        z2_total += float(line.split()[z2_row - 1])

                if isnan(x_total) or isnan(y_total) or isnan(z_total) or isnan(z2_total):
                    continue

                x_values.append(x_total)
                y_values.append(y_total)
                z_values.append(z_total)
                z2_values.append(z2_total)

                epoch = report_directory.split('/')[-1]
                if not epoch.isdigit():
                    epoch = '0'

                annotations.append("Epoch: " + epoch + "\n" +
                                   "Trajectory: " + report_number + "\n" +
                                   "Model: " + str(i + 1))

                labels.append(0)

    if z_max is None:
        z_max = max(z2_values)

    if z_min is None:
        z_min = min(z2_values)

    if z_min == z_max:
        cmap = plt.cm.autumn
    else:
        cmap = plt.cm.plasma

    norm = plt.Normalize(z_min, z_max)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d') 
    sc = ax.scatter(x_values, y_values, z_values, c=z2_values, cmap=cmap,
                        norm=norm)

    ax.margins(0.05)
    ax.set_facecolor('white')
    if font[0]!="":
        ax.set_ylabel(y_name,Dict_of_fonts[font[0]])
        ax.set_xlabel(x_name,Dict_of_fonts[font[0]])
        ax.set_zlabel(z_name,Dict_of_fonts[font[0]])
    else:
        ax.set_ylabel(y_name)
        ax.set_xlabel(x_name)
        ax.set_zlabel(z_name)
    if font[1]!="":
        ax.set_title(title,Dict_of_fonts[font[1]])
    else:
        ax.set_title(title)

    annot = ax.annotate("", xy=(0, 0), xytext=(20, 20),
                        textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    # Activate the colorbar only if the Z2 axis contains data to plot
    if None not in z2_rows:
        cbar = plt.colorbar(sc, drawedges=False)
        if font[0]!="":
            cbar.set_label(z2_name,Dict_of_fonts[font[0]])
        else:
            cbar.set_label(z2_name)

    def update_annot(ind):
        """Update the information box of the selected point"""
        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        annot.set_text(annotations[int(ind["ind"][0])])
        annot.get_bbox_patch().set_facecolor(cmap(norm(
            z2_values[ind["ind"][0]])))

    def hover(event):
        """Action to perform when hovering the mouse on a point"""
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    # Respond to mouse motion
    fig.canvas.mpl_connect("motion_notify_event", hover)

    # Save or display the plot depending on whether an output path was set or not
    if output_path is not None:
        plt.savefig(output_path)
    else:
        plt.show()

def main():
    """Main function

    It is called when this script is the main program called by the interpreter
    """

    # Parse command-line arguments
    reports, x_data, y_data, z_data, z2_data, z_min, z_max, output_path, title, font, size, SP, DP, PP, TP = parseArgs()

    # Parse axis data to label it properly
    x_rows, x_name = parseAxisData(x_data)
    y_rows, y_name = parseAxisData(y_data)
    z_rows, z_name = parseAxisData(z_data)
    z2_rows, z2_name = parseAxisData(z2_data)

    # Generate the plot
    if SP:
        scatterPlot(reports,
                    x_rows=x_rows, y_rows=y_rows, z_rows=z_rows,
                    x_name=x_name, y_name=y_name, z_name=z_name,
                    z_min=z_min, z_max=z_max,
                    output_path=output_path,title=title,font=font,size=size)
    if DP:
        densityPlot(reports,
                    x_rows=x_rows,x_name=x_name,title=title,size=size)
    if PP:
        pointPlot(reports,
                  x_rows=x_rows, y_rows=y_rows,x_name=x_name,y_name=y_name,title=title,size=size)
    if TP:
        ThreeDPlot(reports,
                    x_rows=x_rows, y_rows=y_rows, z_rows=z_rows, z2_rows=z2_rows,
                    x_name=x_name, y_name=y_name, z_name=z_name, z2_name=z2_name,
                    z_min=z_min, z_max=z_max,
                    output_path=output_path,title=title,font=font,size=size)


if __name__ == "__main__":
    """Call the main function"""
    main()
