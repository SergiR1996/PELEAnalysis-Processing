# -*- coding: utf-8 -*-


# Imports
from __future__ import unicode_literals
import os
import glob
import argparse as ap
import numpy as n


# Script information
__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"

# Functions
def parseReports(reports_to_parse, parser):
    """It identifies the reports to add to the plot

    PARAMETERS
    ----------
    reports_to_parse : list of strings
                       all the report files that want to be added to the plot
    parser : ArgumentParser object
             contains information about the command line arguments

    RETURNS
    -------
    parsed_data : tuple of a list and a string
                  the list specifies the report columns that want to be plotted
                  in the axis and the string sets the name of the axis
    """

    reports = []

    for reports_list in reports_to_parse:
        trajectories_found = glob.glob(reports_list)
        if len(trajectories_found) == 0:
            print("Warning: path to report file \'" +
                  "{}".format(reports_list) + "\' not found.")
        for report in glob.glob(reports_list):
            reports.append(report)

    if len(reports) == 0:
        print("Error: list of report files is empty.")
        parser.print_help()
        exit(1)

    return reports

def parseArgs():
    """Parse arguments from command-line

    RETURNS
    -------
    reports : string
              list of report files to look for data
    energy : float
              Cutoff of the binding energy
    sasa :  float
              Cutoff the SASA parameter
    output_path : string
                  output directory where the resulting plot will be saved
    """

    parser = ap.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", required=True, metavar="FILE",
                          type=str, nargs='*', help="path to report files")
    optional.add_argument("-o", "--output", metavar="PATH", type=str,
                          help="output path to save figure", default="PELE_results")
    parser._action_groups.append(optional)
    args = parser.parse_args()

    reports = parseReports(args.input, parser)

    output_path = args.output

    return reports, output_path

def Storeresults(reports):
    """Take the PELE simulation report files and returns the results stored in a dict

    RETURNS
    -------
    Results: dictionary of lists
             dictionary containing the mean of the different quantitative parameters
    """
    Results = {}
    for report in reports:
        reportID = str(os.path.basename(report).split('_')[-1].split('.')[0])
        with open(report,'r') as report_file:
            if '0' not in Results:
                line = report_file.readline()
                Results['0'] = line.split()[3:(len(line.split()))]
                Means = [[] for x in range(3,(len(line.split())))]
            else:
                next(report_file)
            for i, line in enumerate(report_file):
                for element,i in zip(line.split(),range(len(line.split()))):
                    if i>=3:
                        Means[i-3].append(float(element))
            Results[reportID]=[str(n.mean(Means[i])) for i in range(len(Means))]
    return Results

def Outputresults(Results,output_path):
    """Take the PELE simulation results and report the means with a csv table

    RETURNS
    -------
    Results: dictionary of lists
             dictionary containing the mean of the different quantitative parameters
    """
    Output_file = open(output_path,"wt")
    for key,values in sorted(Results.items()):
        Output_file.write(str(key)+','+",".join(values)+"\n")
    Output_file.close()

                

def main():
    """Main function

    It is called when this script is the main program called by the interpreter
    """

    # Parse command-line arguments
    reports, output_path = parseArgs()

    # Store the results in a dictionary
    Results = Storeresults(reports)

    # Save the stored results in an output file
    Outputresults(Results,output_path)


if __name__ == "__main__":
    """Call the main function"""
    main()

