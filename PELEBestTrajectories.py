# -*- coding: utf-8 -*-


# Imports
from __future__ import unicode_literals
import os
import glob
import argparse as ap


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
    output_path : string
                  output directory where the resulting plot will be saved
    """

    parser = ap.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", required=True, metavar="FILE",
                          type=str, nargs='*', help="path to report files")
    optional.add_argument("-o", "--output", metavar="PATH", type=str,
                          help="output path to save figure", default="Important_trajectories")
    parser._action_groups.append(optional)
    args = parser.parse_args()

    reports = parseReports(args.input, parser)

    output_path = args.output

    return reports, output_path


def Storebesttrajectories(reports):
    """It looks on the report files and finds the best trajectories (the minima).

    RETURNS
    -------
    Below50 : Dictionary of lists
		  The best trajectories to the biggest cutoff.
    Below55 : Dictionary of lists
		  The best trajectories to the intermediate cutoff.
    Below60 : Dictionary of lists
		  The best trajectories to the smallest cutoff.
    """
    Below50={}
    Below55={}
    Below60={}
    for report in reports:
        reportID=str(os.path.basename(report).split('_')[-1].split('.')[0])
        with open(report, 'r') as report_file:
            next(report_file)
	    for i, line in enumerate(report_file):
	        if float(line.split()[4])<=float(-50):
		    if round(float(line.split()[5]))%5<=2:
			Distance=(round(float(line.split()[5]))-(round(float(line.split()[5]))%5))
		    else:
		        Distance=(round(float(line.split()[5]))+5-(round(float(line.split()[5]))%5))
		    if Distance not in Below50:
			Below50[Distance]=[]
		    if reportID not in Below50[Distance]:
			Below50[Distance].append(reportID)
	            if float(line.split()[4])<=float(-55):
		        if Distance not in Below55:
			    Below55[Distance]=[]
		        if reportID not in Below55[Distance]:
			    Below55[Distance].append(reportID)
		        if float(line.split()[4])<=float(-60):
		            if Distance not in Below60:
			        Below60[Distance]=[]
		            if reportID not in Below60[Distance]:
			        Below60[Distance].append(reportID)

    return Below50,Below55,Below60


def Reportbesttrajectories(Below50,Below55,Below60,output_path="Important_trajectories"):
    """It looks on the report files and finds the best trajectories (the minima).

    RETURNS
    -------
    output_path : Output file
		  The report output file with the best trajectories.
    """
    Report=open(output_path,"wt")
    Report.write("--- Below -50 kcal/mol ---\n")
    for key,values in Below50.items():
	Report.write("%s Ang: " %key + ", ".join(values)+"\n")
    Report.write("--- Below -55 kcal/mol ---\n")
    for key,values in Below55.items():
	Report.write("%s Ang: " %key + ", ".join(values)+"\n")
    Report.write("--- Below -60 kcal/mol ---\n")
    for key,values in Below60.items():
	Report.write("%s Ang: " %key + ", ".join(values)+"\n")
    Report.close()


def main():
    """Main function

    It is called when this script is the main program called by the interpreter
    """

    # Parse command-line arguments
    reports, output_path = parseArgs()

    # Store the best trajectories
    Report_data_1,Report_data_2,Report_data_3=Storebesttrajectories(reports)

    # Generate the report file with the best trajectories
    Reportbesttrajectories(Report_data_1,Report_data_2,Report_data_3,output_path)


if __name__ == "__main__":
    """Call the main function"""
    main()
