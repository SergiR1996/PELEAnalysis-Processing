# -*- coding: utf-8 -*-


# Global imports
from __future__ import unicode_literals
import os
import glob
import argparse as ap

# Local imports
from PELEParseReports import *

# Script information
__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"


# Functions
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
                  output directory where the resulting best trajectories will be saved
    """

    parser = ap.ArgumentParser(description='Script used to find the most interesting trajectories \
        acoording to a numerical metric against the binding energy or SASA of report files from a PELE simulation')
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", required=True, metavar="FILE",
                          type=str, nargs='*', help="path to report files")
    optional.add_argument("-o", "--output", metavar="PATH", type=str,
                          help="output path to save figure", default="Important_trajectories")
    optional.add_argument("-E", "--energy", metavar="METRIC", type=float,
                          help="Cutoff of the binding energy", default=-50.0)
    optional.add_argument("-S", "--sasa", metavar="METRIC", type=float,
                          help="Cutoff of the SASA parameter", default=0.3)
    optional.add_argument("-M","--metric",metavar="NUMBER",type=int,
                          help="Number of the column where the metric resides", default=5)
    parser._action_groups.append(optional)
    args = parser.parse_args()

    reports = parseReports(args.input, parser)

    output_path = args.output
    energy,sasa,metric = args.energy,args.sasa,args.metric

    return reports, energy, sasa, output_path,metric


def Storebesttrajectories(reports,metric,energy=-50.0,sasa=0.3):
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

    Below50,Below55,Below60,Sasa03={},{},{},{}

    for report in reports:
        reportID=str(os.path.basename(report).split('_')[-1].split('.')[0])
        with open(report, 'r') as report_file:
            next(report_file)
            for i, line in enumerate(report_file):
                if float(line.split()[4])<=energy:
                    if round(float(line.split()[metric]))%5<=2:
                        Distance=(round(float(line.split()[metric]))-(round(float(line.split()[metric]))%5))
                    else:
                        Distance=(round(float(line.split()[metric]))+5-(round(float(line.split()[metric]))%5))
                    if Distance not in Below50:
                        Below50[Distance]=[]
                    if reportID not in Below50[Distance]:
                        Below50[Distance].append(reportID)
                    if float(line.split()[4])<=(energy-5.0):
                        if Distance not in Below55:
                            Below55[Distance]=[]
                        if reportID not in Below55[Distance]:
                            Below55[Distance].append(reportID)
                        if float(line.split()[4])<=(energy-10.0):
                            if Distance not in Below60:
                                Below60[Distance]=[]
                            if reportID not in Below60[Distance]:
                                Below60[Distance].append(reportID)
                if float(line.split()[len(line.split())-1])<=sasa:
                    if round(float(line.split()[metric]))%5<=2:
                        DistSASA=(round(float(line.split()[metric]))-(round(float(line.split()[metric]))%5))
                    else:
                        DistSASA=(round(float(line.split()[metric]))+5-(round(float(line.split()[metric]))%5))
                    if DistSASA not in Sasa03:
                        Sasa03[DistSASA]=[]
                    if reportID not in Sasa03[DistSASA]:
                        Sasa03[DistSASA].append(reportID)

    return Below50,Below55,Below60,Sasa03


def Reportbesttrajectories(Below50,Below55,Below60,Sasa03,energy=-50.0,sasa=0.3,output_path="Important_trajectories"):
    """It looks on the report files and finds the best trajectories (the minima).

    RETURNS
    -------
    output_path : Output file
		  The report output file with the best trajectories.
    """

    Report=open(output_path,"wt")
    Report.write("--- Below %s kcal/mol ---\n" %energy)

    for key,values in Below50.items():
        Report.write("%s Ang: " %key + ", ".join(values)+"\n")
    Report.write("--- Below %s kcal/mol ---\n" %(energy-5.0))
    for key,values in Below55.items():
        Report.write("%s Ang: " %key + ", ".join(values)+"\n")
    Report.write("--- Below %s kcal/mol ---\n" %(energy-10.0))
    for key,values in Below60.items():
        Report.write("%s Ang: " %key + ", ".join(values)+"\n")
    Report.write("--- SASA below %s ---\n" %sasa)
    for key,values in Sasa03.items():
        Report.write("%s Ang: " %key + ", ".join(values)+"\n")
        
    Report.close()


def main():
    """Main function

    It is called when this script is the main program called by the interpreter
    """

    # Parse command-line arguments
    reports, energy, sasa, output_path, metric = parseArgs()

    # Store the best trajectories
    Report_data_1,Report_data_2,Report_data_3,Report_data_4=Storebesttrajectories(reports,metric,energy,sasa)

    # Generate the report file with the best trajectories
    Reportbesttrajectories(Report_data_1,Report_data_2,Report_data_3,Report_data_4,energy,sasa,output_path)


if __name__ == "__main__":
    """Call the main function"""
    main()
