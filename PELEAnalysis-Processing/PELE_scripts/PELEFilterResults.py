# -*- coding: utf-8 -*-


# Global imports
from __future__ import unicode_literals
import os
import glob
import argparse as ap
import pandas as pd
import numpy as n

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
    output_path : string
                  output directory where the csv file will be saved
    column : list of indices
                  list of column indices
    """

    parser = ap.ArgumentParser(description='Script that returns a csv file with the mean of the numerical \
        metrics of the reports file from a PELE simulation')
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", required=True, metavar="FILE",
                          type=str, nargs='*', help="path to report files")
    optional.add_argument("-o", "--output", metavar="PATH", type=str,
                          help="output path to save figure", default="PELE_results")
    optional.add_argument("-C","--column", metavar="LIST",type=str,
                          nargs='*',help="index of the column where the filtering will be applied")
    parser._action_groups.append(optional)
    args = parser.parse_args()

    reports = parseReports(args.input, parser)

    output_path, column = args.output, args.column

    if column is not None:
        column = [int(i)-1 for i in column]

    return reports, output_path, column


def Storeresults(reports,CE, output_path):
    """
    Take the PELE simulation report files and returns the filtered report files

    OUTPUT
    ------
    filtered report files according to some specified filters in some columns
    """

    if not os.path.exists(os.path.join(os.getcwd(), "Filtered_reports")):
        os.mkdir(os.path.join(os.getcwd(), "Filtered_reports"))
    Report_path= os.path.join(os.getcwd(), "Filtered_reports")

    for i,report in enumerate(reports):
        df = pd.read_csv(report,sep="    ")

        df_aux = df[(df[df.columns[CE[0]]].between(0,100)) & (df[df.columns[CE[1]]].between(-30,-10))]
        df_aux.to_csv(os.path.join(Report_path,output_path+"_"+str(i+1)+".out"))

        report = open(os.path.join(Report_path,output_path+"_"+str(i+1)+".out"),"r")
        report_def = open(os.path.join(Report_path,output_path+"f_"+str(i+1)+".out"),"w")

        for line in report:
            report_def.write("    ".join(line.split(",")[1:]))

        os.system("rm {}".format(os.path.join(Report_path,output_path+"_"+str(i+1)+".out")))
        os.system("mv {} {}".format(os.path.join(Report_path,output_path+"f_"+str(i+1)+".out"),
            os.path.join(Report_path,output_path+"_"+str(i+1)+".out")))


def main():
    """Main function

    It is called when this script is the main program called by the interpreter
    """

    # Parse command-line arguments
    reports, output_path, catalytic_event  = parseArgs()

    # Store the filtered report files
    Storeresults(reports,catalytic_event, output_path)


if __name__ == "__main__":
    """Call the main function"""
    main()