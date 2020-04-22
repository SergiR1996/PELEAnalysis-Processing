# -*- coding: utf-8 -*-

# Global imports
from __future__ import unicode_literals
import os
import glob
import math
import argparse as ap
import pandas as pd
import numpy as n

# Local imports
from PELEParseReports import *

# Script information
__author__ = "Sergi Rodà"
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Sergi Rodà"
__email__ = "sergi.rodallordes@bsc.es"

# Functions
def parseArgs():
    """
    Parse arguments from command-line

    RETURNS
    -------
    out_directory : string
              path of the PELE simulation output directory
    """

    parser = ap.ArgumentParser(description='Script that returns a report file with the desired metric from a PELE simulation')
    required = parser.add_argument_group('required arguments')
    optional = parser._action_groups.pop()
    required.add_argument("-i", "--input", required=True, metavar="DIRECTORY",
                          type=str, help="path to report files")
    required.add_argument("-R", "--residue", metavar="INTEGER", type=str,nargs='*',
                          help="residue number of the residues involved in the metric")
    required.add_argument("-M", "--metric", metavar="STRING", type=str,nargs='*',
                          help="desired metric to calculate with the atom names (indicate with _)")
    required.add_argument("-RN", "--res_name", metavar="STRING", type=str,nargs='*',
                          help="name of the two residues involved in the distance metric to avoid confusions")
    optional.add_argument("-CN","--column_name", metavar="STRING",type=str,
                          help="column name of the new metric", default="New_metric")
    optional.add_argument("-RF","--report_format", metavar="STRING",type=str,
                          help="Format of the report file", default="*out")
    optional.add_argument("-TF","--trajectory_format", metavar="STRING",type=str,
                          help="Format of the trajectory file", default="*pdb")
    parser._action_groups.append(optional)
    args = parser.parse_args()

    out_directory, residue, metric, res_name, column_name, report_format, trajectory_format = args.input, args.residue, args.metric, args.res_name, args.column_name, args.report_format, args.trajectory_format

    return out_directory, residue, metric, res_name, column_name, report_format, trajectory_format

def AddMetric(out_directory, residue, metric, res_name, column_name, report_format, trajectory_format):
    """
    Take the PELE simulation trajectory files and returns the report files with the desired metric

    OUTPUT
    ------
    Report files with the desired metric added.
    """

    def ComputeDistance(atom1, atom2):

        r = [atom2[0] - atom1[0], atom2[1] - atom1[1], atom2[2] - atom1[2]]
        return math.sqrt(r[0] ** 2 + r[1] ** 2 + r[2]**2)

    trajectory_list = glob.glob(os.path.join(out_directory,trajectory_format))
    report_list = glob.glob(os.path.join(out_directory,report_format))
    report_list.sort()
    trajectory_list.sort()
    metric_list = []

    for trajectory in trajectory_list:
        with open(trajectory, 'r') as traj_file:
            m_list,i = [],0
            for line in traj_file:

                if (line[0:4] == "ATOM" or line[0:6] == "HETATM") and (line[22:26].strip() == residue[0] and line[12:16].replace(" ","_") == metric[0] and line[17:20] == res_name[0]) or (line[22:26].strip() == residue[1] and line[12:16].replace(" ","_") == metric[1] and line[17:20] == res_name[1]):

                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())

                    if i==0:
                        first_elem = [x,y,z]
                        i+=1
                    else:
                        second_elem = [x,y,z]
                        m_list.append(ComputeDistance(first_elem,second_elem))
                        i=0

            metric_list.append(m_list)

    for report in report_list:
        #print(report)

        if 'j' not in locals():
            j=0
        else:
            j+=1

        with open(report, 'r') as report_file:
            if report.endswith(".out"):
                out_report = open(report.split(".out")[0]+"_metric.out",'w')
            else:
                out_report = open(report+"_metric.out",'w')
            for i,line in enumerate(report_file):
                if i==0:
                    out_report.write(line.strip("\n")+'    '+column_name+"\n")
                else:
                    out_report.write(line.strip("\n")+'    '+str(metric_list[j][i-1])+"\n")
    
    print("%s report files have the desired metric added"%j)

def main():
    """
    Main function

    It is called when this script is the main program called by the interpreter
    """

    # Parse command-line arguments
    out_directory, residue, metric, res_name, column_name, report_format, trajectory_format  = parseArgs()

    # Add the desired metric to the report file
    AddMetric(out_directory, residue, metric, res_name, column_name, report_format, trajectory_format)


if __name__ == "__main__":
    """Call the main function"""
    main()


