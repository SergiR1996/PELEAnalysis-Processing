# -*- coding: utf-8 -*-

# Global imports
from __future__ import unicode_literals
import os
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
    reports : string
              list of report files to look for data
    output_path : string
                  output directory where the csv file will be saved
    metric : index (or list of indices)
                  column index (or list of column indices)
    threshold: integer
                  value(s) of the metric(s) that separates the states
    """

    parser = ap.ArgumentParser(description='Script that outputs the difference in binding energy of two states \
        in the PELE simulation based on a threshold in a metric. It also returns the number of accepted PELE steps \
        from each differentiated state.')
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", required=True, metavar="FILE",
                          type=str, nargs='*', help="path to report files")
    optional.add_argument("-o", "--output", metavar="PATH", type=str,
                          help="output path to save figure", default="Estimated_Diff_GibbsE")
    optional.add_argument("-M","--metric", metavar="LIST",type=str,
                          nargs='*',help="index (or indices) of the column(s) that will be used as the distinguishing metric")
    optional.add_argument("-T","--threshold", metavar="LIST",type=str,
                          nargs='*',help="value(s) of the metric(s) that separates the states")
    parser._action_groups.append(optional)
    args = parser.parse_args()

    reports = parseReports(args.input, parser)

    output_path, metric, threshold = args.output, args.metric, args.threshold

    if metric is not None:
        metric = [int(i)-1 for i in metric]

    return reports, output_path, metric, threshold

def CalculateFreeEnergy(reports,output_path,metric,threshold):
    """
    Take the PELE simulation report files and returns the estimated difference 
    in the free energy of two differentiated states according to a/some metric/s.

    OUTPUT
    ------
    output_path.csv File with the difference in the binding energies and the 
    number of accepted PELE steps of each state.
    """

    Results, G1, G2 = {},[],[]

    for i,report in enumerate(reports):
        rep = pd.read_csv(report,sep="    ")
        for i_row in range(rep.shape[0]):
          if rep.loc[i_row][metric[0]]<=float(threshold[0]):
            G1.append(rep.loc[i_row][4])
          else:
            G2.append(rep.loc[i_row][4])
      
    mean_G1 = n.mean(G1,axis=0); std_G1 = n.std(G1,axis=0)
    mean_G2 = n.mean(G2,axis=0); std_G2 = n.std(G2,axis=0)
    state_1 = len(G1) ; state_2 = len(G2)

    Results["mean,std, and frequency of state 1"] = (mean_G1, std_G1, state_1)
    Results["mean,std, and frequency of state 2"] = (mean_G2, std_G2, state_2)

    df = pd.DataFrame(Results, index = ["Mean", "Standard deviation", "Frequency"])
    df.to_csv(output_path+".csv")

def main():
    """
    Main function

    It is called when this script is the main program called by the interpreter
    """

    # Parse command-line arguments
    reports, output_path, metric, threshold  = parseArgs()

    # Store the results in a CSV file
    CalculateFreeEnergy(reports, output_path, metric, threshold)


if __name__ == "__main__":
    """Call the main function"""
    main()
