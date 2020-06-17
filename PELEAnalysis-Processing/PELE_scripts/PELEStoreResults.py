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
    catalytic_event : list of indices
                  list of column indices
    serine_substrate: float
                  threshold for the serine-substrate distance
    to_drop : list of strings
                  list of column names that want to be dropped
    """

    parser = ap.ArgumentParser(description='Script that returns a csv file with the mean of the numerical \
        metrics of the reports file from a PELE simulation. It also can return the number of catalytic catalytic events \
        and catalytic trajectories in the PELE simulation.')
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", required=True, metavar="FILE",
                          type=str, nargs='*', help="path to report files")
    optional.add_argument("-o", "--output", metavar="PATH", type=str,
                          help="output path to save figure", default="PELE_results")
    optional.add_argument("-CE","--catalytic_event", metavar="LIST",type=str,
                          nargs='*',help="index of the column where catalytic distances reside (they must be 3)")
    optional.add_argument("-SS","--serine_substrate", metavar="THRESHOLD",type=float,
                          help="threshold for the Serine-substrate distance", default = 3.5)
    optional.add_argument("-TD","--to_drop", metavar="LIST",type=str,
                          nargs='*',help="column names that want to be dropped", default=[])
    parser._action_groups.append(optional)
    args = parser.parse_args()

    reports = parseReports(args.input, parser)

    output_path, catalytic_event, to_drop, serine_substrate = args.output, args.catalytic_event, args.to_drop, args.serine_substrate

    if catalytic_event is not None:
        catalytic_event = [int(i)-1 for i in catalytic_event]

    return reports, output_path, catalytic_event, to_drop, serine_substrate

def Storeresults(reports,catalytic_event, output_path, to_drop, serine_substrate):
    """
    Take the PELE simulation report files and returns the results stored in a CSV file

    OUTPUT
    ------
    output_path.csv file with the average metrics of the PELE simulation.
    """

    Results,means_aux, cat_events,cat_trajectories = {},[], 0, 0

    for i,report in enumerate(reports):
        rep = pd.read_csv(report,sep="    ")
        means_aux.append(list(n.mean(rep,axis=0))[3:])
        CATE = rep[(rep[rep.columns[catalytic_event[0]]] <= serine_substrate) & (rep[rep.columns[catalytic_event[1]]] <= 3.5) & ((rep[rep.columns[catalytic_event[2]]] <= 3.5) | (rep[rep.columns[catalytic_event[2]+1]] <= 3.5))]
        CE = CATE.shape[0]
        cat_events += CE
        if CE!=0:
            cat_trajectories += 1
            print("{} --> {}".format(report,list(CATE["numberOfAcceptedPeleSteps"])))
        if i==0:
            column_names = list(rep.columns[3:])


    means = n.mean(means_aux,axis=0)
    std = n.std(means_aux,axis=0)
    means, std = [round(i,3) for i in means],[round(i,3) for i in std]
    means_std = [];column_names.append("cat_events");column_names.append("cat_n_trajectories")

    for i in range(len(means)):
        means_std.append(str(means[i])+"(+-)"+str(std[i]))

    means_std.append(cat_events);means_std.append(cat_trajectories)

    for key,item in zip(column_names,means_std):
        Results[key] = item

    df = pd.DataFrame(Results,index=["Report_summary"])
    df.drop(to_drop, axis = 1, inplace=True)
    df = df.round(3)
    df.to_csv(output_path+".csv")



    #     reportID = int(os.path.basename(report).split('_')[-1].split('.')[0])
    #     with open(report,'r') as report_file:

    #         cat_counting = 0
    #         if '0' not in Results:
    #             line = report_file.readline()
    #             Results = dict.fromkeys(line.split("    ")[3:(len(line.split()))])
    #             Means = [[] for x in range(3,(len(line.split())))]
    #         else:
    #             next(report_file)

    #         for i, line in enumerate(report_file):

    #             if catalytic_event is not None:
    #                 if float(line.split()[catalytic_event[0]])<=3.5 and float(line.split()[catalytic_event[1]])<=3.5 and (float(line.split()[catalytic_event[2]])<=3.5 or float(line.split()[catalytic_event[2]+1])<=3.5):
    #                     cat_counting+=1

    #             for element,i in zip(line.split(),range(len(line.split()))):
    #                 if i>=3:
    #                     Means[i-3].append(float(element))

    #         means_aux.append(n.array([n.mean(Means[i]) for i in range(len(Means))]))
    #         cat_events.append((cat_counting,reportID))

    # cat_events = n.array(cat_events)
    # means = n.mean(means_aux,axis=0)
    # std = n.std(means_aux,axis=0)
    # Results.pop("\n",None)

    # for i,elem in enumerate(Results.keys()):
    #     Results[elem] = str(round(means[i],3))+"(+-)"+str(round(std[i],3))

    # if catalytic_event is not None:
    #     cat_ev = []
    #     Results["cat_events"] = n.sum(cat_events,axis=0)[0]

    #     for elem in cat_events:
    #         if elem[0]!=0:
    #             cat_ev.append(elem)

    #     Results["cat_n_trajectories"] = len(cat_ev)

    # return Results

# def Outputresults(Results,output_path):
#     """Take the PELE simulation results and report the means with a csv table

#     RETURNS
#     -------
#     Results: dictionary of lists
#              dictionary containing the mean of the different quantitative parameters
#     """

#     df = pd.DataFrame(Results,index=["Report_summary"])
#     df.drop(["currentEnergy"], axis = 1, inplace=True)

#     df.to_csv(output_path+".csv")

                

def main():
    """
    Main function

    It is called when this script is the main program called by the interpreter
    """

    # Parse command-line arguments
    reports, output_path, catalytic_event, to_drop, serine_substrate  = parseArgs()

    # Store the results in a CSV file
    Storeresults(reports,catalytic_event, output_path, to_drop, serine_substrate)


if __name__ == "__main__":
    """Call the main function"""
    main()

