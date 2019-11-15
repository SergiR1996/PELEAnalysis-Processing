# -*- coding: utf-8 -*-


# Imports
from __future__ import unicode_literals
import os
import glob

def parseReports(reports_to_parse, parser):
    """It identifies the reports to perform some analysis

    PARAMETERS
    ----------
    reports_to_parse : list of strings
                       all the report files that want to be added to the analysis
    parser : ArgumentParser object
             contains information about the command line arguments

    RETURNS
    -------
    parsed_data : tuple of a list and a string
                  the list specifies the report files that will be added to the analysis.
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
