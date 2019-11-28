# -*- coding: utf-8 -*-

# Global imports
import glob
import argparse as ap
import pandas as pd

# Script information
__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"

def storePDBfilenames(PDBs_to_parse,parser):
    """
    It identifies the reports to add to the protein preparation system

    PARAMETERS
    ----------
    PDBs_to_parse : list of strings
                       all the PDB files that want to be added to the analysis

    RETURNS
    -------
    parsed_data : list of PDB filenames (strings)
    """

    PDBs = []

    for PDBs_list in PDBs_to_parse:
        PDB_found = glob.glob(PDBs_list)
        if len(PDB_found) == 0:
            print("Warning: path to report file \'" +
                  "{}".format(PDBs_list) + "\' not found.")
        for PDB in PDB_found:
            PDBs.append(PDB)

    if len(PDBs) == 0:
        print("Error: list of report files is empty.")
        parser.print_help()
        exit(1)

    return PDBs

def parseArgs():
    """
    Parse arguments from command-line

    RETURNS
    -------
    PDBfiles : list
              list of PDB files
    """

    parser = ap.ArgumentParser()
    required = parser.add_argument_group('required arguments')
    optional = parser._action_groups.pop()
    required.add_argument("-i", "--input", required=True, metavar="FILE",
                          type=str, nargs='*', help="path to PDB files")
    required.add_argument("-c","--csv", required=True, metavar="FILE",
                        type=str,help="path of csv file with the other protein features")
    optional.add_argument("-CP","--computeprops"
                          ,help = "Compute the SiteMap descriptors", action = "store_true")
    args = parser.parse_args()

    PDBfiles = storePDBfilenames(args.input,parser)
    csv = args.csv
    computeprops = args.computeprops

    return PDBfiles, csv, computeprops