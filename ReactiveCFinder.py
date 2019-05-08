# -*- coding: utf-8 -*-


# Imports 
import os,sys
import glob
from Bio.PDB import *
import argparse as ap

# Script information
__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"

def storePDBfilenames(PDBs_to_parse,parser):
    """It identifies the reports to add to the ReactiveCFinder tool

    PARAMETERS
    ----------
    reports_to_parse : list of strings
                       all the PDB files that want to be added to the analysis

    RETURNS
    -------
    parsed_data : list of PDB filenames (strings)
    """

    reports = []

    for reports_list in PDBs_to_parse:
        PDB_found = glob.glob(reports_list)
        if len(PDB_found) == 0:
            print("Warning: path to report file \'" +
                  "{}".format(reports_list) + "\' not found.")
        for report in PDB_found:
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
              list of PDB files
    """

    parser = ap.ArgumentParser()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", required=True, metavar="FILE",
                          type=str, nargs='*', help="path to PDB sfiles")
    args = parser.parse_args()
    PDBfiles = storePDBfilenames(args.input,parser)

    return PDBfiles

def ReactiveCFinder(PDBfilename):
    """ Take a PDB file and store the reactive C on a dictionary

    RETURNS
    -------
    Results : dict
              dictionary of ester C, which are the reactive ones
    """
    parser,Results = PDBParser(),{}
    structure = parser.get_structure(PDBfilename,PDBfilename)
    for atom1 in structure.get_atoms():
        for atom2 in structure.get_atoms():
                if "C" in atom1.get_name() and "O" in atom2.get_name():
                    if abs(atom1-atom2)<=1.25:
                        for atom3 in structure.get_atoms():
                            for atom4 in structure.get_atoms():
                                if "O" in atom3.get_name() and 1.25<abs(atom1-atom3)<1.4 \
                                        and "C" in atom4.get_name() and 1.25<abs(atom3-atom4)<1.4 and \
                                        atom4.get_name()!=atom1.get_name():
                                    if PDBfilename not in Results:
                                        Results[PDBfilename]=[]
                                    if atom1.get_fullname() not in Results[PDBfilename]:
                                        Results[PDBfilename].append(atom1.get_fullname())

    return Results




def main():
    """Main function

    It is called when this script is the main program called by the interpreter
    """

    # Store PDB filenames on a list
    PDB_List = parseArgs()

    # Get the reactive C and store them in the output file.
    Dict_file=open("Ligand_file","wt")
    i=0
    for PDB in PDB_List:
        Results = ReactiveCFinder(PDB)
        if i==0:
            for PDB,reactive_C in Results.items():
                Dict_file.write("Protein dict = {" + '"' + str(PDB) + '"' + ": " + str(reactive_C) + ", ")
        elif i==(len(PDB_List)-1):
            for PDB,reactive_C in Results.items():
                Dict_file.write('"' + str(PDB) + '"' + ": " + str(reactive_C) + "}")
        else:
            for PDB,reactive_C in Results.items():
                Dict_file.write('"' + str(PDB) + '"' + ": " + str(reactive_C) + ", ")
        i+=1



if __name__ == "__main__":
    """Call the main function"""
    main()

