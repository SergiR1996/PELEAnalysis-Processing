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
    """It identifies the reports to add to the ActiveSiteFinder tool

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
    reports : list
              list of PDB files
    """

    parser = ap.ArgumentParser()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", required=True, metavar="FILE",
                          type=str, nargs='*', help="path to PDB sfiles")
    args = parser.parse_args()
    PDBfiles = storePDBfilenames(args.input,parser)

    return PDBfiles

def PDBopener(PDBfilename):
    """ Take a PDB file, open it and store the position of Ser, His, Asp and Glu residues

    RETURNS
    -------
    Ser_list, His_list, Acid_list : lists
              list of the Ser, His, Asp and Glu residues
    """
    parser = PDBParser()
    structure = parser.get_structure(PDBfilename,PDBfilename)
    residues = structure.get_residues()
    Ser_list,His_list,Acid_list=[],[],[]
    for residue in residues:
        if residue.get_resname() == "SER":
            Ser_list.append(residue)
        elif residue.get_resname() == "HIS":
            His_list.append(residue)
        elif residue.get_resname() == "ASP" or residue.get_resname() == "GLU":
                Acid_list.append(residue)

    return Ser_list,His_list,Acid_list

def ActiveSiteFinder(PDBfilename,Ser_list,His_list,Acid_list):
    Results = {}
    for residue_H in His_list:
        for residue_S in Ser_list:
            if abs(residue_S["HG"]-residue_H["NE2"])<=3:
                for residue_A in Acid_list:
                    if residue_A.get_resname()=="ASP":
                        if (residue_A["OD1"]-residue_H["HD1"])<=2.4 or (residue_A["OD2"]-residue_H["HD1"])<=2.4:
                            Results[PDBfilename] = [residue_S.id[1],residue_H.id[1],residue_A.id[1],0]
                    elif residue_A.get_resname()=="GLU":
                        if (residue_A["OE1"] - residue_H["HD1"]) <= 2.4 or (residue_A["OE2"] - residue_H["HD1"]) <= 2.4:
                            Results[PDBfilename] = [residue_S.id[1],residue_H.id[1],residue_A.id[1],1]

    return Results




def main():
    """Main function

    It is called when this script is the main program called by the interpreter
    """

    # Store PDB filenames on a list
    PDB_List = parseArgs()

    # Get the Ser, His, Asp and/or Glu residues and store the ones that are in the active site in the output file.
    Dict_file=open("Dict_file","wt")
    i=0
    for PDB in PDB_List:
        Ser_List,His_List,Acid_List = PDBopener(PDB)
        Results = ActiveSiteFinder(PDB,Ser_List,His_List,Acid_List)
        if i==0:
            for PDB,residues in Results.items():
                Dict_file.write("Protein dict = {"+ '"' + str(PDB) + '"' + ": " + str(residues) + ", ")
        elif i==(len(PDB_List)-1):
            for PDB,residues in Results.items():
                Dict_file.write('"' + str(PDB) + '"' + ": " + str(residues) +"}")
        else:
            for PDB,residues in Results.items():
                Dict_file.write('"' + str(PDB) + '"' + ": " + str(residues) + ", ")
        i+=1



if __name__ == "__main__":
    """Call the main function"""
    main()

