# -*- coding: utf-8 -*-


# Imports
import os,sys,re
import glob
import argparse as ap

# Script information
__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"

def storePDBfilenames(PDBs_to_parse, parser):
    """It identifies the PDB files to add to the Checkmissingligands tool

    PARAMETERS
    ----------
    PDBs_to_parse : list of strings
                       all the PDB files that want to be added to the analysis

    RETURNS
    -------
    parsed_data : list of PDB filenames (strings)
    """

    PDBs = []

    for PDB_list in PDBs_to_parse:
        PDB_found = glob.glob(PDB_list)
        if len(PDB_found) == 0:
            print("Warning: path to PDB file \'" +
                  "{}".format(PDB_list) + "\' not found.")
        for PDB in PDB_found:
            PDBs.append(PDB)

    if len(PDBs) == 0:
        print("Error: list of PDB files is empty.")
        parser.print_help()
        exit(1)

    return PDBs

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
                          type=str, nargs='*', help="path to PDB files")
    required.add_argument("-f", "--file", required=True, metavar="FILE",
                          type=str, nargs='*', help="path to csv file")
    args = parser.parse_args()
    PDBfiles = storePDBfilenames(args.input,parser)
    CSVfile = args.file

    return PDBfiles,CSVfile[0]

def preprocessligpdb(Ligand_list):
    """It takes the list of ligand PDB filenames and process them for the matching

    PARAMETERS
    ----------
    Ligand_list : list of strings
                       list of ligand PDB filenames

    RETURNS
    -------
    lig_aux : list of ligand PDB filenames processed
    """

    lig_aux = []

    for ligand in Ligand_list:
        lig_aux.append(ligand.replace("ligands/","").replace(".pdb",""))

    return lig_aux

def Storeligandnames(csv_file):
    """It identifies the names of the ligands in the csv file

    PARAMETERS
    ----------
    csv_file : filename of the csv file with the ligands

    RETURNS
    -------
    lig_list : list of ligand names (list of strings)
    """

    Lig = open(csv_file,"rt")
    lig_aux = []

    for ligand in Lig:
        lig_aux.append(ligand.replace(" ","_").replace("\n","").lower())

    return lig_aux

def Checkmissingligands(Lig_list,Lig_aux):
    """It checks the ligand PDB files that are missing

    PARAMETERS
    ----------
    Lig_list : List with the ligand PDB filenames
    Lig_aux : List with the ligand names

    RETURNS
    -------
    Missing_lig : List ot the name of the missing ligand PDB files

    """


    Missing_lig = []

    for ligand in Lig_aux:
        if ligand not in Lig_list:
            Missing_lig.append(ligand)

    return Missing_lig



def main():
    """Main function

    It is called when this script is the main program called by the interpreter
    """

    # Store the ligands PDB files and the CSV file with the ligand names
    Ligand_list,csv_file = parseArgs()

    # Preprocess the ligand pdb filenames to match them with the csv file
    Ligand_list = preprocessligpdb(Ligand_list)

    # Store the ligand names in a list of strings
    Lig_aux = Storeligandnames(csv_file)

    # Display the missing PDB ligand files
    print(Checkmissingligands(Ligand_list,Lig_aux))



if __name__=="__main__":
    """Call the main function"""
    main()