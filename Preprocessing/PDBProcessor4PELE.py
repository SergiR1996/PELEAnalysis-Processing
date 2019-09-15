# -*- coding: utf-8 -*-

import os
import glob
import argparse as ap

# Script information

__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"

"""This program takes PDB files and
 process them to be able to be used
 in a PELE simulation"""

def storePDBfilenames(PDBs_to_parse, parser):
    """It identifies the PDB files to add to the PDBProcessor4PELE tool

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
    PDBs: list
            list of PDB files
    """

    parser = ap.ArgumentParser(description='Script used to clean PDB files for PELE simulation')
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", required=True, metavar="FILE",
                          type=str, nargs='*', help="path to PDB files")

    args = parser.parse_args()

    PDBs=storePDBfilenames(args.input, parser)

    return PDBs

def PDB_processing(PDB_filename):
    """Opens the PDB file, modifies its content and overwrites it
    in order to be used in a PELE simulation.

    RETURNS
    -------
    PDB modified file
    """

    PDB_original,PDB_modified=open("%s"%PDB_filename,"rt"),open("%s_modified.pdb"%PDB_filename[:-4],"wt")

    for line in PDB_original:
        if line.find("TER")!=-1:
            PDB_modified.write("TER\n")
        elif line.find("CONECT")!=-1:
            PDB_modified.write("TER\n")
            break
        elif line.find("ATOM")!=-1 or line.find("HETATM")!=-1:
            if line.find("H2  HOH")!=-1:
               PDB_modified.write(line+"TER\n")
            else:
               PDB_modified.write(line)
        else:
            continue

    PDB_modified.close()
    PDB_original.close()
    os.system("rm %s" %PDB_filename)
    os.system("mv %s_modified.pdb %s" %(PDB_filename[:-4],PDB_filename))



def main():
    """Main function

    It is called when this script is the main program called by the interpreter
    """

    PDBs=parseArgs()

    for PDB_filename in PDBs:
        PDB_processing(PDB_filename)

if __name__ == "__main__":
    """Call the main function"""
    main()

