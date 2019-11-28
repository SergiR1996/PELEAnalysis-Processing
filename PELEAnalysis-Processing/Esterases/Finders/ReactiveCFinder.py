# -*- coding: utf-8 -*-


# Imports 
import os,sys
import glob
import math
import argparse as ap

# Script information
__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"

def storePDBfilenames(PDBs_to_parse, parser):
    """It identifies the PDB files to add to the ReactiveCFinder tool

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
    reports : string
              list of PDB files
    """

    parser = ap.ArgumentParser(description='Script used to find the reactive C of an ester in a ligand PDB file')
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

    def ComputeDistance(atom1,atom2):
        """Computes the module or distance between the two points"""

        r = [atom2[0] - atom1[0], atom2[1] - atom1[1], atom2[2] - atom1[2]]
        return math.sqrt(r[0] ** 2 + r[1] ** 2 + r[2]**2)

    PDB = open(PDBfilename,"r")
    lines = PDB.readlines()
    Results = {}

    for atom1 in lines:
        for atom2 in lines:
            if "HETATM" in atom1 and "HETATM" in atom2 and "C" in atom1[12:16].strip() and "O" in atom2[12:16].strip():
                x1,y1,z1 = float(atom1[30:38].strip()),float(atom1[38:46].strip()),float(atom1[46:54].strip())
                x2,y2,z2 = float(atom2[30:38].strip()),float(atom2[38:46].strip()),float(atom2[46:54].strip())
                coord1,coord2 = [x1,y1,z1],[x2,y2,z2]
                if ComputeDistance(coord1,coord2)<=1.25:
                    for atom3 in lines:
                        for atom4 in lines:
                            if "HETATM" in atom3 and "HETATM" in atom4:
                                x3,y3,z3 = float(atom3[30:38].strip()),float(atom3[38:46].strip()),float(atom3[46:54].strip())
                                x4,y4,z4 = float(atom4[30:38].strip()),float(atom4[38:46].strip()),float(atom4[46:54].strip())
                                coord3,coord4 = [x3,y3,z3],[x4,y4,z4]
                                if "O" in atom3[12:16].strip() and "C" in atom4[12:16].strip() \
                                and 1.2<ComputeDistance(coord1,coord3)<1.525 and 1.2<ComputeDistance(coord1,coord4)<1.525 \
                                and atom4[12:16].strip()!=atom1[12:16].strip():
                                    if PDBfilename not in Results:
                                        Results[PDBfilename]=[]
                                    if atom1[12:16].strip() not in Results[PDBfilename]:
                                        Results[PDBfilename].append(atom1[12:16].strip())

    return Results

def main():
    """Main function

    It is called when this script is the main program called by the interpreter
    """

    # Store PDB filenames on a list
    PDB_List = parseArgs()

    # Get the reactive C and store them in the output file.
    Dict_file=open("Ligand_file","wt")

    for PDB in PDB_List:
        Results = ReactiveCFinder(PDB)
        print(Results)
        if PDB==PDB_List[0]:
            for PDB,reactive_C in Results.items():
                Dict_file.write("Ligand_dict = {" + '"' + str(PDB) + '"' + ": " + str(reactive_C) + ", ")
        elif PDB==PDB_List[-1]:
            for PDB,reactive_C in Results.items():
                Dict_file.write('"' + str(PDB) + '"' + ": " + str(reactive_C) + "}")
        else:
            for PDB,reactive_C in Results.items():
                Dict_file.write('"' + str(PDB) + '"' + ": " + str(reactive_C) + ", ")



if __name__ == "__main__":
    """Call the main function"""
    main()

