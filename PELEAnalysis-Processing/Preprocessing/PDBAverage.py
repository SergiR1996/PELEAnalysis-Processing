# -*- coding: utf-8 -*-

# Global imports
import os, time
import glob
import argparse as ap
import numpy as np

# Local imports
from StorePDBFilenames import *

# Script information
__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"

def parseArgs():
    """
    Parse arguments from command-line.

    RETURNS
    -------
    PDBs: list
            list of PDB files
    """

    parser = ap.ArgumentParser(description='Script used to average the coordinates from a set\
        of PDB files')
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", required=True, metavar="PDB FILE",
                          type=str, nargs='*', help="path to PDB files")
    optional.add_argument("-o", "--output", metavar="PDB FILE", type=str,
                          help="the output PDB filename", default="average")
    optional.add_argument("-B", "--B_factor", help="Output the calculated\
        B_factor out of the average of coordinates", action="store_true")
    parser._action_groups.append(optional)

    args = parser.parse_args()

    PDBs=storePDBfilenames(args.input, parser)

    return PDBs, args.output, args.B_factor

def getBFactor(atid_coords):
    """
    Estimate the B-factor from the root mean square fluctuation using
    the average from the trajectories or list 
    of PDB structurally aligned files.

    RETURNS
    -------
    B_factor: string
            Number of the temperature factor in string format
    """

    MSF = 0

    for list_of_coordinates in atid_coords.values():
        for coordinate in list_of_coordinates:
            MSF += (abs(coordinate - np.average(list_of_coordinates)) ** 2) / len(list_of_coordinates)

    B_factor = str(np.round(MSF * (8 * np.pi ** 2) / 3, 2))
    if len(B_factor.split(".")[1])<2:
        B_factor += (2 - len(B_factor.split(".")[1])) * "0"

    return B_factor

def averageStructure(PDBs, output, B_factor):
    """
    Function that calculates the average coordinates and 
    outputs the "average" PDB file. The B-factor can be 
    added to the atoms or not.

    OUTPUT
    ------
    The average PDB file with the calculated B-factor or not
    """
    
    start = time.time()
    D_av = {}
    for PDB_filename in PDBs:
        for atom in open(PDB_filename):
            atid = atom[6:11]
            if atid not in D_av:
                D_av[atid] = {"x": [], "y": [], "z": []}
            if atom[0:4] == "ATOM":
                x = float(atom[30:38].strip(" "))
                y = float(atom[38:46].strip(" "))
                z = float(atom[46:54].strip(" "))
                D_av[atid]["x"].append(x)
                D_av[atid]["y"].append(y)
                D_av[atid]["z"].append(z)


    average_pdb = open(output + ".pdb","w")
    for atom in open(PDBs[0]):
        atid = atom[6:11]
        if atom[0:4] == "ATOM":
            x = str(round(np.average(D_av[atid]["x"]),3))
            y = str(round(np.average(D_av[atid]["y"]),3))
            z = str(round(np.average(D_av[atid]["z"]),3))
            if len(x.split(".")[1])<3:
                x += (3 - len(x.split(".")[1])) * "0"
            if len(y.split(".")[1])<3:
                y += (3 - len(y.split(".")[1])) * "0"
            if len(z.split(".")[1])<3:
                z += (3 - len(z.split(".")[1])) * "0"
            if B_factor:
                B_factor_value = getBFactor(D_av[atid])
                average_pdb.write(atom[0:30]+"{:>8}".format(x)+"{:>8}".format(y)+"{:>8}".format(z)
                +atom[54:60]+"{:>6}".format(B_factor_value)+atom[66:])
            else:
                average_pdb.write(atom[0:30]+"{:>8}".format(x)+"{:>8}".format(y)+"{:>8}".format(z)
                +atom[54:])
        elif atom[0:3] == "TER" or atom[0:3] == "END" or atom[0:6] == "HETATM":
            average_pdb.write("TER")
            break
        else:
            continue

    end = time.time()
    print("\nThe main code needed {} seconds to compute the average PDB file \n".format(end-start))

def main():
    """
    Main function

    It is called when this script is the main program called by the interpreter.
    """

    # Parse command-line arguments
    PDBs, output, B_factor = parseArgs()

    # Average PDB structure from input PDB files
    averageStructure(PDBs, output, B_factor)

if __name__=="__main__":
	"""Call the main function"""
	main()
