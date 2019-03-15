# -*- coding: utf-8 -*-

import os,sys

# Script information

__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"

"""This program takes PDB files and
 process them to be able to be used
 in a PELE simulation"""

def help():
    """It prints the help text"""
    if sys.argv[1]=="-h" or sys.argv[1]=="--help":
       print("This python scripts modifies PDB files, the first command-line argument must be the filename of the input file without the .pdb extension format \n")
       exit()

def parseArgs():
    """Parse arguments from command-line

    RETURNS
    -------
    PDB_filename: string
    The data concerning the structure of the complex
    """
    PDB_filename=sys.argv[1]
    return PDB_filename

def PDB_processing(PDB_filename):
    """Opens the PDB file, modifies its content and overwrites it
    in order to be used in a PELE simulation.

    RETURNS
    -------
    PDB modified file
    """
    PDB_original,PDB_modified=open("%s.pdb"%PDB_filename,"rt"),open("%s_modified.pdb"%PDB_filename,"wt")
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
    os.system("rm %s.pdb" %PDB_filename)
    os.system("mv %s_modified.pdb %s.pdb" %(PDB_filename,PDB_filename))



def main():
    """Main function

    It is called when this script is the main program called by the interpreter
    """
    help()
    PDB_filename=parseArgs()
    PDB_processing(PDB_filename)

if __name__ == "__main__":
    """Call the main function"""
    main()

