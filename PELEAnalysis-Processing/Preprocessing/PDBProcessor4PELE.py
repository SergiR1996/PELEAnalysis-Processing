# -*- coding: utf-8 -*-

# Global imports
import os,re
import glob
import argparse as ap

# Script information

__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"

# Global variables
Protein_list = ["GLY","ALA","VAL","LEU","ILE","SER","THR","ARG","LYS","PHE","TYR","TRP","ASP","GLU","ASN","GLN","PRO","CYS","CYX","MET","HIS","HID","HIE","HIP","HOH"]

"""This program takes PDB files and
 process them to be able to be used
 in a PELE simulation"""

def storePDBfilenames(PDBs_to_parse, parser):
    """
    It identifies the PDB files to add to the PDBProcessor4PELE tool

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
    """
    Parse arguments from command-line

    RETURNS
    -------
    PDBs: list
            list of PDB files
    """

    parser = ap.ArgumentParser(description = "Script used to clean PDB files for PELE simulation")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    required.add_argument("-i", "--input", required = True, metavar = "FILE",
                          type = str, nargs = '*', help = "path to PDB files")
    optional.add_argument("-O", "--output", metavar = "STRING",
                          type = str, help = "filename for the output proecessed PDB file", default = "")
    parser._action_groups.append(optional)
    args = parser.parse_args()

    PDBs = storePDBfilenames(args.input, parser)
    Output = args.output

    return PDBs, Output

def PDB_processing(PDB_filename, Output):
    """
    Opens the PDB file, modifies its content and overwrites it
    in order to be used in a PELE simulation.

    PARAMETERS
    ----------
    PDB_filename : string
                      filename of the input PDB file that wants to be processed
    Output : string
                      filename of the output PDB file after processing

    RETURNS
    -------
    PDB modified file
    """

    # Depending on if a output filename is specified on the command line or not, the writable files will be different.
    if Output != "":
        PDB_original, PDB_modified = open("{}".format(PDB_filename),"rt"), open("{}.pdb".format(Output),"wt")

    else:
        PDB_original, PDB_modified = open("{}".format(PDB_filename),"rt"), open("{}_modified.pdb".format(PDB_filename[:-4]),"wt")

    Lines = PDB_original.readlines() # Store the lines of the input PDB file in the Lines variable

    i = 0 # To iterate over the lines of the PDB file
    water_i,water_j = 1,0 # To count the number of water molecules and rewrite the PDB atom names
    while i < len(Lines):

        line = Lines[i]
        if line.find("TER") != -1 or line.find("END") != -1: # Preserve the TER or END lines
            PDB_modified.write("TER\n")

        elif line.find("CONECT") != -1: # If CONECT is found in the line, it means the rest of the file is useless for PELE.
            PDB_modified.write("TER\n")
            break

        elif line.find("ATOM") != -1 or line.find("HETATM") != -1: # These conditional statements control the TER between molecules and also the remaining coordinate lines.

            if (line.strip()[21:22] != Lines[i+1].strip()[21:22]) and (line.strip()[21:22] != "" and Lines[i+1].strip()[21:22] != ""):
                PDB_modified.write(line+"TER\n")
                i += 1
                continue
            if (line.find("2HW  HOH") != -1 or line.find("H2  HOH") != -1) and Lines[i+1].find("CONECT") != -1:
                PDB_modified.write(line[0:12] + "2HW  HOH W{:>4}".format(water_i) + line[26:] + "TER\n")
                break
            elif line.find("2HW  HOH") != -1 or line.find("H2  HOH") != -1:
                PDB_modified.write(line[0:12] + "2HW  HOH W{:>4}".format(water_i) + line[26:] + "TER\n")
                water_i += 1; water_j = 0
            else:

                if line.strip()[17:20] not in Protein_list:
                    PDB_modified.write(line[0:21]+"L   1"+line[26:])
                    i += 1 
                    continue
                if line.strip()[17:20] =="HOH" and Lines[i-1].strip()[17:20] not in Protein_list:
                    PDB_modified.write("TER\n"+line[0:12]+" OW  HOH W{:>4}".format(water_i)+line[26:])
                    water_j += 1
                elif line.strip()[17:20] =="HOH" and water_j == 0:
                    PDB_modified.write(line[0:12]+" OW  HOH W{:>4}".format(water_i)+line[26:])
                    water_j += 1
                elif line.strip()[17:20] =="HOH" and water_j == 1:
                    PDB_modified.write(line[0:12]+"1HW  HOH W{:>4}".format(water_i)+line[26:])
                else:
                    PDB_modified.write(line)
        else:
            pass

        i += 1

    PDB_modified.close()
    PDB_original.close()
    
    # Depending on if a output filename is specified on the command line or not, the input PDB file will be preserved or overwritten.
    if Output != "":
        print("{} has been succesfully processed and written to {}.pdb" .format((PDB_filename, Output)))

    else:
        os.system("rm {}" .format(PDB_filename))
        os.system("mv {}_modified.pdb {}" .format((PDB_filename[:-4], PDB_filename)))
        print("{} has been succesfully processed and overwritten" .format(PDB_filename))



def main():
    """
    Main function

    It is called when this script is the main program called by the interpreter
    """

    PDBs, Output = parseArgs()

    for PDB_filename in PDBs:
        PDB_processing(PDB_filename, Output)

if __name__ == "__main__":
    """Call the main function"""
    main()