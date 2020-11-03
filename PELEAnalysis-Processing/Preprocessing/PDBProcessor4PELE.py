# -*- coding: utf-8 -*-

# Global imports
import os,re
import glob
import argparse as ap

# Local imports
from StorePDBFilenames import *

# Script information
__author__ = "Sergi Rodà"
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Sergi Rodà"
__email__ = "sergi.rodallordes@bsc.es"

# Global variables
Protein_list = ["GLY","ALA","VAL","LEU","ILE","SER","THR","ARG","LYS","PHE","TYR","TRP","ASP","GLU","ASN",
"GLN","PRO","CYS","CYX","MET","HIS","HID","HIE","HIP","HOH","GLH","ASH","LYN"]

"""This program takes PDB files and
 process them to be able to be used
 in a PELE simulation"""

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
    optional.add_argument("-RN", "--residue_name", metavar = "LIST",
                          nargs="*",type = str, help = "residue name of modified amino acids", default = [])    
    parser._action_groups.append(optional)
    args = parser.parse_args()

    PDBs = storePDBfilenames(args.input, parser)
    Output = args.output

    return PDBs, Output, args.residue_name

def PDB_processing(PDB_filename, Output, Residue_name):
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

    Non_aminoacid_dict, L_number = {},0

    # Depending on if a output filename is specified on the command line or not, the writable files will be different.
    if Output != "":
        PDB_original, PDB_modified = open("{}".format(PDB_filename),"rt"), open("{}.pdb".format(Output),"wt")

    else:
        PDB_original, PDB_modified = open("{}".format(PDB_filename),"rt"), open("{}_modified.pdb".format(PDB_filename[:-4]),"wt")

    Lines = PDB_original.readlines() # Store the lines of the input PDB file in the Lines variable
    i = 0 # To iterate over the lines of the PDB file
    water_i,water_j = 0,0 # To count the number of water molecules and rewrite the PDB atom names
    while i < len(Lines):

        line = Lines[i]
        if (line[0:3] == "TER" != -1 or line[0:3] == "END"): # Preserve the TER or END lines
            PDB_modified.write("TER\n")

        elif line.find("CONECT") != -1: # If CONECT is found in the line, it means the rest of the file is useless for PELE.
            PDB_modified.write("TER\n")
            break

        elif (line[0:4] == "ATOM"  or line[0:6] == "HETATM"): # These conditional statements control the TER between molecules and also the remaining coordinate lines.

            #if (line[21:22].strip() != Lines[i+1][21:22].strip()) and (line[21:22].strip() != "" and Lines[i+1][21:22].strip() != ""):
            #    print(Lines[i+1][21:22].strip())
            #    print(line.strip("\n"))
            #    PDB_modified.write(line+"TER\n")
            #    i += 1
            #    continue
            if (line.find("2HW  HOH") != -1 or line.find("H2  HOH") != -1) and Lines[i+1].find("CONECT") != -1:
                PDB_modified.write(line[0:12] + "2HW  HOH W{:>4}".format(water_i) + line[26:] + "TER\n")
                break
            elif line.find("2HW  HOH") != -1 or line.find("H2  HOH") != -1:
                PDB_modified.write(line[0:12] + "2HW  HOH W{:>4}".format(water_i) + line[26:] + "TER\n")
                water_i += 1; water_j = 0
            else:
                if line[17:20].strip() not in Protein_list:
                    if line[17:20].strip() not in Non_aminoacid_dict:
                        if line[17:20].strip() in Residue_name and len(Residue_name)!=0:
                            PDB_modified.write(line)
                        else:
                            Non_aminoacid_dict[line[17:20].strip()] = 1
                            L_number+=1
                            if Lines[i-1][0:3]!="TER":
                                PDB_modified.write("TER\n")
                            PDB_modified.write(line[0:21]+"L{:>4}".format(L_number)+line[26:])
                    elif line[17:20].strip() != Lines[i-1][17:20].strip():
                        Non_aminoacid_dict[line.strip()[17:20]] += 1
                        L_number+=1
                        if Lines[i-1][0:3]!="TER":
                            PDB_modified.write("TER\n")
                        PDB_modified.write(line[0:21]+"L{:>4}".format(L_number)+line[26:])
                    else:
                        PDB_modified.write(line[0:21]+"L{:>4}".format(L_number)+line[26:])
                    i += 1 
                    continue         
                if line[17:20].strip() =="HOH" and Lines[i-1][17:20].strip()!="HOH":
                    water_i += 1
                    if Lines[i-1][0:3]!="TER":
                        PDB_modified.write("TER\n")
                    PDB_modified.write(line[0:12]+" OW  HOH W{:>4}".format(water_i)+line[26:])
                    water_j += 1
                elif line[17:20].strip() =="HOH" and water_j == 0:
                    PDB_modified.write(line[0:12]+" OW  HOH W{:>4}".format(water_i)+line[26:])
                    water_j += 1
                elif line[17:20].strip() =="HOH" and water_j == 1:
                    PDB_modified.write(line[0:12]+"1HW  HOH W{:>4}".format(water_i)+line[26:])
                else:
                    if line[12:16] == " HXT":
                        PDB_modified.write(line[0:12]+" OXT"+line[16:77]+"O1-\n")
                    else:
                        PDB_modified.write(line)
        else:
            pass

        i += 1

    PDB_modified.close()
    PDB_original.close()

    return water_i, Non_aminoacid_dict

def Printing_summary(PDB_filename, Output, water_i, Non_aminoacid_dict):
    """
    Prints a summary of the processing of the PDB file prior to the PELE
    simulation. It includes the number of water molecules found and the 
    unconventional molecules found.

    PARAMETERS
    ----------
    PDB_filename : string
                      filename of the input PDB file that wants to be processed
    Output : string
                      filename of the output PDB file after processing
    water_i : integer
                      number of water molecules found in the PDB input file
    Non_aminoacid_dict : dictionary
                      dicitonary with the key:item pair being the unconventional AA and the number of instances                

    RETURNS
    -------
    PDB modified file
    """

    Screen_ticks = 10
    
    # Depending on if a output filename is specified on the command line or not, the input PDB file will be preserved or overwritten.
    if Output != "":
        print("\n{} has been succesfully processed and written to {}.pdb\n\n".format(PDB_filename, Output) + "-"*Screen_ticks + "SUMMARY" + "-"*Screen_ticks)
        if water_i == 0:
            print("\nNo water molecules were found in the PDB file")
        elif water_i == 1:
            print("\n{} water molecule was found in the PDB file".format(water_i))
        else:
            print("\n{} water molecules were found in the PDB file".format(water_i))
        if Non_aminoacid_dict != {}:
            print("\nThe following ligands and cofactors were found in the PDB file\n")
            for ligand in Non_aminoacid_dict.keys():
                print(ligand + ": " + str(Non_aminoacid_dict[ligand]) + " molecule/s\n")

    else:
        os.system("rm {}" .format(PDB_filename))
        os.system("mv {}_modified.pdb {}" .format(PDB_filename[:-4], PDB_filename))
        print("\n{} has been succesfully processed and overwritten\n\n" .format(PDB_filename) + "-"*Screen_ticks + "SUMMARY" + "-"*Screen_ticks)
        print("\n{} water molecules were found in the PDB file" .format(water_i))
        if Non_aminoacid_dict != {}:
            print("\nThe following ligands and cofactors were found in the PDB file\n")
            for ligand in Non_aminoacid_dict.keys():
                print(ligand + ": " + str(Non_aminoacid_dict[ligand]) + " molecule/s\n")



def main():
    """
    Main function

    It is called when this script is the main program called by the interpreter
    """

    PDBs, Output, Residue_name = parseArgs()

    for PDB_filename in PDBs:
        water_i, Non_aminoacid_dict = PDB_processing(PDB_filename, Output, Residue_name)
        Printing_summary(PDB_filename, Output, water_i, Non_aminoacid_dict)

if __name__ == "__main__":
    """Call the main function"""
    main()
