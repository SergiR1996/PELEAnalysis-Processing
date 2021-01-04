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
    HXT_terminal, Gly_to_other_residue, other_residue_to_Gly, ASH_list, GLH_list = False, False, False, False, False

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
                    # Conditional statement used when C-ter residue contains the HXT atom name and type automatically given by Maestro/Schrödinger
                    if line[12:16] == " HXT":
                        #print(line)
                        PDB_modified.write(line[0:12]+" OXT"+line[16:77]+"O1-\n")
                        HXT_terminal = True
                    # Conditional statement used when residues are mutated from Gly to another thing with Maestro/Schrödinger
                    elif line[12:16] == " HA2" and line[17:20] != "GLY":
                        PDB_modified.write(line[0:12]+" HA "+line[16:])
                        if Gly_to_other_residue == False:
                            Gly_to_other_residue = list()
                        if line[22:26].strip() not in Gly_to_other_residue:
                            Gly_to_other_residue.append(line[22:26].strip())
                    # Conditional statement used when residues are mutated to Gly with Maestro/Schrödinger
                    elif line[12:16] == " HA " and line[17:20] == "GLY":
                        PDB_modified.write(line[0:12]+" HA2"+line[16:])
                        if other_residue_to_Gly == False:
                            other_residue_to_Gly = list()
                        if line[22:26].strip() not in other_residue_to_Gly:
                            other_residue_to_Gly.append(line[22:26].strip())
                    # Conditional statement used when Asp residue is protonated in the carboxyl group
                    elif line[12:16] == " HD2" and line[17:20] == "ASP":
                        PDB_modified.write(line[0:17]+"ASH"+line[20:])
                        if ASH_list == False:
                            ASH_list = list()
                        if line[22:26].strip() not in ASH_list:
                            ASH_list.append(line[22:26].strip())
                    # Conditional statement used when Glu residue is protonated in the carboxyl group
                    elif line[12:16] == " HE2" and line[17:20] == "GLU":
                        PDB_modified.write(line[0:17]+"GLH"+line[20:])
                        if GLH_list == False:
                            GLH_list = list()
                        if line[22:26].strip() not in GLH_list:
                            GLH_list.append(line[22:26].strip())
                    else:
                        if HXT_terminal == True:
                            print(line)
                        PDB_modified.write(line)
        else:
            pass

        i += 1

    PDB_modified.close()
    PDB_original.close()

    return water_i, Non_aminoacid_dict, HXT_terminal, Gly_to_other_residue, other_residue_to_Gly, ASH_list, GLH_list

def Printing_summary(PDB_filename, Output, water_i, Non_aminoacid_dict, HXT_terminal, Gly_to_other_residue, other_residue_to_Gly, ASH_list, GLH_list):
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
    HXT_terminal : boolean
                      The HXT atom name is present in the C-terminal residue of the protein chain/s (when True)
    Gly_to_other_residue : boolean / list
                      The list of new residues mutated have the HA2 atom name to the HA atom name due to the previous Gly residue (when True)
    ASH_list : boolean / list
                      The list of ASH residues (ASP with HD2) in the protein chain/s (when True)
    GLH_list : boolean / list
                      The list of GLH residues (GLU with HD2) in the protein chain/s (when True)     

    RETURNS
    -------
    PDB modified file
    """

    Screen_ticks = 10
    
    # Depending on if a output filename is specified on the command line or not, the input PDB file will be preserved or overwritten.
    if Output != "":
        print("\n{} has been succesfully processed and written to {}.pdb\n\n".format(PDB_filename, Output) + "-"*Screen_ticks + "SUMMARY" + "-"*Screen_ticks)
    else:
        os.system("rm {}" .format(PDB_filename))
        os.system("mv {}_modified.pdb {}" .format(PDB_filename[:-4], PDB_filename))
        print("\n{} has been succesfully processed and overwritten\n\n" .format(PDB_filename) + "-"*Screen_ticks + "SUMMARY" + "-"*Screen_ticks)
    if water_i == 0:
        print("\nNo water molecules were found in the PDB file")
    elif water_i == 1:
        print("\n{} water molecule was found in the PDB file".format(water_i))
    else:
        print("\n{} water molecules were found in the PDB file".format(water_i))
    if HXT_terminal:
        print("\nThe C-terminal residue contained the HXT atom name automatically put by Maestro/Schrödinger. Be aware! (Modified by the code)")
    if Gly_to_other_residue:
        print("\nThe following residues contained the typical HA2 atom name from a GLY residue but were another residue: {}. (Modified by the code)".format(Gly_to_other_residue))
    if other_residue_to_Gly:
        print("\nThe following residues contained the typical HA atom name from any residue, apart from GLY and they are a GLY residue: {}. (Modified by the code)".format(other_residue_to_Gly))
    if ASH_list:
        print("\nThe following residues are really ASH residues: {}. (Modified manually by the user)".format(ASH_list))
    if GLH_list:
        print("\nThe following residues are really GLH residues: {}. (Modified manually by the user)".format(GLH_list))
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
        water_i, Non_aminoacid_dict, HXT_terminal, Gly_to_other_residue, other_residue_to_Gly, ASH_list, GLH_list = PDB_processing(PDB_filename, Output, Residue_name)
        Printing_summary(PDB_filename, Output, water_i, Non_aminoacid_dict, HXT_terminal, Gly_to_other_residue, other_residue_to_Gly, ASH_list, GLH_list)

if __name__ == "__main__":
    """Call the main function"""
    main()
