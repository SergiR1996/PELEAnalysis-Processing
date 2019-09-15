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
    reports : list
              list of PDB files
    """

    parser = ap.ArgumentParser(description='Script used to find the catalytic residues of an esterase in a PDB file')
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", required=True, metavar="FILE",
                          type=str, nargs='*', help="path to PDB sfiles")
    args = parser.parse_args()
    PDBfiles = storePDBfilenames(args.input,parser)

    return PDBfiles

class ActiveSite():

    def __init__(self,filename = ""):

        self.__filename = filename

    def PDBopener(self,residue_name):
        """ Take a PDB file, open it and store the position of Ser, His, Asp, Glu, Lys and Tyr residues

        RETURNS
        -------
        residues : dictionary
                  dictionary with the coordinates of the atoms in the catalytic residues
        """

        residues = {}
        PDB = open(self.__filename,'r')

        if residue_name == "SER": atom_name = "HG"
        elif residue_name == "HIS_1":
            residue_name = "HIS"
            atom_name = "NE2"
        elif residue_name == "HIS_2":
            residue_name = "HIS"
            atom_name = "HD1"
        elif residue_name == "ASP": atom_name = "OD2"
        elif residue_name == "GLU": atom_name = "OE2"
        elif residue_name == "LYS": atom_name = "NZ"
        elif residue_name == "TYR": atom_name = "OH"

        for line in PDB:
            if (line[0:6].strip() == "ATOM") and (line[17:20].strip() == residue_name) and (line[12:16].strip() == atom_name):
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    residue_number = int(line[22:26])
                    residues["{}_{}".format(residue_name, residue_number)] = [x,y,z]

        return residues

    def ActiveSiteFinder(self):
        """
        This method finds the catalytric triad of a given esterase (in a pdb) according to pre-fixed 
        thresholds for the H-bond distances in the catalytic residues
        """

        def ComputeDistance(atom1,atom2):
            """Computes the module or distance between the two points"""

            r = [atom2[0] - atom1[0], atom2[1] - atom1[1], atom2[2] - atom1[2]]
            return math.sqrt(r[0] ** 2 + r[1] ** 2 + r[2]**2)

        def ThresholdFinder(index1,index2,threshold,l):
            """Stores residues that statisfy the condition limited by the threshold"""

            for key1, value1 in result[index1].items():
                for key2, value2 in result[index2].items():
                    if (ComputeDistance(value1, value2) <= threshold):
                        l.append("{}-{}".format(key1,key2))

        def ActiveSiteReporter(l1,l2):
            """Takes the stored residues in the ThresholdFinder function and fins the ones that belong to the catalytic triad"""

            catalytic_residues = None 

            for elem1 in l1:
                nucleophile = elem1.split("-")[0]
                base1 = elem1.split("-")[-1]
                for elem2 in l2:
                    base2 = elem2.split("-")[0]
                    acid = elem2.split("-")[-1]
                    if base1==base2: catalytic_residues = "{}-{}-{}".format(nucleophile, base1, acid)

            return catalytic_residues

        residues = ["SER", "HIS_1", "HIS_2", "ASP", "GLU", "LYS", "TYR"]
        result = list(map(self.PDBopener,residues))
        catalytic_residues,active_site_type = None,None
        ser_his,his_asp,his_glu,ser_lys,lys_tyr = [],[],[],[],[]

        ThresholdFinder(0,1,4.5,ser_his)
        ThresholdFinder(2,3,4.5,his_asp)
        ThresholdFinder(2,4,4.5,his_glu)
        ThresholdFinder(0,5,4.0,ser_lys)
        ThresholdFinder(5,6,4.5,lys_tyr)

        cr1 = ActiveSiteReporter(ser_his,his_asp)
        cr2 = ActiveSiteReporter(ser_his,his_glu)
        cr3 = ActiveSiteReporter(ser_lys,lys_tyr)

        # Classify the type of active site
        if "ASP" in cr1 and cr1!=None:
            active_site_type = 0
            print("catalytic triad of {} is {}".format(self.__filename,cr1))
            return cr1, ["SER", int(cr1.split("-")[0].split("_")[-1])],active_site_type

        elif "GLU" in cr2:
            active_site_type = 1
            print("catalytic triad of {} is {}".format(self.__filename,cr2))
            return cr2, ["SER", int(cr2.split("-")[0].split("_")[-1])],active_site_type

        elif "LYS" in cr3:
            active_site_type = 2
            print("catalytic triad of {} is {}".format(self.__filename,cr3))
            return cr3, ["SER", int(cr3.split("-")[0].split("_")[-1])],active_site_type

def main():
    """Main function

    It is called when this script is the main program called by the interpreter
    """

    PDB_List = parseArgs()

    Dict_file=open("PDict_file","wt")

    for PDB in PDB_List:
        AS = ActiveSite(PDB)
        Results = AS.ActiveSiteFinder()
        if PDB==PDB_List[0]:
            Dict_file.write("Protein_dict = {"+ '"' + str(Results[0]) + '"' + ": " + str(Results[2]) + ", ")
        elif PDB==PDB_List[-1]:
            Dict_file.write('"' + str(Results[0]) + '"' + ": " + str(Results[2]) +"}")
        else:
            Dict_file.write('"' + str(Results[0]) + '"' + ": " + str(Results[2]) + ", ")


if __name__ == "__main__":
    """Call the main function"""
    main()

