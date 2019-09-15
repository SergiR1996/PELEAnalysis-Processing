# -*- coding: utf-8 -*-


# Imports 
import os,sys
import glob
import argparse as ap
import multiprocessing as mp # Module to parallelize job.

# Local imports
import PDBProcessor4PELE as PPP

# Script information
__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"

class Protein():

    def __init__(self,filename = "",pH = 7.0,JP = 0):

        self.__filename = filename
        self.__pH = pH
        self.__JP = JP

    @property
    def filename(self):
        return self.__filename

    def storePDBfilenames(PDBs_to_parse, parser):
	    """It identifies the PDB files to add to the ProteinPreparationSystem tool

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

    def parseArgs(self):

        """Parse arguments from command-line

        RETURNS
        -------
        PDBfiles : list
                  list of PDB files
        pH : float
                  specific value of the pH of the system
    	JP : int (0,1)
              Binary value to indicate if HETATM wants to be erased
        """

        parser = ap.ArgumentParser(description="Script used to prepare PDB files with prepwizard")
        required = parser.add_argument_group('required arguments')
        required.add_argument("-i","--input",required=True,metavar="STRING",
                                type=str,nargs="*",help="Location of input and output files")
        required.add_argument("-H","--pH",required=True,metavar="FLOAT",
                              type=float,nargs="*",help="pH of the system")
        required.add_argument("-JP","--justprotein",required=True,metavar="INTEGER",
                          type=int,nargs="*",help="Use the protein or all atoms",default=0)
        args = parser.parse_args()

        self.__filename = storePDBfilenames(args.input,parser)
        self.__pH = args.pH[0]
        self.__JP = args.justprotein

        return self.__filename,float(self.__pH),self.__JP

    def ProteinPreparationSystem(self, file = None):

        """ Take a PDB file and preparate the protein system using prepwizard of Schrodinger utilities.

        RETURNS
        -------
        Results : file
                  PDB file with the added hydrogens and the pH-sensitive groups properly modified according
                  to the pH of the system.
        """

        os.system("$SCHRODINGER/utilities/prepwizard %s %s -epik_pH %s -propka_pH %s" %(file,
            file,self.__pH,self.__pH))

        if self.__JP==0: # Preprocess as usual.
            PPP.PDB_processing(file)

        if self.__JP==1: # If just the protein is wanted, execute this part of the code.
            PDB_original, PDB_modified = open("%s" % file, "rt"), open("%s_modified.pdb" % file[:-4], "wt")
            for line in PDB_original:
                if line.find("ATOM")!=-1:
                    PDB_modified.write(line)
                elif line.find("TER")!=-1:
                    PDB_modified.write(line)
                    break
            PDB_modified.close()
            PDB_original.close()
            os.system("rm %s" % file)
            os.system("mv %s_modified.pdb %s" % (file[:-4], file))

    def preparesystem(self):

        """Main function

        It is called when this script is the main program called by the interpreter
        """

        # Get the protein in the PDB file and prepare it by adding H and performing the good protonations using parallelization.
        pool = mp.Pool(6)
        pool.map(self.ProteinPreparationSystem, self.__filename)
        pool.terminate()

if __name__ == "__main__":
    """Call the main function"""
	PDB_file = Protein()
	PDB_file.parseArgs()
	PDB_file.preparesystem()

