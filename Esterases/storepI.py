# -*- coding: utf-8 -*-

from ActiveSitepKa import *
from Esterase import *
import glob
import argparse as ap
import pandas as pd
# Script information
__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"

def storePDBfilenames(PDBs_to_parse,parser):
    """It identifies the PDB files to add to the protein preparation system

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
            print("Warning: path to report file \'" +
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
    PDBfiles : list
              list of PDB files
    """

    parser = ap.ArgumentParser(description='Script used to compute the pI values \
        of the an esterase and its active site')
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", required=True, metavar="FILE",
                          type=str, nargs='*', help="path to PDB files")
    args = parser.parse_args()

    PDBfiles = storePDBfilenames(args.input,parser)

    return PDBfiles

def main():
	"""Main function

	It is called when this script is the main program called by the interpreter
	"""

	# Store PDB filenames on a list

	PDB_List = parseArgs()
	esterase_names = [line.split(",")[0] for line in open("final_promiscuity.csv","rt")][1:]
	df = pd.DataFrame(columns=["esterase", "pi_folded", "pi_unfolded", "pi_active_site"])

	for idx, PDB in enumerate(PDB_List):
		try:
			if str(PDB.split("/")[-1].split(".")[0]) in esterase_names:
				
				esterase = Esterase(PDB)
				cat, ser = esterase._DetectActiveSite()
				P = pKa(PDB,ser[1])
				PO = P.computepI()
				df.loc[idx, "esterase"] = PDB.split("/")[-1].split(".")[0]
				df.loc[idx,"pi_folded"] = PO[0]
				df.loc[idx,"pi_unfolded"] = PO[1]
				df.loc[idx,"pi_active_site"] = round(PO[2],3)

		except Exception as e:
			print("the error is {}".format(e))
			continue

	df.set_index("esterase",inplace=True); df.sort_index(inplace=True)
	df.to_csv("pI.csv", sep=",")

if __name__ == "__main__":
	"""Call the main function"""
	main()


