# -*- coding: utf-8 -*-

# Global imports
import glob
import argparse as ap
import pandas as pd

# Local imports
from ActiveSitepKa import *
from PDBOpener import *
from Esterase import *

# Script information
__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"

def main():
    """
    Main function

    It is called when this script is the main program called by the interpreter
    """

    # Store PDB filenames on a list
    PDB_List, csv, _ = parseArgs()
    input_df = pd.read_csv(csv);esterase_names = list(input_df["esterase"])
    output_df = pd.DataFrame(columns=["esterase", "pI_folded", "pI_unfolded", "pI_active_site"])

    for idx, PDB in enumerate(PDB_List):

        try:
            if str(PDB.split("/")[-1].split(".")[0]) in esterase_names:

                esterase = Esterase(PDB)
                cat, ser, active_type = esterase._DetectActiveSite()
                P = pKa(PDB,ser[1]);P.computepI()
                output_df.loc[idx, "esterase"] = PDB.split("/")[-1].split(".")[0]
                output_df.loc[idx,"pI_folded"] = P.pI[0]
                output_df.loc[idx,"pI_unfolded"] = P.pI[1]
                output_df.loc[idx,"pI_active_site"] = round(P.pI[2],3)

        except Exception as e:
            print("the error is {}".format(e))
            continue

    output_df.set_index("esterase",inplace=True); output_df.sort_index(inplace=True)
    output_df.to_csv("pI.csv", sep=",")

if __name__ == "__main__":
	"""Call the main function"""
	main()


