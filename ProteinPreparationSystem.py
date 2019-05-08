# -*- coding: utf-8 -*-


# Imports 
import os,sys
import glob
import argparse as ap

# Local imports
import PDBProcessor4PELE as PPP

# Script information
__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"

def storePDBfilenames(PDBs_to_parse,parser):
    """It identifies the reports to add to the protein preparation system

    PARAMETERS
    ----------
    reports_to_parse : list of strings
                       all the PDB files that want to be added to the analysis

    RETURNS
    -------
    parsed_data : list of PDB filenames (strings)
    """

    reports = []

    for reports_list in PDBs_to_parse:
        PDB_found = glob.glob(reports_list)
        if len(PDB_found) == 0:
            print("Warning: path to report file \'" +
                  "{}".format(reports_list) + "\' not found.")
        for report in PDB_found:
            reports.append(report)

    if len(reports) == 0:
        print("Error: list of report files is empty.")
        parser.print_help()
        exit(1)

    return reports

def parseArgs():
    """Parse arguments from command-line

    RETURNS
    -------
    PDBfiles : list
              list of PDB files
    pH : float
              specific value of the pH of the system
    """

    parser = ap.ArgumentParser()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", required=True, metavar="FILE",
                          type=str, nargs='*', help="path to PDB files")
    required.add_argument("-H","--pH",required=True,metavar="FLOAT",
                          type=float,nargs="*",help="pH of the system")
    args = parser.parse_args()

    PDBfiles = storePDBfilenames(args.input,parser)
    pH = args.pH

    return PDBfiles,float(pH[0])

def ProteinPreparationSystem(PDBfilename,pH):
    """ Take a PDB file and preparate the protein system using prepwizard of Schrodinger utilities.

    RETURNS
    -------
    Results : file
              PDB file with the added hydrogens and the pH-sensitive groups properly modified according
              to the pH of the system.
    """

    os.system("$SCHRODINGER/utilities/prepwizard %s %s -epik_pH %s -propka_pH %s" %(PDBfilename,PDBfilename,pH,pH))
    PPP.PDB_processing(PDBfilename[:-4])




def main():
    """Main function

    It is called when this script is the main program called by the interpreter
    """

    # Store PDB filenames on a list
    PDB_List,pH = parseArgs()

    # Get the protein in the PDB file and prepare it by adding H and performing the good protonations.
    for PDB in PDB_List:
        ProteinPreparationSystem(PDB,pH)



if __name__ == "__main__":
    """Call the main function"""
    main()

