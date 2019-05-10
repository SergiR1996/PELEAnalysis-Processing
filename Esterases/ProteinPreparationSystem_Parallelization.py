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
    JP : int (0,1)
              Binary value to indicate if HETATM wants to be erased
    """

    parser = ap.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-H","--pH",required=True,metavar="FLOAT",
                          type=float,nargs="*",help="pH of the system")
    optional.add_argument("-JP","--justprotein",required=True,metavar="INTEGER",
                          type=int,nargs="*",help="Use the protein or all atoms",default=0)
    parser._action_groups.append(optional)
    args = parser.parse_args()

    pH = args.pH
    JP = args.justprotein

    return float(pH[0]),JP[0]

def ProteinPreparationSystem(PDBfilename):
    """ Take a PDB file and preparate the protein system using prepwizard of Schrodinger utilities.

    RETURNS
    -------
    Results : file
              PDB file with the added hydrogens and the pH-sensitive groups properly modified according
              to the pH of the system.
    """

    pH, JP = parseArgs() # Parse arguments inside the function to enable the call of parallelization.
    
    os.system("$SCHRODINGER/utilities/prepwizard %s %s -epik_pH %s -propka_pH %s" %(PDBfilename,PDBfilename,pH,pH))

    if JP==0: # Preprocess as usual.
        PPP.PDB_processing(PDBfilename[:-4])

    if JP==1: # If just the protein is wanted, execute this part of the code.
        PDB_original, PDB_modified = open("%s" % PDBfilename, "rt"), open("%s_modified" % PDBfilename, "wt")
        for line in PDB_original:
            if line.find("ATOM")!=-1:
                PDB_modified.write(line)
            elif line.find("TER")!=-1:
                PDB_modified.write(line)
                break
        PDB_modified.close()
        PDB_original.close()
        os.system("rm %s" % PDBfilename)
        os.system("mv %s_modified %s" % (PDBfilename, PDBfilename))



def main():
    """Main function

    It is called when this script is the main program called by the interpreter
    """

    # Store PDB filenames on a list
    PDB_List = glob.glob("*.pdb")

    # Get the protein in the PDB file and prepare it by adding H and performing the good protonations using parallelization.
    pool = mp.Pool(6) # In this part, the number of processors/cores is set.
    pool.map(ProteinPreparationSystem, PDB_List) # The function is called by the list of PDB files iteratively using parallelization.
    pool.terminate() # The parallelization and use of cores is killed once the job has finished.




if __name__ == "__main__":
    """Call the main function"""
    main()

