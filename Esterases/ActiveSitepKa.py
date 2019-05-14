# -*- coding: utf-8 -*-


# Imports
import os,sys,re
import glob
import numpy as n

# Script information
__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"


def computepKa(PDB):
    """Take the PDB file and calculate the pKa of titrable residues using propka

    PARAMETERS
    ----------
    PDB : string
                       PDB file that wants to be added to the analysis

    RETURNS
    -------
    Results : dict of titrable residues with the calculated pKa

    pI_folded: The isoelectric point of the protein in the folded state

    pI_unfolded: The isoelectric point of the protein in the unfolded state
    """

    Results,index_pKa1,index_pKa2,pI_folded,pI_unfolded = {},0,0,0,0
    os.system("propka31 %s"%PDB)
    pKa_file = open("%s.pka -q"%PDB[:-4])

    for line in pKa_file:
        if "SUMMARY OF THIS PREDICTION" in line:
            index_pKa1=1
            continue
        if index_pKa1!=0:
            index_pKa2=index_pKa1
            index_pKa1=0
            continue
        if index_pKa2!=0:
            Results[line[3:6]+"_"+line[7:10]] = [int(line[7:10]),float(line[16:21])]
        if "N+" in line and index_pKa2!=0:
            Results[line[3:6]+"_"+line[7:10]] = [int(line[7:10]), float(line[16:21])]
            index_pKa2=0
        if "The pI is " in line:
            pI_folded, pI_unfolded = float(line[11:15]), float(line[30:34])

    return Results,pI_folded,pI_unfolded


def Neighbouratoms(PDB,Ser_residue,Results):
    """Take the atoms near the active site to store the pKa values

    PARAMETERS
    ----------
    PDB : string
                       PDB file that wants to be added to the analysis

    Ser_residue : int
                       Index or number referring to the catalytic Ser residue

    Results : dict
                       dict of titrable residues with the calculated pKa

    RETURNS
    -------
    Aux_results : dict of titrable residues with the calculated pKa in the active site
    """

    Aux_results = {}

    # Get the coordinates of the Ser residue to look for the neighbour titrable residues
    PDB_file = open(PDB, "rt")
    for line in PDB_file:
        if line[17:20]=="SER" and int(Ser_residue)==int(line[23:26]) and "OG" in line:
            x,y,z = float(line[30:38]),float(line[38:46]),float(line[46:54])

    # Get the neighbour residues and store them with the pKa value
    PDB_file = open(PDB, "rt")
    for line in PDB_file:
        if "TER" in line:
            pass
        else:
            x_aux, y_aux, z_aux = float(line[30:38]), float(line[38:46]), float(line[46:54])
            if n.sqrt((x-x_aux)**2+(y-y_aux)**2+(z-z_aux)**2)<=float(10):
                print(line[17:20]+"_"+line[23:26])
                if line[17:20]+"_"+line[23:26] in Results:
                    Aux_results[line[17:20]+"_"+line[23:26]] = Results[line[17:20]+"_"+line[23:26]]

    return Aux_results

def returnvalues(pI_folded,pI_unfolded,pKa_results):
    """
    Take the calculated values in the previous functions and return the pI values

    PARAMETERS
    ----------
    pI_folded : float
                       Calculated pI value in the folded state

    pI_unfolded : float
                       Calculated pI value in the unfolded state

    pKa_results : dict
                       Dict of all pKa values from the titrable residues around the active site

    RETURNS
    -------
    pI_folded : float
                       Calculated pI value in the folded state

    pI_unfolded : float
                       Calculated pI value in the unfolded state

    pI_active_site : float
                       Calculated pI value from the active site

    """

    values=[]
    for value in list(pKa_results.values()):
        values.append(value[1])
    pI = n.mean(values)

    return pI,pI_folded,pI_unfolded




def main(PDB,Ser_residue):
    """Main function

    It is called when this script is the main program called by the interpreter
    """

    # Compute pKa of titrable residues from the PDB file
    Results = computepKa(PDB)

    #Store the pKa of the residues around the active site
    pI=Neighbouratoms(PDB,Ser_residue,Results[0])

    returnvalues(Results[1],Results[2],pI)

if __name__ == "__main__":
    """Call the main function"""
    main(PDB,Ser_residue)
