# -*- coding: utf-8 -*-


# Imports
import os,sys,re
import subprocess
import glob
import numpy as n

# Script information
__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"

class pKa:

    def __init__(self,PDB,Ser_residue):
        self.__PDB = PDB
        self.__Ser_residue = Ser_residue
        self.__Results = {}
        self.__pI_folded = 0
        self.__pI_unfolded = 0
        self.__pI_active_site = 0

    @property
    def PDB(self):

        return self.__PDB

    @property
    def Ser_residue(self):

        return self.__Ser_residue

    @property
    def pI(self):

        return self.__pI_folded,self.__pI_unfolded,self.__pI_active_site

    def propka(self):
        """Take the PDB file and calculate the pKa of titrable residues using propka

        PARAMETERS
        ----------
        PDB : string
                PDB file that wants to be added to the analysis

        OUTPUT
        ------
        Results : dict of titrable residues with the calculated pKa

        pI_folded: The isoelectric point of the protein in the folded state

        pI_unfolded: The isoelectric point of the protein in the unfolded state
        """

        index_pKa1,index_pKa2 = 0,0

        try:
            subprocess.call("propka31" +"%s" % self.__PDB + "-q")
            print("Computing pI values...")

        except:
            answer = input("propka is not installed. Do you want to install it? [Y/N]")
            if answer.lower() == "y" or answer.lower() == "yes":
                #subprocess.call(["git","clone","https://github.com/jensengroup/propka-3.1"])
                os.system("cd propka-3.1");os.system("python setup.py install --user")
                subprocess.call("propka31" +"%s" % self.__PDB + "-q")
            elif answer.lower() == "n" or answer.lower() == "No":
                exit()
            else:
                pass

        else:
            os.system("rm *.propka_input")
            pKa_file = open("%s.pka"%self.__PDB[self.__PDB.rindex("/")+1:-4])
            for line in pKa_file:
                if "SUMMARY OF THIS PREDICTION" in line:
                    index_pKa1=1
                    continue
                if index_pKa1!=0:
                    index_pKa2=index_pKa1
                    index_pKa1=0
                    continue
                if index_pKa2!=0:
                    self.__Results[line[3:6]+"_"+line[7:10]] = [int(line[7:10]),float(line[16:21])]
                if "N+" in line and index_pKa2!=0:
                    self.__Results[line[3:6]+"_"+line[7:10]] = [int(line[7:10]), float(line[16:21])]
                    index_pKa2=0
                if "The pI is " in line:
                    self.__pI_folded, self.__pI_unfolded = float(line[10:15]), float(line[29:34])
            os.system("rm *.pka")

    def Neighbouratoms(self):
        """Take the atoms near the active site to compute the pI around this area

        PARAMETERS
        ----------
        PDB : string
                           PDB file that wants to be added to the analysis

        Ser_residue : int
                           Index or number referring to the catalytic Ser residue

        Results : dict
                           dict of titrable residues with the calculated pKa

        OUTPUT
        ------
        pI_active_site : pI of the active site and surroundings (10 Å)
        """

        Aux_results,values = {},[]

        # Get the coordinates of the Ser residue to look for the neighbour titrable residues
        PDB_file = open(self.__PDB, "rt")
        for line in PDB_file:
            if line[17:20]=="SER" and int(self.__Ser_residue)==int(line[23:26]) and "OG" in line:
                x,y,z = float(line[30:38]),float(line[38:46]),float(line[46:54])

        # Get the neighbour residues and store them with the pKa value
        PDB_file = open(self.__PDB, "rt")
        for line in PDB_file:
            if "TER" in line:
                pass
            elif "ATOM" in line:
                x_aux, y_aux, z_aux = float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())
                if n.sqrt((x-x_aux)**2+(y-y_aux)**2+(z-z_aux)**2)<=float(10):
                    if line[17:20]+"_"+line[23:26] in self.__Results:
                        Aux_results[line[17:20]+"_"+line[23:26]] = self.__Results[line[17:20]+"_"+line[23:26]]
            else:
                pass
        self.__Results = Aux_results
        for value in list(Aux_results.values()):
            values.append(value[1])

        self.__pI_active_site = n.mean(values)

    def computepI(self):
        """It executes the methods of the class sequentially,
        returning the 3 computed values of pI"""

        self.propka()
        self.Neighbouratoms()
        self.pI()
