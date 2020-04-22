# -*- coding: utf-8 -*-

# Global imports
import re,sys

# Local imports
from SequenceFunctionTools import *

# Global variables
AminoAcids = {"A":0,"C":0,"G":0,"T":0,"N":0,"S":0,"F":0,"Y":0,"W":0,"L":0,"D":0,"E":0,"K":0,"R":0,"I":0,"V":0,"P":0,"Q":0,"H":0,"M":0}

# Script information
__author__ = "Sergi Rodà"
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Sergi Rodà"
__email__ = "sergi.rodallordes@bsc.es"

class SequenceOperator():
    """
    Class that takes a FASTA/Q file containing a DNA/protein sequence and executes
    different user operations
    """
   
    def __init__(self,Input_filename,Output_filename,trim_right=0,trim_left=0,Adaptor="",Input_filename2=""):
        self.__Input_filename = Input_filename
        self.__Output_filename = Output_filename
        self.__trim_right = trim_right
        self.__trim_left = trim_left
        self.__Adaptor = Adaptor
        self.__Input_filename2 = Input_filename2

    def _OpenInputandOutputFiles(self):
        """
        This function handles the opening of input and output files

        RETURNS
        -------
        Fasta_q: _io.TextIOWrapper
                The opened input file in read mode
        result: _io.TextIOWrapper
                The opened output file in write mode
        """

        f = open(self.__Input_filename, "rt")
        if ".fasta" == self.__Input_filename[-6:] or ".fa" == self.__Input_filename[-3:] or ".fsa" == self.__Input_filename[-4:]:
            Fasta_q = open_fasta_file(f)
        else: 
            Fasta_q = open_fastq_file(f)
        result = open(self.__Output_filename,"wt")
        f.close()

        return Fasta_q, result

    def _PrintSummaryResults(self, Fasta_q, Trim = "", operation = "", Number_of_adaptors = ""):
        """
        This function takes the input file and it returns a summary of the bases/residues 
        processed and the ones that have been modified.

        OUTPUT
        ------
        The method outputs a summary of the processing executed in the sequence
        """

        Protein_nature = len(re.findall(r'[DEFHIKLMNPQRSVWY]',Fasta_q[0][1])) != 0
        if Protein_nature:
            Results = Count_bases_and_percentages(Fasta_q, AminoAcids)
            if operation == "trim":
                Results_trimmed = Count_bases_and_percentages(Trim, AminoAcids)
        else: 
            Results = Count_bases_and_percentages(Fasta_q)


        print ("File '%s' has been successfully processed ('%s')" %(self.__Input_filename,self.__Output_filename))
        print ("Summary:\n"+(len(Results[0])-len(str(len(Fasta_q))))*" "+"%s reads processed"%("{:,}".format(len(Fasta_q)).replace(",",".")))
        if Protein_nature:
            print ("%s" % Results[0] + " residues processed (%s%% A, %s%% C, %s%% G, %s%% T, %s%% N, %s%% S, %s%% F, %s%% Y, %s%% W, %s%% L, %s%% D, %s%% E, %s%% K, %s%% R, %s%% I, %s%% V, %s%% P, %s%% Q, %s%% H, %s%% M)" % Results[1])
            if operation == "trim":
                print ((len(Results[0])-len(Results_trimmed[0]))*" "+"%s"%Results_trimmed[0]+" residues trimmed   "
                        "(%s%% A, %s%% C, %s%% G, %s%% T, %s%% N, %s%% S, %s%% F, %s%% Y, %s%% W, %s%% L, %s%% D, %s%% E, %s%% K, %s%% R, %s%% I, %s%% V, %s%% P, %s%% Q, %s%% H, %s%% M)"% Results_trimmed[1])
        else:
            print ("%s" % Results[0] + " bases processed (%s%% A, %s%% C, %s%% G, %s%% T, %s%% N)" % Results[1])
            if operation == "trim":
                print ((len(Results[0])-len(Results_trimmed[0]))*" "+"%s"%Results_trimmed[0]+" bases trimmed   "
                        "(%s%% A, %s%% C, %s%% G, %s%% T, %s%% N)"% Results_trimmed[1])
        if type(Number_of_adaptors)!= str:
            print ((len(Results[0])-len(str(Number_of_adaptors)))*" "+"%s adaptors found" %Number_of_adaptors)

    def ReverseComplementFasta(self):
        """
        This function opens the Input_filename (FASTA/Q) and returns an output with the Output_filename and with
        the reverse-complement sequences. A Summary of the reads and bases processed is also implemented.

        OUTPUT
        ------
        This method returns the reverse complemented DNA sequence of the input FASTA/Q file.
        A summary of the processed bases is also returned.
        """

        Fasta_q, result = self._OpenInputandOutputFiles()
        for i in range(len(Fasta_q)):
            result.write(Fasta_q[i][0]),result.write(reverse_DNA(Fasta_q[i][1])+"\n")
            if Fasta_q[i][0][0] == "@":
                result.write(Fasta_q[i][2]),result.write(Fasta_q[i][3][::-1] + "\n")
        result.close()
        self._PrintSummaryResults(Fasta_q)


    def TrimFasta(self):
        """
        This function opens the Input_filename (FASTA/Q) and returns an output with the Output_filename and with
        the sequences trimmed to the left and the right in a number of bases depending on the trim_left and trim_right
        parameters. A Summary of the reads, bases and trimmed bases processed is also implemented.

        OUTPUT
        ------
        This method returns the trimmed DNA/Protein sequence of the input FASTA/Q file.
        A summary of the processed bases is also returned.
        """

        Fasta_q, result = self._OpenInputandOutputFiles(); Trim = []
        for i in range(len(Fasta_q)):
            if self.__trim_left>=len(Fasta_q[i][1]) or self.__trim_right>=len(Fasta_q[i][1]) or (self.__trim_left+self.__trim_right)>=len(Fasta_q[i][1]):
                result.write("")
            else:
                result.write(Fasta_q[i][0])
                result.write(Fasta_q[i][1][int(self.__trim_left):len(Fasta_q[i][1])-int(self.__trim_right)]+"\n")
                Trim.append([Fasta_q[i][0],(Fasta_q[i][1][0:int(self.__trim_left)]+
                                          Fasta_q[i][1][len(Fasta_q[i][1])-int(self.__trim_right):len(Fasta_q[i][1])])])
                if Fasta_q[i][0][0] == "@":
                    result.write(Fasta_q[i][2])
                    result.write(Fasta_q[i][3][int(self.__trim_left):len(Fasta_q[i][1]) - int(self.__trim_right)] + "\n")
        result.close()
        self._PrintSummaryResults(Fasta_q, Trim = Trim, operation = "trim")

    def AdaptorRemovalFasta(self):
        """
        This function opens the input filename (FASTA/Q) and returns an output with the Output_filename and with
        the sequences trimmed from the left until the Adaptor last base if the Adaptor is found in the sequence.
        A Summary of the reads, bases processed and found adaptors is also implemented.

        OUTPUT
        ------
        This method returns the DNA/Protein sequence of the input FASTA/Q file with the user adaptor removed.
        A summary of the processed bases is also returned.
        """

        Fasta_q, result = self._OpenInputandOutputFiles(); Number_of_adaptors = 0
        for i in range(len(Fasta_q)):
                if Fasta_q[i][1].find(self.__Adaptor.upper())!=-1:
                    if len(self.__Adaptor)==len(Fasta_q[i][1]):
                        result.write("")
                    else:
                        result.write(Fasta_q[i][0])
                        result.write(Fasta_q[i][1][Fasta_q[i][1].find(self.__Adaptor.upper())+len(self.__Adaptor):]+"\n")
                        if Fasta_q[i][0][0] == "@": result.write(Fasta_q[i][2]),result.write(Fasta_q[i][3][Fasta_q[i][1].find(self.__Adaptor.upper()) + len(self.__Adaptor):] + "\n")
                else:
                    result.write(Fasta_q[i][0])
                    result.write(Fasta_q[i][1] + "\n")
                    if Fasta_q[i][0][0] == "@": result.write(Fasta_q[i][2]),result.write(Fasta_q[i][3] + "\n")
                Number_of_adaptors+=len(re.findall(self.__Adaptor.upper(),Fasta_q[i][1]))
        result.close()
        self._PrintSummaryResults(Fasta_q, Number_of_adaptors = Number_of_adaptors)

    def Alignment(self):
        """
        This functions take 2 sequences and their alignment CIGAR and
        prints the nucleotide alignment in a more visual format 
        (it adds the T operation).

        OUTPUT
        ------
        The alignment between the two input sequences        
        """

        f1,f2=open(self.__Input_filename, "rt"),open(self.__Input_filename2, "rt")
        if ".fasta" == self.__Input_filename[-6:] or ".fa" == self.__Input_filename[-3:] or ".fsa" == self.__Input_filename[-4:]: 
            faq_1,faq_2=open_fasta_file(f1),open_fasta_file(f2)
        else: 
            faq_1,faq_2=open_fastq_file(f1),open_fastq_file(f2)

        if len(re.findall(r'[DEFHIKLMNPQRSVWY]',faq_1[0][1])) != 0:
            Align_Protein(self.__Output_filename,faq_1,faq_2)
        else: 
            Align_DNA(self.__Output_filename,faq_1,faq_2)

        print ("File '%s' has been successfully aligned with '%s'" %(self.__Input_filename,self.__Input_filename2))
