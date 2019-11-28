# -*- coding: utf-8 -*-

import re,sys # Imports the required modules to build the functions for the python tool. The
# re module is used for the checking functions and the adaptor-removal operation. The sys module
# is used to work with the command line arguments.

Blosum = {('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ('V', 'T'): 0, ('Q', 'Q'): 5, ('N', 'A'): -2, ('Z', 'Y'): -2, ('W', 'R'): -3, ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 8, ('S', 'H'): -1, ('H', 'D'): -1, ('L', 'N'): -3, ('W', 'A'): -3, ('Y', 'M'): -1, ('G', 'R'): -2, ('Y', 'I'): -1, ('Y', 'E'): -2, ('B', 'Y'): -3, ('Y', 'A'): -2, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 7, ('G', 'N'): 0, ('E', 'C'): -4, ('Y', 'Q'): -1, ('Z', 'Z'): 4, ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -2, ('T', 'N'): 0, ('P', 'P'): 7, ('V', 'I'): 3, ('V', 'S'): -2, ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -2, ('V', 'Q'): -2, ('K', 'K'): 5, ('P', 'D'): -1, ('I', 'H'): -3, ('I', 'D'): -3, ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2, ('P', 'H'): -2, ('F', 'Q'): -3, ('Z', 'G'): -2, ('X', 'L'): -1, ('T', 'M'): -1, ('Z', 'C'): -3, ('X', 'H'): -1, ('D', 'R'): -2, ('B', 'W'): -4, ('X', 'D'): -1, ('Z', 'K'): 1, ('F', 'A'): -2, ('Z', 'W'): -3, ('F', 'E'): -3, ('D', 'N'): 1, ('B', 'K'): 0, ('X', 'X'): -1, ('F', 'I'): 0, ('B', 'G'): -1, ('X', 'T'): 0, ('F', 'M'): 0, ('B', 'C'): -3, ('Z', 'I'): -3, ('Z', 'V'): -2, ('S', 'S'): 4, ('L', 'Q'): -2, ('W', 'E'): -3, ('Q', 'R'): 1, ('N', 'N'): 6, ('W', 'M'): -1, ('Q', 'C'): -3, ('W', 'I'): -3, ('S', 'C'): -1, ('L', 'A'): -1, ('S', 'G'): 0, ('L', 'E'): -3, ('W', 'Q'): -2, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 0, ('N', 'R'): 0, ('H', 'C'): -3, ('Y', 'N'): -2, ('G', 'Q'): -2, ('Y', 'F'): 3, ('C', 'A'): 0, ('V', 'L'): 1, ('G', 'E'): -2, ('G', 'A'): 0, ('K', 'R'): 2, ('E', 'D'): 2, ('Y', 'R'): -2, ('M', 'Q'): 0, ('T', 'I'): -1, ('C', 'D'): -3, ('V', 'F'): -1, ('T', 'A'): 0, ('T', 'P'): -1, ('B', 'P'): -2, ('T', 'E'): -1, ('V', 'N'): -3, ('P', 'G'): -2, ('M', 'A'): -1, ('K', 'H'): -1, ('V', 'R'): -3, ('P', 'C'): -3, ('M', 'E'): -2, ('K', 'L'): -2, ('V', 'V'): 4, ('M', 'I'): 1, ('T', 'Q'): -1, ('I', 'G'): -4, ('P', 'K'): -1, ('M', 'M'): 5, ('K', 'D'): -1, ('I', 'C'): -1, ('Z', 'D'): 1, ('F', 'R'): -3, ('X', 'K'): -1, ('Q', 'D'): 0, ('X', 'G'): -1, ('Z', 'L'): -3, ('X', 'C'): -2, ('Z', 'H'): 0, ('B', 'L'): -4, ('B', 'H'): 0, ('F', 'F'): 6, ('X', 'W'): -2, ('B', 'D'): 4, ('D', 'A'): -2, ('S', 'L'): -2, ('X', 'S'): 0, ('F', 'N'): -3, ('S', 'R'): -1, ('W', 'D'): -4, ('V', 'Y'): -1, ('W', 'L'): -2, ('H', 'R'): 0, ('W', 'H'): -2, ('H', 'N'): 1, ('W', 'T'): -2, ('T', 'T'): 5, ('S', 'F'): -2, ('W', 'P'): -4, ('L', 'D'): -4, ('B', 'I'): -3, ('L', 'H'): -3, ('S', 'N'): 1, ('B', 'T'): -1, ('L', 'L'): 4, ('Y', 'K'): -2, ('E', 'Q'): 2, ('Y', 'G'): -3, ('Z', 'S'): 0, ('Y', 'C'): -2, ('G', 'D'): -1, ('B', 'V'): -3, ('E', 'A'): -1, ('Y', 'W'): 2, ('E', 'E'): 5, ('Y', 'S'): -2, ('C', 'N'): -3, ('V', 'C'): -1, ('T', 'H'): -2, ('P', 'R'): -2, ('V', 'G'): -3, ('T', 'L'): -1, ('V', 'K'): -2, ('K', 'Q'): 1, ('R', 'A'): -1, ('I', 'R'): -3, ('T', 'D'): -1, ('P', 'F'): -4, ('I', 'N'): -3, ('K', 'I'): -3, ('M', 'D'): -3, ('V', 'W'): -3, ('W', 'W'): 11, ('M', 'H'): -2, ('P', 'N'): -2, ('K', 'A'): -1, ('M', 'L'): 2, ('K', 'E'): 1, ('Z', 'E'): 4, ('X', 'N'): -1, ('Z', 'A'): -1, ('Z', 'M'): -1, ('X', 'F'): -1, ('K', 'C'): -3, ('B', 'Q'): 0, ('X', 'B'): -1, ('B', 'M'): -3, ('F', 'C'): -2, ('Z', 'Q'): 3, ('X', 'Z'): -1, ('F', 'G'): -3, ('B', 'E'): 1, ('X', 'V'): -1, ('F', 'K'): -3, ('B', 'A'): -2, ('X', 'R'): -1, ('D', 'D'): 6, ('W', 'G'): -2, ('Z', 'F'): -3, ('S', 'Q'): 0, ('W', 'C'): -2, ('W', 'K'): -3, ('H', 'Q'): 0, ('L', 'C'): -1, ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -4, ('W', 'S'): -3, ('S', 'E'): 0, ('H', 'E'): 0, ('S', 'I'): -2, ('H', 'A'): -2, ('S', 'M'): -1, ('Y', 'L'): -1, ('Y', 'H'): 2, ('Y', 'D'): -3, ('E', 'R'): 0, ('X', 'P'): -2, ('G', 'G'): 6, ('G', 'C'): -3, ('E', 'N'): 0, ('Y', 'T'): -2, ('Y', 'P'): -3, ('T', 'K'): -1, ('A', 'A'): 4, ('P', 'Q'): -1, ('T', 'C'): -1, ('V', 'H'): -3, ('T', 'G'): -2, ('I', 'Q'): -3, ('Z', 'T'): -1, ('C', 'R'): -3, ('V', 'P'): -2, ('P', 'E'): -1, ('M', 'C'): -1, ('K', 'N'): 0, ('I', 'I'): 4, ('P', 'A'): -1, ('M', 'G'): -3, ('T', 'S'): 1, ('I', 'E'): -3, ('P', 'M'): -2, ('M', 'K'): -1, ('I', 'A'): -1, ('P', 'I'): -3, ('R', 'R'): 5, ('X', 'M'): -1, ('L', 'I'): 2, ('X', 'I'): -1, ('Z', 'B'): 1, ('X', 'E'): -1, ('Z', 'N'): 0, ('X', 'A'): 0, ('B', 'R'): -1, ('B', 'N'): 3, ('F', 'D'): -3, ('X', 'Y'): -1, ('Z', 'R'): 0, ('F', 'H'): -1, ('B', 'F'): -3, ('F', 'L'): 0, ('X', 'Q'): -1, ('B', 'B'): 4}

# Script information
__author__ = "Sergi Rodà"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Sergi Rodà"
__email__ = "sergi.rodallordes@bsc.es"

def open_fasta_file(fasta):
    """This function takes a simple or multi-FASTA file and returns a list
    containing the IDs and the sequences of each DNA chain in the file.

    PARAMETERS
    ----------
    fasta : string
        filename of a FASTA file

    RETURNS
    -------
    File : tuple of strings
        the tuple contains the ID of the read/sequence and the DNA sequence
    """

    ID,Sequence,File=[],[],[]

    for line in fasta: # This for loop iterates for all the lines in the FASTA file. It checks
        # if the line is the ID of sequence [">" symbol contained] or it is the sequence itself.
        if ">" in line:
            ID.append(line)
            Sequence.append("")
        else:
            line=line.strip("\n")
            Sequence[-1]=Sequence[-1]+line

    for i in range(len(ID)): File.append([ID[i],Sequence[i]])

    return File

def open_fastq_file(fastq):
    """This function takes a FASTQ file and returns a list containing the IDs,the sequences,
    the "+" symbols and the quality of each DNA chain in the file.

    PARAMETERS
    ----------
    fastq : string
             filename of a FASTQ file

    RETURNS
    -------
    File : tuple of strings
                  the tuple contains the ID of the read/sequence and the
                  DNA sequence (and the plus and the quality)
    """

    ID, Sequence,plus,quality,File,Counter = [], [], [], [], [], 1

    for line in fastq: # This for loop iterates for all the lines in the FASTQ file. It goes in a inner 4
        # conditional statements being the first line the ID of the sequence, the second one, the sequence.
        # The third and fourth one representing the "+" symbol and the quality of the sequence.
        if Counter%4==1: ID.append(line)
        elif Counter%4==2:
            line = line.strip("\n")
            Sequence.append(line)
        elif Counter%4==3: plus.append(line)
        elif Counter%4==0:
            line = line.strip("\n")
            quality.append(line)
        Counter+=1

    for i in range(len(ID)): File.append([ID[i],Sequence[i],plus[i],quality[i]])

    return File

def reverse_DNA(DNA_string):
    """This function takes a DNA string and returns the reverse-complement sequence. It uses the
    Nucleotides dictionary to change the nucleotides with and iterative for loop.

    PARAMETERS
    ----------
    DNA_string : string
             DNA sequence of the FASTA/Q file

    RETURNS
    -------
    The reverse-complement of the DNA_string.
    """

    Nucleotides={"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

    return "".join(Nucleotides[DNA_string[i]] for i in range(len(DNA_string)-1,-1,-1))

def Count_bases_and_percentages(Fasta_q):
    """This function takes whether the list of the open_fasta_file or open_fastq_file and 
       returns the number of bases processed and the relative abundance of each nucleotide.

    PARAMETERS
    ----------
    Fasta_q : Tuple of strings
             Tuple of properties of the different lines in the FASTA/Q file

    RETURNS
    -------
    The total number of bases in the FASTA/Q file and the relative abundance of each nucleotide.
    """

    Bases,Nucleotides=0,{"A":0,"C":0,"G":0,"T":0,"N":0}
    for i in range(len(Fasta_q)):
        Bases += float(len(Fasta_q[i][1]))
        for Nucleotide in Nucleotides: Nucleotides[Nucleotide] += float(Fasta_q[i][1].count(Nucleotide))

    return "{:,}".format(int(Bases)).replace(",","."),\
           tuple([int((Nucleotides[Nucleotide] / Bases) * 100) for Nucleotide in Nucleotides]) #The "{:,}".format(int(Bases)).replace(",",".")     part is used to return the number with . when the number has
    # more than 3 digits. What it basically does, it takes the number and replaces the implicit , with a .
    # (due to the english format of numbers, which they use ,). The relative abundance of nucleotides must be
    # changed to a tuple to later printing it with the % format options.

def edit_modified_distance_dp(pattern,text):
    """This function takes 'pattern' and 'text' nucleotide sequences and compares them forming the DP matrix,
    which contains the path of the optimal alignment. It uses a modified version of the edit distance approach,
    which contemplates the presence of transpositions of 2 bases (its score is also 1).

    PARAMETERS
    ----------
    pattern : string
             DNA sequence of the first input file
    text : string
             DNA sequence of the second input file

    RETURNS
    -------
    dp_matrix : list of lists (matrix)
         Dynamic Programming matrix of the alignment
    """

    dp_matrix = [[0 for _ in range(len(text)+1)] for _ in range(len(pattern)+1)]
    for v in range(len(pattern)+1): dp_matrix[v][0] = v
    for h in range(len(text)+1): dp_matrix[0][h] = h
    # Compute DP Matrix
    for h in range(1,len(text)+1):
        for v in range(1,len(pattern)+1):
            dp_matrix[v][h]=min(dp_matrix[v - 1][h - 1] + (0 if pattern[v - 1] == text[h - 1] else 1),
            dp_matrix[v][h - 1] + 1,dp_matrix[v - 1][h] + 1)
            if v > 1 and h > 1 and pattern[v-1] == text[h - 2] and pattern[v - 2] == text[h-1]: # The tranposition condition.
                dp_matrix[v][h] = min(dp_matrix[v][h], dp_matrix[v - 2][h - 2] + 1)
    
    return dp_matrix

def backtrace_matrix(pattern,text,dp_matrix):
    """This function takes the calculated DP matrix and computes back the CIGAR using also
     the sequences. The transposition (T) operation is added into the possible operations.
    
    PARAMETERS
    ----------
    pattern : string
             DNA sequence of the first input file
    text : string
             DNA sequence of the second input file
    dp_matrix : list of lists (matrix)
             Dynamic Programming matrix of the alignment

    RETURNS
    -------
    cigar : list of strings
         CIGAR with the instructions of the alignment (M,X,D,I,T)
    """
    
    v,h,cigar= len(pattern),len(text),[]
    while v>0 and h>0:
        if dp_matrix[v][h] == dp_matrix[v-1][h] +1:
            v -= 1;cigar.insert(0,"D")
        elif dp_matrix[v][h] == dp_matrix[v][h-1] +1:
            h -= 1;cigar.insert(0,"I")
        elif v > 1 and h > 1 and pattern[v - 1] == text[h - 2] and pattern[v - 2] == text[h - 1]:
            v -= 1;h -= 1;cigar.insert(0,"T")
        else:
            v -= 1;h -= 1
            if pattern[v] == text[h]: cigar.insert(0,"M")
            else: cigar.insert(0,"X")
    if v>0:
        for _ in range(v): cigar.insert(0,"D")
    if h>0:
        for _ in range(h): cigar.insert(0,"I")

    return cigar

def Align_DNA(Output_filename,faq_1,faq_2):
    """
    Performs the alignment of DNA sequences

    PARAMETERS
    ----------
    Output_filename : string
             Filename of the alignment file
    faq_1 : tuple of strings
             DNA sequence of the first input file
    faq_2 : tuple of strings
             DNA sequence of the second input file

    RETURNS
    -------
    Alingment stored in the Output_filename
    """
    alignment,DP,BK=open(Output_filename,"wt"),[],[]
    for i in range(len(faq_1)):
        DP.append(edit_modified_distance_dp(faq_1[i][1],faq_2[i][1]))
        BK.append(backtrace_matrix(faq_1[i][1],faq_2[i][1],DP[i]))
        (pattern_txt,j),operation_txt,(text_txt,k) = ("",0),[],("",0)
        for op in BK[i]:
            if op == "M":
                pattern_txt += faq_1[i][1][j];j += 1
                operation_txt.append("|")
                text_txt += faq_2[i][1][k];k += 1
            elif op == "X":
                pattern_txt += faq_1[i][1][j];j += 1
                operation_txt.append(" ")
                text_txt += faq_2[i][1][k];k += 1
            elif op=="T":
                pattern_txt += faq_1[i][1][j];j += 1
                operation_txt[-1]="*";operation_txt.append("*")
                text_txt += faq_2[i][1][k];k += 1
            elif op == "I":
                pattern_txt += "-"
                operation_txt.append(" ")
                text_txt += faq_2[i][1][k];k += 1
            elif op == "D":
                pattern_txt += faq_1[i][1][j];j += 1
                operation_txt.append(" ")
                text_txt += "-"
        alignment.write("%s (Identity: %s%%, Mismatches: %s%%, Gaps: %s%%): \n" %(i,float(100*(BK[i].count("M")+BK[i].count("T")))/len(BK[i]),
            float(100*BK[i].count("X"))/len(BK[i]),
            float(100*(BK[i].count("I")+BK[i].count("D")))/len(BK[i])))
        alignment.write(pattern_txt+"\n");alignment.write("".join(operation_txt)+"\n"),alignment.write(text_txt+"\n")

def distance_dp_protein(pattern,text):
    """This function takes 'pattern' and 'text' protein sequences and compares them forming
    the dynamic programming (DP) matrix, which contains the path of the optimal alignment.
    It uses the same scoring methodology as the previous function.

    PARAMETERS
    ----------
    pattern : string
             Protein sequence of the first input file
    text : string
             Protein sequence of the second input file

    RETURNS
    -------
    dp_matrix : list of lists (matrix)
         Dynamic Programming matrix of the alignment
    """

    dp_matrix = [[0 for _ in range(len(text)+1)] for _ in range(len(pattern)+1)]
    for v in range(len(pattern)+1): dp_matrix[v][0] = -(2*v) # The first column is defined by the deletion penalty.
    for h in range(len(text)+1): dp_matrix[0][h] = -(4*h) # The first row is defined by the insertion penalty.
    # Compute DP Matrix
    for h in range(1,len(text)+1):
        for v in range(1,len(pattern)+1):
            try : score=dp_matrix[v-1][h-1] +Blosum[pattern[v-1], text[h-1]]
            except: score=dp_matrix[v - 1][h - 1] + Blosum[text[h - 1], pattern[v - 1]]
            dp_matrix[v][h] = max (score,dp_matrix[v][h-1] -4,dp_matrix[v-1][h] -2)

    return dp_matrix

def backtrace_protein_matrix(pattern,text,dp_matrix):
    """This function takes the calculated DP matrix and computes back the CIGAR using also
     the sequences. Once the CIGAR is obtained, we can output the alignment.
    
    PARAMETERS
    ----------
    pattern : string
             Protein sequence of the first input file
    text : string
             Protein sequence of the second input file
    dp_matrix : list of lists (matrix)
             Dynamic Programming matrix of the alignment

    RETURNS
    -------
    cigar : list of strings
         CIGAR with the instructions of the alignment (M,X,Y,Z,D,I)
    """

    v,h,cigar = len(pattern),len(text),[]
    while v>0 and h>0: # While loop to iterate through all the elements of the DP matrix.
        if dp_matrix[v][h] == dp_matrix[v-1][h] -2:
            v -= 1;cigar.insert(0,"D")
        elif dp_matrix[v][h] == dp_matrix[v][h-1] -4:
            h -= 1;cigar.insert(0,"I")
        else:
            v -= 1;h -= 1
            try: score = Blosum[pattern[v],text[h]]
            except: score = Blosum[text[h],pattern[v]]
            if pattern[v] == text[h]: cigar.insert(0,"M")
            elif score>0: cigar.insert(0,"X") # This differentiation in the mismatches is used to see the importance of
            # the changed amino acid in comparison with the pattern/reference.
            elif score==0: cigar.insert(0,"Y")
            else: cigar.insert(0,"Z")
    if v>0: # If the pattern sequence is bigger than the text sequence.
        for _ in range(v): cigar.insert(0,"D")
    if h>0:  # If the text sequence is bigger than the pattern sequence.
        for _ in range(h): cigar.insert(0,"I")
    
    return cigar

def Align_Protein(Output_filename,faq_1,faq_2):
    """
    Performs the alignment of DNA sequences

    PARAMETERS
    ----------
    Output_filename : string
             Filename of the alignment file
    faq_1 : tuple of strings
             Protein sequence of the first input file
    faq_2 : tuple of strings
             Protein sequence of the second input file

    RETURNS
    -------
    Alingment stored in the Output_filename
    """
    alignment,DP,BK=open(Output_filename,"wt"),[],[]
    for i in range(len(faq_1)):
        DP.append(distance_dp_protein(faq_1[i][1],faq_2[i][1]))
        BK.append(backtrace_protein_matrix(faq_1[i][1],faq_2[i][1],DP[i]))
        (pattern_txt,j),operation_txt,(text_txt,k) = ("",0),[],("",0)
        for op in BK[i]:
            if op == "M":
                pattern_txt += faq_1[i][1][j];j += 1
                operation_txt.append("|")
                text_txt += faq_2[i][1][k];k += 1
            elif op == "X":
                pattern_txt += faq_1[i][1][j];j += 1
                operation_txt.append(":")
                text_txt += faq_2[i][1][k];k += 1
            elif op=="Y":
                pattern_txt += faq_1[i][1][j];j += 1
                operation_txt.append(".")
                text_txt += faq_2[i][1][k];k += 1
            elif op == "Z":
                pattern_txt += faq_1[i][1][j];j += 1
                operation_txt.append(" ")
                text_txt += faq_2[i][1][k];k += 1
            elif op == "D":
                pattern_txt += faq_1[i][1][j];j += 1
                operation_txt.append(" ")
                text_txt += "-"
            elif op == "I":
                pattern_txt += "-"
                operation_txt.append(" ")
                text_txt += faq_2[i][1][k];k += 1
        alignment.write("%s (Identity: %s%%, Similarity: %s%%, Gaps: %s%%): \n" %(i,float(100*BK[i].count("M"))/len(BK[i]),
            float(100*(BK[i].count("M")+BK[i].count("X")+BK[i].count("Y")))/len(BK[i]),
            float(100*(BK[i].count("I")+BK[i].count("D")))/len(BK[i])))
        alignment.write(pattern_txt+"\n");alignment.write("".join(operation_txt)+"\n"),alignment.write(text_txt+"\n")

class FASTA_Q():
   
    def __init__(self,Input_filename,Output_filename,trim_right=0,trim_left=0,Adaptor="",Input_filename2=""):
        self.__Input_filename = Input_filename
        self.__Output_filename = Output_filename
        self.__trim_right = trim_right
        self.__trim_left = trim_left
        self.__Adaptor = Adaptor
        self.__Input_filename2 = Input_filename2

    def ReverseComplementFasta(self):
        """This function opens the Input_filename (FASTA/Q) and returns an output with the Output_filename and with
        the reverse-complement sequences. A Summary of the reads and bases processed is also implemented."""
        f = open(self.__Input_filename, "rt") # The input file is opened, whether it is a FASTA or a FASTQ file.
        if ".fasta" in self.__Input_filename:Fasta_q,result=open_fasta_file(f),open(self.__Output_filename,"wt") # If FASTA file, use the open_fasta_file function and overwrite the output file.
        else: Fasta_q,result=open_fastq_file(f),open(self.__Output_filename,"wt") # If a FASTQ file, use the open_fastq_file.
        for i in range(len(Fasta_q)):
            result.write(Fasta_q[i][0]),result.write(reverse_DNA(Fasta_q[i][1])+"\n")
            if Fasta_q[i][0][0] == "@": # If a FASTQ file add the plus and quality lines.
                result.write(Fasta_q[i][2]),result.write(Fasta_q[i][3][::-1] + "\n")
        f.close(),result.close() #Once the output file is modified, a summary is created with the proper format in the assignment.
        Results=Count_bases_and_percentages(Fasta_q) # It counts the bases and the relative abundance of nucleotides in
            # the sequences of the input file.
        print ("File '%s' has been successfully reversed-complemented ('%s')" %(self.__Input_filename,self.__Output_filename))
        print ("Summary:\n"+(len(Results[0])-len(str(len(Fasta_q))))*" "+"%s reads processed"%("{:,}".format(len(Fasta_q)).replace(",",".")))
        print ("%s" % Results[0] + " bases processed (%s%% A, %s%% C, %s%% G, %s%% T, %s%% N)" % Results[1])
        # The spaces between numbers is defined by the number of digits and their correlation
        # (for example, if bases contains 5 digits and the reads contain 1 digit, 4 spaces must be entered in the read part).

    def TrimFasta(self):
        """This function opens the Input_filename (FASTA/Q) and returns an output with the Output_filename and with
        the sequences trimmed to the left and the right in a number of bases depending on the trim_left and trim_right
        parameters. A Summary of the reads, bases and trimmed bases processed is also implemented."""
        f,Trim = open(self.__Input_filename, "rt"),[]  # The input and output files are opened, whether they are a FASTA or a FASTQ file.
        if ".fasta" in self.__Input_filename:Fasta_q,result=open_fasta_file(f),open(self.__Output_filename,"wt") # If FASTA file, use the open_fasta_file function and overwrite the output file.
        else: Fasta_q,result=open_fastq_file(f),open(self.__Output_filename,"wt") # If a FASTQ file, use the open_fastq_file.
        for i in range(len(Fasta_q)): # It trims the sequences in the FASTA/FASTQ file.
            if self.__trim_left>=len(Fasta_q[i][1]) or self.__trim_right>=len(Fasta_q[i][1]) or (self.__trim_left+self.__trim_right)>=len(Fasta_q[i][1]):
                result.write("") # If the trim parameters are bigger than the input sequence, then it is skipped in the output file.
            else:
                result.write(Fasta_q[i][0])
                result.write(Fasta_q[i][1][int(self.__trim_left):len(Fasta_q[i][1])-int(self.__trim_right)]+"\n")
                Trim.append([Fasta_q[i][0],(Fasta_q[i][1][0:int(self.__trim_left)]+
                                          Fasta_q[i][1][len(Fasta_q[i][1])-int(self.__trim_right):len(Fasta_q[i][1])])])
                if Fasta_q[i][0][0] == "@":  # If a FASTQ file add the plus and the quality trimmed lines.
                    result.write(Fasta_q[i][2])
                    result.write(Fasta_q[i][3][int(self.__trim_left):len(Fasta_q[i][1]) - int(self.__trim_right)] + "\n")
        f.close(),result.close()
        Results,Results_trimmed=Count_bases_and_percentages(Fasta_q),Count_bases_and_percentages(Trim)
        print ("File '%s' has been successfully hard-trimmed ('%s')" %(self.__Input_filename,self.__Output_filename))
        print ("Summary:\n"+(len(Results[0])-len(str(len(Fasta_q))))*" "+"%s reads processed"% ("{:,}".format(len(Fasta_q)).replace(",",".")))
        print ("%s" % Results[0] + " bases processed (%s%% A, %s%% C, %s%% G, %s%% T, %s%% N)" % Results[1])
        print ((len(Results[0])-len(Results_trimmed[0]))*" "+"%s"%Results_trimmed[0]+" bases trimmed   "
                "(%s%% A, %s%% C, %s%% G, %s%% T, %s%% N)"% Results_trimmed[1])

    def AdaptorRemovalFasta(self):
        """This function opens the input filename (FASTA/Q) and returns an output with the Output_filename and with
        the sequences trimmed from the left until the Adaptor last base if the Adaptor is found in the sequence.
        A Summary of the reads, bases processed and found adaptors is also implemented."""
        f,Number_of_adaptors = open(self.__Input_filename, "rt"),0 # The input file is opened, whether is a FASTA or a FASTQ file.
        if ".fasta" in self.__Input_filename:Fasta_q,result=open_fasta_file(f),open(self.__Output_filename,"wt") # If FASTA file, use the open_fasta_file function and overwrite the output file.
        else: Fasta_q,result=open_fastq_file(f),open(self.__Output_filename,"wt") # If a FASTQ file, use the open_fastq_file.
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
        f.close(),result.close()
        Results = Count_bases_and_percentages(Fasta_q)
        print ("File '%s' has been successfully processed ('%s')" %(self.__Input_filename,self.__Output_filename))
        print ("Summary:\n"+(len(Results[0])-len(str(len(Fasta_q))))*" "+"%s reads processed"
                   % ("{:,}".format(len(Fasta_q)).replace(",",".")))
        print("%s" % Results[0] + " bases processed (%s%% A, %s%% C, %s%% G, %s%% T, %s%% N)" % Results[1])
        print ((len(Results[0])-len(str(Number_of_adaptors)))*" "+"%s adaptors found" %Number_of_adaptors)

    def Alignment(self):
        """This functions take 2 sequences and their alignment CIGAR and
        prints the nucleotide alignment in a more visual format (it adds the T operation)."""

        f1,f2=open(self.__Input_filename, "rt"),open(self.__Input_filename2, "rt")
        if ".fasta" in self.__Input_filename: faq_1,faq_2=open_fasta_file(f1),open_fasta_file(f2)
        else: faq_1,faq_2=open_fastq_file(f1),open_fastq_file(f2)

        if len(re.findall(r'[DEFHIKLMNPQRSVWY]',faq_1[0][1])) != 0:
            Align_Protein(self.__Output_filename,faq_1,faq_2)
        else: 
            Align_DNA(self.__Output_filename,faq_1,faq_2)

        print ("File '%s' has been successfully aligned with '%s'" %(self.__Input_filename,self.__Input_filename2))
