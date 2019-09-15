# -*- coding: utf-8 -*-

import sys,time
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
import itertools
import argparse as ap
import multiprocessing as mp

# Script information
__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"


class MutateScore():

	def __init__(self):
		
		self.__ref_file, self.__filename, self.__ref_residues_indices, self.__Atom_types, self.__output = self.parseArgs()

	def parseArgs(self):
		"""
		Parse arguments from command-line
		"""

		parser = ap.ArgumentParser(description='Script used to compute the local RMSD od the specified \
			residues of a reference PDB file against all the residues of a target PDB file')
		optional = parser._action_groups.pop()
		required = parser.add_argument_group('required arguments')
		parser.add_argument("reference", metavar="FILE",type=str, help="path to reference PDB file")
		parser.add_argument("input", metavar="FILE",type=str, help="path to input PDB file")
		required.add_argument("-r","--residues",required=True,metavar="STRING",
								type=int,nargs='*',help="reference residue indices")
		optional.add_argument("-AT","--Atom_types",metavar="STRING",type=str,nargs='*',
			help="Atom types for the RMSD calculation",default=["CA","N","O"])
		optional.add_argument("-O","--output",metavar="STRING",type=str,
			help="Output filename",default="scores.txt")
		parser._action_groups.append(optional)
		args = parser.parse_args()

		self.__ref_file = args.reference
		self.__filename =args.input
		self.__ref_residues_indices = args.residues
		self.__Atom_types = args.Atom_types
		self.__output = args.output

		return self.__ref_file, self.__filename, self.__ref_residues_indices, self.__Atom_types, self.__output

	@property
	def filename(self):
		return self.__filename

	@property
	def ref_residues_indices(self):
		return self.__ref_residues_indices
	

	def GetCoordinates(self,file,ref_coord = False):
		"""
		This method returns the coordinates of a PDB file and
		returns them. The reference coordinates are stored or 
		not according to a boolean value.

		PARAMETERS
        ----------
        file : string
                PDB filename
        ref_coord: bool
        		Boolean value to indicate whether they are
        		reference residues or not

        RETURNS
        ------
        res_coord: list of floats
		"""

		PDB = open(file)
		coordinates,res_coord,aux_coord,counter = [],[],[],0
		for line in PDB:
			if (line.split()[0] == "ATOM") and (line.split()[2] in self.__Atom_types):
				if ref_coord:
					if int(line.split()[5]) in self.__ref_residues_indices:
						x = float(line.split()[6])
						y = float(line.split()[7])
						z = float(line.split()[8])
						# res_name = line.split()[2]+"_"+line.split()[3]+"_"+line.split()[5]
						# coordinates.append([x,y,z])
						res_coord.append([x,y,z])
						# counter +=1
						# if counter%3 == 0:
						# 	res_coord.append(aux_coord)
						# 	aux_coord = []
				else:
					x = float(line.split()[6])
					y = float(line.split()[7])
					z = float(line.split()[8])
					# res_name = line.split()[2]+"_"+line.split()[3]+"_"+line.split()[5]
					# coordinates.append([x,y,z])
					aux_coord.append([x,y,z])
					counter +=1
					if counter%len(self.__Atom_types) == 0:
						res_coord.append(aux_coord)
						aux_coord = []

		return res_coord

	def RMSD(self,ref_coordinates,coordinates):
		"""
		This method computes the root-mean-square deviation (RMSD) for
		the target coordinates against the reference coordinates.

		PARAMETERS
        ----------
        ref_coordinates : array of floats
                Coordinates of the reference atoms
        coordinates: array of floats
        		Coordinates of the target atoms

        RETURNS
        ------
        rmsd: float
        		Value of the RMSD
		"""

		D,N,RMS = len(ref_coordinates[0]),len(ref_coordinates),0.0
		for v, w in zip(ref_coordinates, coordinates):
				RMS += sum([(v[i] - w[i])**2.0 for i in range(D)])
		rmsd = np.sqrt((RMS/N))

		return rmsd

	def Translate(self,coordinates):
		"""
		This method computes the centroid of some coordinates and
		substracts it from them, translating to the center of 
		coordinates.

		PARAMETERS
        ----------
        coordinates: array of floats
        		Coordinates of the atoms

        RETURNS
        ------
        centered_coordinates: array of floats
        		Centered coordinates of the atoms
		"""

		centered_coordinates = []
		C = np.mean(coordinates,axis=0)
		centered_coordinates = (coordinates-C)

		return centered_coordinates

	def Kabsch(self,coordinates,ref_coordinates):
		"""
		https://en.wikipedia.org/wiki/Kabsch_algorithm
		"""

		# Computation of the covariance matrix
		C = np.dot(np.transpose(coordinates), ref_coordinates)

		# Computation of the optimal rotation matrix
		V, S, W = np.linalg.svd(C)
		d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

		if d:
			S[-1] = -S[-1]
			V[:, -1] = -V[:, -1]

		# Create Rotation matrix U
		U = np.dot(V, W)

		return U

	def Rotate(self,coordinates,ref_coordinates):
		"""
		This method computes the optimal rotate matrix
		to rotate the coordinates according to the 
		reference coordinates.

		PARAMETERS
        ----------
        coordinates: array of floats
        		Coordinates of the target atoms
        ref_coordinates: array of floats
        		Coordinates of the reference atoms

        RETURNS
        ------
        coordinates: array of floats
        		Rotated coordinates of the atoms
		"""

		U = self.Kabsch(coordinates,ref_coordinates)

		# Rotate P
		coordinates = np.dot(coordinates,U)

		return coordinates

	def DecompressList(self,coordinates):
		"""
		This method decompress the coordinates saved 
		as a list of lists into a list

		PARAMETERS
        ----------
        coordinates: array of floats
        		Coordinates of the atoms

        RETURNS
        ------
        new_coordinates: array of floats
        		Rotated coordinates of the atoms	
		"""

		new_coordinates = []
		for sublist in coordinates:
			for item in sublist:
				new_coordinates.append(item)

		return new_coordinates

	def FindResidues(self,file,coordinates):
		"""
		This method finds  the residues that are contained in 
		the coordinates of the specified atoms.

		PARAMETERS
        ----------
        file: string
        		PDB filename
        coordinates: array of floats
        		Coordinates of the atoms

        RETURNS
        ------
        "".join(Residues): string
        		The name of the residues that contain the coordinates	
		"""

		def RoundFloat(number):
			"""
			This method appends up to 3 decimals to
			all coordinates in order to be converted
			to float in the FindResidues method.

			PARAMETERS
	        ----------
	        number: float
	        		coordinate in one of the axis

	        RETURNS
	        ------
	        A: string
	        		The string of the coodinates with 3 decimals	
			"""

			A = str(round(number,3))
			while len(A.split(".")[1]) < 3:
				A+="0"

			return A

		Residues,res_index,val = [],[],0
		PDB = open(file)
		lines = PDB.readlines();PDB.close()
		for i in range(len(coordinates)):
			for line in lines:
				if (line.split()[0] == "ATOM") and (line.split()[2] in self.__Atom_types):
					x,y,z = RoundFloat(coordinates[i][0][0]),RoundFloat(coordinates[i][0][1]),RoundFloat(coordinates[i][0][2])
					if (line.split()[6] == x) and (line.split()[7] == y) and (line.split()[8] == z):
						if line.split()[5] not in res_index:
							Residues.append(line.split()[3]+"_"+line.split()[5]+" ")

		return "".join(Residues)

	def ComputeScore(self,combination):
		"""
		This method computes the RMSD for all the permutations 
		in a combination of coordinates of some target 
		residues.

		PARAMETERS
        ----------
        combination: list of lists of floats
        		Coordinates of the atoms of the target residues

        RETURNS
        ------
        results: list of lists of a integer and a string
        		The RMSD and the residues of the combination	
		"""

		ref_coordinates = self.GetCoordinates(self.__ref_file, True)
		ref_coordinates = self.Translate(ref_coordinates)
		ref_coordinates = np.array(ref_coordinates)
		results = []

		for permutation in itertools.permutations(combination,len(self.__ref_residues_indices)):
			combination_cent = self.Translate(self.DecompressList(permutation))
			final_coordinates = self.Rotate(combination_cent,ref_coordinates)
			RMSD = self.RMSD(ref_coordinates,final_coordinates)
			Aux = [RMSD,self.FindResidues(self.__filename,permutation)]
			if Aux[0] < 1.0:
			 	print(Aux+"\n")
			results.append(Aux)

		return results

	def main(self):
		"""
		Main function
	
    	It is called when this script is the main program called by the interpreter
    	"""

		output = open(self.__output,"w")
		results = []

		start = time.time()
		print("RMSD of all combinations is starting to be computed \n")

		pool = mp.Pool(6)
		results.append(pool.map(
			A.ComputeScore,itertools.combinations(A.GetCoordinates(A.filename),len(A.ref_residues_indices))))
		pool.terminate()

		results = self.DecompressList(self.DecompressList((results)))
		
		results.sort(key= lambda x : x[0])

		for elem in results:
			output.write("\nRMSD: {}".format(elem[0])+" // Residues: {} \n".format(elem[1]))

		end = time.time()
		print("The main code needed {} seconds to compute all scores for all the combinations \n".format(end-start))


if __name__=="__main__":
	"""Call the main function"""
	A = MutateScore()
	A.main()