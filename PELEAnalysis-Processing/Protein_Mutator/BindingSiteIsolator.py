# -*- coding: utf-8 -*-

# Global imports
import sys,glob
import numpy as np
import argparse as ap

# Script information
__author__ = "Sergi Rodà"
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Sergi Rodà"
__email__ = "sergi.rodallordes@bsc.es"

class BindingSiteIsolator(object):


	def __init__(self):

		self.__filename, self.__resid, self.__radius, self.__NC = self.parseArgs() 
		self.__lines = self.__PDBParser()

	def __PDBParser(self):
		"""
		This method returns the lines of the PDB file
		"""

		PDB = open(self.__filename)
		lines = PDB.readlines(); PDB.close()

		return lines

	def parseArgs(self):
		"""
		Parse arguments from command-line
		"""

		parser = ap.ArgumentParser(description='Script used to output the PDB file around some specified reisdues')
		optional = parser._action_groups.pop()
		required = parser.add_argument_group('required arguments')
		parser.add_argument("input", metavar="FILE",type=str, help="path to PDB file")
		required.add_argument("-r","--residues",required=True,metavar="STRING",
								type=int,nargs='*',help="residue indices")
		optional.add_argument("-R","--radius",metavar="INTEGER",type=int,
								help="radius to the centroid",default=10)
		optional.add_argument("-NC","--not_centroid",help="Not centroid based",action = "store_true")
		parser._action_groups.append(optional)
		args = parser.parse_args()

		self.__filename =args.input
		self.__resid = args.residues
		self.__radius = args.radius
		self.__NC = args.not_centroid

		return self.__filename, self.__resid, self.__radius, self.__NC

	def GetBindingSiteCoordinates(self):
		"""
		This method takes the specified residues in the BindingSiteIsolator class
		instance and return the centroid of their coordinates

		RETURNS
		-------
		coords : 3D coordinates of the centroid (list of floats)
		"""

		coords = []

		for line in self.__lines:
			if (line[0:4] == "ATOM"  or line[0:6] == "HETATM"):
				if int(line[22:26].strip()) in self.__resid:
					x = float(line[30:38].strip())
					y = float(line[38:46].strip())
					z = float(line[46:54].strip())
					coords.append([x,y,z])

		# Find the centroid of all the coordinates of the selected residues
		coords = list(np.mean(coords,axis = 0))

		return coords

	def GetNeighboringAtoms(self,coords):
		"""
		This method finds the neighboring residues around the centroid of the coordinates
		found in the GetBindingSiteCoordinates method according to a specified radius in 
		the BindingSiteIsolator class

		PARAMETERS
		----------
		coords : 3D coordinates of the centroid (list of floats)

		RETURNS
		-------
		neighboring_residues : list of residues around the centroid of the binding site (list of string and int)
		"""

		def ComputeDistance(atom1, atom2):
			"""Computes the module or distance between the two points"""

			r = [abs(atom2[0] - atom1[0]), abs(atom2[1] - atom1[1]), abs(atom2[2] - atom1[2])]
			return np.sqrt(r[0] ** 2 + r[1] ** 2 + r[2]**2)

		neighboring_residues = []
		for line in self.__lines:
			if (line[0:4] == "ATOM"  or line[0:6] == "HETATM"):
				x = float(line[30:38].strip())
				y = float(line[38:46].strip())
				z = float(line[46:54].strip())
				resid_name = line[17:20].strip()
				resid_num = line[22:26].strip()
				atom_coords = [x,y,z]
				if ComputeDistance(atom_coords, [coords[0],coords[1],coords[2]]) <= self.__radius:
					residue = [str(resid_name), int(resid_num)]
					if residue not in neighboring_residues:
						neighboring_residues.append(residue)

		return neighboring_residues

	def CreateBindingPDBFile2(self):
		"""
		This method outputs a PDB file with just the residues specified
		in the BindingSiteIsolator class instance

		RETURNS
		-------
		PDB file of the specified residues
		"""

		output = open("{}_binding_site.pdb".format(self.__filename[:-4]), "w")

		for line in self.__lines:
			if (line[0:4] == "ATOM"  or line[0:6] == "HETATM"):
				if int(line[22:26].strip()) in self.__resid:
					output.write(line)

	def CreateBindingPDBFile(self,neighboring_residues):
		"""
		This method outputs a PDB file with the neighboring residues found in the 
		GetNeighboringAtoms method

		PARAMETERS
		----------
		neighboring_residues : list of residues around the centroid of the binding site (list of string and int)

		RETURNS
		-------
		PDB file of the neighboring residues
		"""

		output = open("{}_binding_site.pdb".format(self.__filename[:-4]), "w")

		for line in self.__lines:
			if (line[0:4] == "ATOM"  or line[0:6] == "HETATM"):
				resid = line[17:20].strip()
				number = line[22:26].strip()
				for elem in neighboring_residues:
					if (int(elem[1]) == int(number)) and (elem[0] == resid):
						output.write(line)

	def main(self):
		"""
		Main function

		It is called when this script is the main program called by the interpreter
		"""

		neighboring_residues = self.GetNeighboringAtoms(self.GetBindingSiteCoordinates())
		
		if self.__NC:
			self.CreateBindingPDBFile2()
		else:
			self.CreateBindingPDBFile(neighboring_residues)
			

if __name__=="__main__":
	"""Call the main function"""
	A = BindingSiteIsolator()
	A.main()
