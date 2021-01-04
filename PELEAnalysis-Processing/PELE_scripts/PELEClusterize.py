# -*- coding: utf-8 -*-

# Global imports
import argparse as ap
import numpy as np
from sklearn import cluster

# Local variable
Backbone_atoms = ["_CA_","_C__","_N__","_O__","_OXT","_HA_","_HA2","_HA3","_H__","_H1_","_H2_","_H3_"]

# Script information
__author__ = "Sergi Rodà"
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Sergi Rodà"
__email__ = "sergi.rodallordes@bsc.es"

def parseArgs():
	"""
	Parse arguments from command-line

	RETURNS
	-------
	traj : string
					list of trajectory files to analyze
	reference : string
					Path to the PDB reference file
	n_processors : integer
					number of processors to make the analysis
	cluster_width : float
					cluster width used in MeanShift algorithm
	non_ligand : boolean
					the boolean to decide where the ligand is clusterized or not
	residue_number : list of strings
					the residue number of the clusterized protein residues
	only_sidechain : boolean
					clusterize only the side chain atoms of the residues if true
	output : string
					filename of the output file					
	"""

	parser = ap.ArgumentParser(description='Script that returns the density of the clusters of the ligand in a PELE trajectory')
	optional = parser._action_groups.pop()
	required = parser.add_argument_group('required arguments')
	parser.add_argument("traj", metavar="FILE",
							type=str, help="path to trajectory file")
	required.add_argument("-R", "--reference", required=True, metavar="FILE",
							type=str, help="path to PDB reference file")
	optional.add_argument("-n", "--n_processors", metavar="INTEGER",
							type=int, help="number of processors to make the analysis", default = None)
	optional.add_argument("-CW", "--cluster_width", metavar="FLOAT",
							type=float, help="cluster width in $\AA$", default = 1.5)
	optional.add_argument("-NL","--non_ligand",
							help = "Clusterize against a residue of the protein instead of the ligand", action = "store_true")
	optional.add_argument("-RN", "--residue_number", metavar="INTEGER",
							type=str, nargs='*', help="the residue numbers of the clusterized protein residues",default = [])
	optional.add_argument("-OS","--only_sidechain",
							help = "Clusterize only the atoms of the side chain of the residue/s", action = "store_true")	
	optional.add_argument("-o","--output", metavar="PATH", type=str,help="filename of the output file", default="centroid")
	parser._action_groups.append(optional)
	args = parser.parse_args()

	return args.traj, args.reference, args.n_processors, args.cluster_width, args.non_ligand, args.residue_number, args.only_sidechain, args.output

def GetCoordinates(file, non_ligand, residue_number, only_sidechain, reference = False):
	"""
	This function returns the coordinates of a PDB file.

	PARAMETERS
	----------
	file : string
			PDB filename
	non_ligand : boolean
					the boolean to decide where the ligand is clusterized or not
	residue_number : list of strings
					the residue number of the clusterized protein residues
	only_sidechain : boolean

	RETURNS
	------
	res_coord: list of floats
	"""

	PDB = open(file)
	res_coord, res_coord_aux = [], []
	for line in PDB:
		if line[0:5]=="MODEL":
			if len(res_coord_aux) == 0:
				pass
			else:
				res_coord.append(list(np.mean(res_coord_aux,axis=0)))
				res_coord_aux = []
		if (line[21:22] == "L") and (non_ligand == False):
			x = float(line[30:38].strip())
			y = float(line[38:46].strip())
			z = float(line[46:54].strip())
			res_coord_aux.append([x,y,z])

		elif (non_ligand == True) and (line[22:26].strip() in residue_number) and (only_sidechain == False):
			x = float(line[30:38].strip())
			y = float(line[38:46].strip())
			z = float(line[46:54].strip())
			res_coord_aux.append([x,y,z])

		elif (non_ligand == True) and (line[22:26].strip() in residue_number) and (only_sidechain == True) and not (line[12:16].replace(" ","_") in Backbone_atoms):
			x = float(line[30:38].strip())
			y = float(line[38:46].strip())
			z = float(line[46:54].strip())
			res_coord_aux.append([x,y,z])

	if reference:
		res_coord = res_coord_aux

	return res_coord

def ClusterizeAtoms(traj, reference, non_ligand, residue_number, cluster_width, n_processors, only_sidechain, output):
	"""
	This function takes the arguments from the command line and clusterizes the position of the ligand 
	(or the selected protein residues) during the PELE simulation according to the reference structure.

	PARAMETERS
	----------
	All command line arguments

	OUTPUT
	------
	A PDB file with the centroids of the average positions of the atom during the trajectory
	and the density is stored in the B-factor.
	"""

	# Retrieve trajectory data for cluster analysis
	Coords = GetCoordinates(traj, non_ligand, residue_number, only_sidechain)

	# Retrieve reference data for cluster analysis
	Ref_coords = GetCoordinates(reference, non_ligand, residue_number, only_sidechain, True)

	# Clustering
	Estimator = cluster.MeanShift(bandwidth = cluster_width, n_jobs = n_processors, cluster_all = True)
	Results = Estimator.fit_predict(Coords)

	Ref_clusters = []
	for ref_coord in Ref_coords:
				Ref_clusters += Estimator.predict([ref_coord]).tolist()
	        
	# Clustering analysis
	aux = np.bincount(Results)
	aux_count = np.nonzero(aux)[0]
	density = {}
	for i in aux_count:
		density[i]=aux[i]/np.sum(aux)
	
	# Write centroids to PDB
	centroids = Estimator.cluster_centers_
	n_clusters = len(centroids)
	with open(output+".pdb", 'w') as f:
		for i, centroid in enumerate(centroids):
			f.write("ATOM    {:3d}  CEN BOX A {:3d} {:>11.3f}{:>8.3f}{:>8.3f}  1.00{:>5.2f}\n".format(i+1, i+1, *centroid, density[i]))

def main():
	"""
	Main function

	It is called when this script is the main program called by the interpreter
	"""

	traj, reference, n_processors, cluster_width, non_ligand, residue_number, only_sidechain, output = parseArgs()

	ClusterizeAtoms(traj, reference, non_ligand, residue_number, cluster_width, n_processors, only_sidechain, output)

if __name__=="__main__":
	main()