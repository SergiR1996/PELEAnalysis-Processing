# -*- coding: utf-8 -*-

# Global imports
import argparse as ap
import numpy as np
from sklearn import cluster

# Local imports 
from MDAnalysisTools import *

def parseArgs():
    """
    Parse arguments from command-line

    RETURNS
    -------
    reports : string
              list of report files to look for data
    output_path : string
                  output directory where the csv file will be saved
    catalytic_event : list of indices
                  list of column indices
    to_drop : list of strings
                  list of column names that want to be dropped
    """

    parser = ap.ArgumentParser(description='Script that returns the density of the clusters of a selected atom in a trajectory')
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    parser.add_argument("traj", metavar="FILE",
                          type=str, help="path to trajectory file")
    parser.add_argument("top", metavar="FILE",
                          type=str, help="path to topology file")
    required.add_argument("-R", "--reference", required=True, metavar="FILE",
                          type=str, help="path to PDB reference file")
    optional.add_argument("-n", "--n_processors", metavar="INTEGER",
                          type=int, help="number of processors to make the analysis", default = None)
    optional.add_argument("-CW", "--cluster_width", metavar="FLOAT",
                          type=float, help="cluster width in $\AA$", default = 1.5)
    optional.add_argument("-RW", "--ref_atoms", metavar="LIST",
                          type=str,nargs='*', help="the reference atoms",default = 2)
    optional.add_argument("-SW", "--sim_atoms", metavar="LIST",
                          type=str,nargs='*', help="the simulation atoms", default = 2)
    optional.add_argument("-AN","--atom_name", metavar="STRING",
    					  type=str, help="the PDB atom name of the reference PDB file", default = "_CA_")
    optional.add_argument("-o","--output", metavar="PATH", type=str,help="filename of the output file", default="centroid")
    parser._action_groups.append(optional)
    args = parser.parse_args()

    traj, top, reference, n_processors, cluster_width, atom_ref_ids, atom_sim_ids, atom_name, output = args.traj, args.top, args.reference, args.n_processors, args.cluster_width, args.ref_atoms, args.sim_atoms, args.atom_name, args.output

    return traj, top, reference, n_processors, cluster_width, atom_ref_ids, atom_sim_ids, atom_name, output

def GetCoordinates(file,atom_ref_ids, atom_name):
	"""
	This function returns the coordinates of a PDB file.

	PARAMETERS
    ----------
    file : string
            PDB filename

    RETURNS
    ------
    res_coord: list of floats
	"""

	PDB = open(file)
	res_coord = []
	for line in PDB:
		if (line.strip()[0:4] == "ATOM") and (line.strip()[12:16].replace(" ","_") == atom_name) and (line.strip()[22:26].replace(" ","") in atom_ref_ids):
			x = float(line.strip()[30:38])
			y = float(line.strip()[38:46])
			z = float(line.strip()[46:54])
			res_coord.append([x,y,z])

	return res_coord

def ClusterizeAtoms(traj,top,reference,atom_sim_ids,atom_ref_ids, cluster_width, n_processors, atom_name, output):
	"""
	This function takes the arguments from the command line and clusterizes the position of the desired atoms 
	during the MD simulation according to the reference structure.

	PARAMETERS
	----------
	All command line arguments

	OUTPUT
	------
	A PDB file with the centroids of the average positions of the atom during the trajectory
	and the density is stored in the B-factor.
	"""


	# Open trajectory file with topology and extract interesting atom coordinates
	traj_aux = OpenFiles(traj, top)
	trajectory = traj_aux.load_trajectory()

	Atom_indices = ""
	for elem in atom_sim_ids:
		if elem==atom_sim_ids[0] and len(atom_sim_ids)==1:
			Atom_indices+="(resSeq {})".format(elem)
		elif elem==atom_sim_ids[0] and len(atom_sim_ids)==2:
			Atom_indices+="(resSeq {} or ".format(elem)
		elif elem==atom_sim_ids[0]:
			Atom_indices+="(resSeq {}".format(elem)
		elif elem==atom_sim_ids[len(atom_sim_ids)-1]:
			Atom_indices+="resSeq {})".format(elem)
		else:
			Atom_indices+=" or resSeq {} or ".format(elem)

	Atoms = trajectory.topology.select("name %s and %s" %(atom_name.strip("_"),Atom_indices))
	Atom_coordinates = []
	for elem in Atoms:
		for model in trajectory.xyz[:,elem,:]:
			Atom_coordinates.append(model*10) # The multiplier by 10 is to convert to Angstrom units.

	# Retrieve reference data for cluster analysis
	Ref_coords = GetCoordinates(reference, atom_ref_ids, atom_name)

	# Clustering
	Estimator = cluster.MeanShift(bandwidth = cluster_width, n_jobs = n_processors, cluster_all = True)
	Results = Estimator.fit_predict(Atom_coordinates)

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

    traj, top, reference, n_processors, cluster_width, atom_ref_ids, atom_sim_ids, atom_name, output = parseArgs()

    ClusterizeAtoms(traj, top, reference, atom_sim_ids, atom_ref_ids, cluster_width, n_processors, atom_name, output)

if __name__=="__main__":
	main()
