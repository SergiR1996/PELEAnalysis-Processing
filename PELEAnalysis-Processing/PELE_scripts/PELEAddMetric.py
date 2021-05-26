# -*- coding: utf-8 -*-

# Global imports
from __future__ import unicode_literals
import os,re
import glob
import math
import argparse as ap
import mdtraj as md
import multiprocessing as mp

# Script information
__author__ = "Sergi Rodà"
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Sergi Rodà"
__email__ = "sergi.rodallordes@bsc.es"

class MetricCalculator(): 
    
    def __init__(self):

        self.input, self.atom_index, self.residue, self.metric, self.res_name, self.column_name, \
        self.report_format, self.trajectory_format, self.topology, \
        self.output_name, self.verbose, self.n_processors = self.parseArgs()

    def parseArgs(self):
        """
        Parse arguments from command-line

        RETURNS
        -------
        out_directory : string
                  path of the PELE simulation output directory
        """

        parser = ap.ArgumentParser(description='Script that returns a report file with the desired metric from a PELE simulation')
        required = parser.add_argument_group('required arguments')
        optional = parser._action_groups.pop()
        required.add_argument("-i", "--input", required=True, metavar="DIRECTORY",
                              type=str, nargs='*', help="path to report files")
        optional.add_argument("-A", "--atom_index", metavar="INTEGER", type=str,nargs='*',
                              help="atom index of the atoms involved in the metric", default=[])
        optional.add_argument("-R", "--residue", metavar="INTEGER", type=str,nargs='*',
                              help="residue number of the residues involved in the metric", default=[])
        optional.add_argument("-M", "--metric", metavar="STRING", type=str,nargs='*',
                              help="desired metric to calculate with the atom names (indicate with _)")
        optional.add_argument("-RN", "--res_name", metavar="STRING", type=str,nargs='*',
                              help="name of the n residues involved in the metric to avoid confusions")
#        required.add_argument("-CH", "--chain_id", metavar="STRING", type=str,nargs='*',
#                              help="number of the chains where the n residues involved in the metric to avoid confusions")
        optional.add_argument("-CN","--column_name", metavar="STRING",type=str,
                              help="column name of the new metric", default="New_metric")
        optional.add_argument("-RF","--report_format", metavar="STRING",type=str,
                              help="Format of the report file", default="*out")
        optional.add_argument("-TF","--trajectory_format", metavar="STRING",type=str,
                              help="Format of the trajectory file", default="*pdb")
        optional.add_argument("-T", "--topology", metavar="FILE", type=str,
                              help="Topology file for the PELE trajectory in case it is in xtc format")
#        optional.add_argument("-TM","--type_of_metric", metavar="STRING",type=str,
#                              help="Type of analysis to perform", default="distance")
        optional.add_argument("-ON", "--output_name", metavar="STRING",type=str,
                              help="Extension in the name of the new report files", default="_metric.out")
        optional.add_argument("-V","--verbose",
                              help="Print some messages while executing the code to help the user", action="store_true")
        optional.add_argument("-NP","--n_processors", metavar="INTEGER",type=int,
                              help="number of processors to execute the code", default = 4)
        parser._action_groups.append(optional)
        args = parser.parse_args()

        return args.input, args.atom_index, args.residue, args.metric, args.res_name, args.column_name, args.report_format, \
        args.trajectory_format, args.topology, args.output_name, args.verbose, args.n_processors

    def DecompressList(self,l_of_lists):
        """
        This method decompress a 
        list of lists into a list

        PARAMETERS
        ----------
        l_of_lists: list of lists
                List of lists

        RETURNS
        -------
        new_list: list
                New list    
        """

        new_list = []
        for sublist in l_of_lists:
            for item in sublist:
                new_list.append(item)

        return new_list

    # def ComputeDistance(self, atom1, atom2):
    #     """
    #     This method computes the distance between two atoms

    #     PARAMETERS
    #     ----------
    #     atom1: XYZ coordinates of first atom
    #             List of floats
    #     atom2: XYZ coordinates of second atom
    #             List of floats

    #     RETURNS
    #     -------
    #     distance between the 2 atoms (float)  
    #     """

    #     r = [atom2[0] - atom1[0], atom2[1] - atom1[1], atom2[2] - atom1[2]]
    #     return math.sqrt(r[0] ** 2 + r[1] ** 2 + r[2]**2)

    # def ComputeAngle(self, atom1, atom2, atom3):
    #     """
    #     This method computes the angle between three atoms

    #     PARAMETERS
    #     ----------
    #     atom1: XYZ coordinates of first atom
    #             List of floats
    #     atom2: XYZ coordinates of second atom
    #             List of floats
    #     atom3: XYZ coordinates of third atom
    #             List of floats

    #     RETURNS
    #     -------
    #     angle between the 3 atoms (float)  
    #     """

    #     at12 = [atom2[0] - atom1[0], atom2[1] - atom1[1], atom2[2] - atom1[2]]
    #     at23 = [atom3[0] - atom2[0], atom3[1] - atom2[1], atom3[2] - atom2[2]]

    #     at_dot_product = (at12[0]*at23[0])+(at12[1]*at23[1])+(at12[2]*at23[2])
    #     at_12_length = math.sqrt(at12[0] ** 2 + at12[1] ** 2 + at12[2]**2)
    #     at_23_length = math.sqrt(at23[0] ** 2 + at23[1] ** 2 + at23[2]**2)

    #     return math.degrees(math.acos((at_dot_product)/(at_12_length*at_23_length)))

    def AddDistance(self, trajectory):
        """
        Take the PELE simulation trajectory files and returns the list of values of the desired distance metric

        RETURNS
        -------
        metric_list: list
                List of values of the desired distance metric
        """

        if ".xtc" in trajectory:
            traj = md.load_xtc(trajectory, self.topology)
        else:
            traj = md.load_pdb(trajectory)

        if len(self.atom_index) != 0:
            metric_list = 10*md.compute_distances(traj,[[int(self.atom_index[0]),int(self.atom_index[1])]])
        else:
            Atom_pair_1 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(self.residue[0],self.metric[0].strip("_"),self.res_name[0])))
            Atom_pair_2 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(self.residue[1],self.metric[1].strip("_"),self.res_name[1])))
            metric_list = 10*md.compute_distances(traj,[[Atom_pair_1,Atom_pair_2]])

            # with open(trajectory, 'r') as traj_file:
            #     m_list,i = [],0
            #     for line in traj_file:

            #         if (line[0:4] == "ATOM" or line[0:6] == "HETATM") and (line[22:26].strip() == self.residue[0] and line[12:16].replace(" ","_") == self.metric[0] \
            #         and line[17:20] == self.res_name[0] and line[21:22] == self.chain_name[0]) or (line[22:26].strip() == self.residue[1] \
            #         and line[12:16].replace(" ","_") == self.metric[1] and line[17:20] == self.res_name[1] and line[21:22] == self.chain_name[1]):

            #             x = float(line[30:38].strip())
            #             y = float(line[38:46].strip())
            #             z = float(line[46:54].strip())

            #             if i==0:
            #                 first_elem = [x,y,z]
            #                 i+=1
            #             else:
            #                 second_elem = [x,y,z]
            #                 m_list.append(self.ComputeDistance(first_elem,second_elem))
            #                 i=0

            #         metric_list.append(m_list)

        return metric_list

    def AddAngle(self, trajectory):
        """
        Take the PELE simulation trajectory files and returns the list of values of the desired angle metric

        RETURNS
        -------
        metric_list: list
                List of values of the desired angle metric
        """

        if ".xtc" in trajectory:
            traj = md.load_xtc(trajectory, self.topology)
        else:
            traj = md.load_pdb(trajectory)

        if len(self.atom_index) != 0:
            metric_list = md.compute_angles(traj,[[int(self.atom_index[0]),int(self.atom_index[1]), int(self.atom_index[2])]])
        else:
            Atom_pair_1 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(self.residue[0],self.metric[0].strip("_"),self.res_name[0])))
            Atom_pair_2 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(self.residue[1],self.metric[1].strip("_"),self.res_name[1])))
            Atom_pair_3 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(self.residue[2],self.metric[2].strip("_"),self.res_name[2])))
            metric_list = md.compute_angles(traj,[[Atom_pair_1,Atom_pair_2, Atom_pair_3]])

        return metric_list

    def AddDihedral(self, trajectory):
        """
        Take the PELE simulation trajectory files and returns the list of values of the desired dihedral metric

        RETURNS
        -------
        metric_list: list
                List of values of the desired dihedral metric
        """

        if ".xtc" in trajectory:
            traj = md.load_xtc(trajectory, self.topology)
        else:
            traj = md.load_pdb(trajectory)

        if len(self.atom_index) != 0:
            metric_list = md.compute_dihedrals(traj,[[int(self.atom_index[0]),int(self.atom_index[1]), int(self.atom_index[2]), int(self.atom_index[3])]])
        else:
            Atom_pair_1 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(self.residue[0],self.metric[0].strip("_"),self.res_name[0])))
            Atom_pair_2 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(self.residue[1],self.metric[1].strip("_"),self.res_name[1])))
            Atom_pair_3 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(self.residue[2],self.metric[2].strip("_"),self.res_name[2])))
            Atom_pair_4 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(self.residue[3],self.metric[3].strip("_"),self.res_name[3])))
            metric_list = md.compute_dihedrals(traj,[[Atom_pair_1,Atom_pair_2, Atom_pair_3, Atom_pair_4]])

        return metric_list

    def AddMetric(self):
        """
        Take the PELE simulation trajectory files and returns the report files with the desired metric

        OUTPUT
        ------
        Report files with the desired metric added.
        """
     
        def atoi(text):
            return int(text) if text.isdigit() else text
        
        def natural_keys(text):
            return [atoi(c) for c in re.split(r'(\d+)', text)]

        trajectory_list, report_list = [],[]

        for path in self.input: 
            trajectory_list += glob.glob(os.path.join(path,self.trajectory_format))
            report_list += glob.glob(os.path.join(path,self.report_format))
        report_list.sort(key=natural_keys)
        trajectory_list.sort(key=natural_keys)
        Metric_results = []

        pool = mp.Pool(self.n_processors)
        if len(self.atom_index) == 2 or len(self.residue) == 2:
            Metric_results.append(pool.map(self.AddDistance,trajectory_list))
        elif len(self.atom_index) == 3 or len(self.residue) == 3:
            Metric_results.append(pool.map(self.AddAngle,trajectory_list))
        elif len(self.atom_index) == 4 or len(self.residue) == 4:
            Metric_results.append(pool.map(self.AddDihedral,trajectory_list))
        pool.terminate()

        Metric_results = self.DecompressList(Metric_results)

        for report in report_list:
            if self.verbose:
                print(report)

            if 'j' not in locals():
                j=0
            else:
                j+=1

            if self.verbose:
                print(j)

            with open(report, 'r') as report_file:
                if report.endswith(".out"):
                    out_report = open(report.split(".out")[0]+self.output_name,'w')
                else:
                    out_report = open(report+self.output_name,'w')
                for i,line in enumerate(report_file):
                    if i==0:
                        out_report.write(line.strip("\n")+'    '+self.column_name+"\n")
                    else:
                        if len(self.atom_index) > 2 or len(self.residue) > 2:
                            out_report.write(line.strip("\n")+'    '+str(math.degrees(Metric_results[j][i-1][0]))+"\n")
                        else:
                            out_report.write(line.strip("\n")+'    '+str(Metric_results[j][i-1][0])+"\n")

        if self.verbose:
            print("{} report files have the desired metric added".format(j+1))

if __name__ == "__main__":
    """Call the main function"""
    Metriccalculator = MetricCalculator()
    Metriccalculator.AddMetric()