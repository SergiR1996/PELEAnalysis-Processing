# -*- coding: utf-8 -*-

# Global imports
from __future__ import unicode_literals
import os,re
import glob
import math
import argparse as ap
import mdtraj as md

# Local imports
from PELEParseReports import *

# Script information
__author__ = "Sergi Rodà"
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Sergi Rodà"
__email__ = "sergi.rodallordes@bsc.es"

# Functions
def parseArgs():
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
                          type=str, help="path to report files")
    required.add_argument("-R", "--residue", metavar="INTEGER", type=str,nargs='*',
                          help="residue number of the residues involved in the metric")
    required.add_argument("-M", "--metric", metavar="STRING", type=str,nargs='*',
                          help="desired metric to calculate with the atom names (indicate with _)")
    required.add_argument("-RN", "--res_name", metavar="STRING", type=str,nargs='*',
                          help="name of the two residues involved in the distance metric to avoid confusions")
    required.add_argument("-CH", "--chain_name", metavar="STRING", type=str,nargs='*',
                          help="name of the chains where the two residues involved in the distance metric to avoid confusions")
    optional.add_argument("-CN","--column_name", metavar="STRING",type=str,
                          help="column name of the new metric", default="New_metric")
    optional.add_argument("-RF","--report_format", metavar="STRING",type=str,
                          help="Format of the report file", default="*out")
    optional.add_argument("-TF","--trajectory_format", metavar="STRING",type=str,
                          help="Format of the trajectory file", default="*pdb")
    optional.add_argument("-T", "--topology", metavar="FILE", type=str,
                          help="Topology file for the PELE trajectory in case it is in xtc format")
    parser._action_groups.append(optional)
    args = parser.parse_args()

    return args.input, args.residue, args.metric, args.res_name, args.chain_name, args.column_name, args.report_format, args.trajectory_format, args.topology

def AddMetric(out_directory, residue, metric, res_name, chain_name, column_name, report_format, trajectory_format, topology):
    """
    Take the PELE simulation trajectory files and returns the report files with the desired metric

    OUTPUT
    ------
    Report files with the desired metric added.
    """

    def ComputeDistance(atom1, atom2):

        r = [atom2[0] - atom1[0], atom2[1] - atom1[1], atom2[2] - atom1[2]]
        return math.sqrt(r[0] ** 2 + r[1] ** 2 + r[2]**2)
    
    def atoi(text):
        return int(text) if text.isdigit() else text
    
    def natural_keys(text):
        return [atoi(c) for c in re.split(r'(\d+)', text)]

    trajectory_list = glob.glob(os.path.join(out_directory,trajectory_format))
    report_list = glob.glob(os.path.join(out_directory,report_format))
    report_list.sort(key=natural_keys)
    trajectory_list.sort(key=natural_keys)
    metric_list = []

    for trajectory in trajectory_list:
        if ".xtc" in trajectory:
            traj = md.load_xtc(trajectory, topology)
            Atom_pair_1 = int(traj.topology.select("resSeq {} and name {} and resname {} and chainid {}".format(residue[0],metric[0].strip("_"),res_name[0],chain_name[0])))
            Atom_pair_2 = int(traj.topology.select("resSeq {} and name {} and resname {} and chainid {}".format(residue[1],metric[1].strip("_"),res_name[1],chain_name[1])))
            metric_list.append(10*md.compute_distances(traj,[[Atom_pair_1,Atom_pair_2]]))
        else:
            with open(trajectory, 'r') as traj_file:
                m_list,i = [],0
                for line in traj_file:

                    if (line[0:4] == "ATOM" or line[0:6] == "HETATM") and (line[22:26].strip() == residue[0] and line[12:16].replace(" ","_") == metric[0] and line[17:20] == res_name[0] and line[21:22] == chain_name[0]) or (line[22:26].strip() == residue[1] and line[12:16].replace(" ","_") == metric[1] and line[17:20] == res_name[1] and line[21:22] == chain_name[1]):

                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())

                        if i==0:
                            first_elem = [x,y,z]
                            i+=1
                        else:
                            second_elem = [x,y,z]
                            m_list.append(ComputeDistance(first_elem,second_elem))
                            i=0

                metric_list.append(m_list)

    for report in report_list:
        print(report)
        #print(report)

        if 'j' not in locals():
            j=0
        else:
            j+=1

        print(j)

        with open(report, 'r') as report_file:
            if report.endswith(".out"):
                out_report = open(report.split(".out")[0]+"_metric.out",'w')
            else:
                out_report = open(report+"_metric.out",'w')
            for i,line in enumerate(report_file):
                if i==0:
                    out_report.write(line.strip("\n")+'    '+column_name+"\n")
                else:
                    if ".xtc" in trajectory:
                        out_report.write(line.strip("\n")+'    '+str(metric_list[j][i-1][0])+"\n")
                    else:
                        out_report.write(line.strip("\n")+'    '+str(metric_list[j][i-1])+"\n")
    
    print("{} report files have the desired metric added".format(j+1))

def main():
    """
    Main function

    It is called when this script is the main program called by the interpreter
    """

    # Parse command-line arguments
    out_directory, residue, metric, res_name, chain_name, column_name, report_format, trajectory_format, topology  = parseArgs()

    # Add the desired metric to the report file
    AddMetric(out_directory, residue, metric, res_name, chain_name, column_name, report_format, trajectory_format, topology)


if __name__ == "__main__":
    """Call the main function"""
    main()


