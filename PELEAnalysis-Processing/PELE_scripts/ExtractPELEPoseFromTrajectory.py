# -*- coding: utf-8 -*-

# Global imports
from __future__ import unicode_literals
import os
import argparse as ap
import mdtraj as md

# Script information
__author__ = "Sergi Rodà"
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Sergi Rodà"
__email__ = "sergi.rodallordes@bsc.es"

# Functions
def parseArgs():
    """
    Parse arguments from command-line.

    RETURNS
    -------
    trajectory : string
              trajectory filename where to extract the pose
    metric : index (or list of indices)
                  column index (or list of column indices)
    output : string
                  output filename where the accepted PELE step will be saved
    """

    parser = ap.ArgumentParser(description='Script that extracts a/some accepted PELE step/s from the input trajectory file.')
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", required=True, metavar="FILE",
                          type=str, help="trajectory filename")
    required.add_argument("-N","--number", metavar="INTEGER",type=int,
                          help="number of the accepted PELE step that will be extracted from the trajectory file")
    optional.add_argument("-SN","--second_number", metavar="INTEGER",type=int,
                          help="number of the final accepted PELE step that will be extracted from the trajectory file (when a chunk of the trajectory wants to be extracted, SN > N)", default=0)    
    optional.add_argument("-o", "--output", metavar="FILE", type=str,
                          help="output path to save the pose", default="traj")
    optional.add_argument("-PN", "--particular_name",
                          help="specify if just the string given in the output flag will be used as the output filename", action="store_true")    
    optional.add_argument("-T", "--topology", metavar="FILE", type=str,
                          help="Topology file for the PELE trajectory in case it is in xtc format")
    parser._action_groups.append(optional)
    args = parser.parse_args()

    return args.input, args.number, args.second_number, args.output, args.particular_name, args.topology

def ExtractPoseFromTrajectory(trajectory, number, second_number, output, particular_name, topology):
    """
    This function extracts the lines of the trajectory file that
    refer to the desired accepted PELE step/s and saves it in an
    output file.

    OUTPUT
    ------
    The PDB file with the accepted PELE step/s
    """

    if ".xtc" in trajectory:
        traj = md.load_xtc(trajectory, topology)
        if second_number==0:
            if particular_name:
                traj[number-1].save_pdb("{}.pdb".format(output))
            else:
                traj[number-1].save_pdb("{}_{}_{}.pdb".format(output,trajectory.split("/")[-1].split("_")[-1].split(".")[0],number))
        else:
            if particular_name:
                traj[number-1:second_number].save_pdb("{}.pdb".format(output))
            else:
                traj[number-1:second_number].save_pdb("{}_{}_from_{}_to_{}.pdb".format(output,trajectory.split("/")[-1].split("_")[-1].split(".")[0],number,second_number))

    else:
        write_lines_down = False

        if particular_name:
            output_file = open("{}.pdb".format(output),"w")
        else:
            if second_number==0:
                output_file = open("{}_{}_{}.pdb".format(output,trajectory.split("/")[-1].split("_")[-1].split(".")[0],number),"w")
            else:
                output_file = open("{}_{}_from_{}_to_{}.pdb".format(output,trajectory.split("/")[-1].split("_")[-1].split(".")[0],number,second_number),"w")

        with open(trajectory, 'r') as traj_file:
            for line in traj_file:
                if line[0:5]=="MODEL":
                    if int(line.split()[1])==number:
                        write_lines_down = True
                    elif int(line.split()[1])==number+1 and second_number==0:
                        write_lines_down = False
                    elif int(line.split()[1])==second_number+1:
                        write_lines_down = False
                if write_lines_down == True:
                    output_file.write(line)

def main():
    """
    Main function

    It is called when this script is the main program called by the interpreter.
    """

    # Parse command-line arguments
    trajectory, number, second_number, output, particular_name, topology = parseArgs()

    if second_number!=0 and second_number < number:
    	print("Remember that to make a chunk of the trajectory, the second or last number must be logically higher than the first given number (refering to the first accepted PELE step)")
    	print("For more information on how to employ the code, use the -h or --help flag")
    	exit()

    # Extract the accepted PELE step/s from the trajectory file
    ExtractPoseFromTrajectory(trajectory, number, second_number, output, particular_name, topology)


if __name__ == "__main__":
    """Call the main function"""
    main()