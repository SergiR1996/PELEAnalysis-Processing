# -*- coding: utf-8 -*-

# Global Imports 
import os,sys
import glob
import argparse as ap
import multiprocessing as mp # Module to parallelize job.
import argparse

# Schrödinger software import
from schrodinger import structure
from schrodinger.application.bioluminate import protein
from schrodinger.structutils import minimize

# Script information
__author__ = "Sergi Rodà"
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Sergi Rodà"
__email__ = "sergi.rodallordes@bsc.es"

class ProteinMutator():

    def __init__(self):

        self.__filename, self.__mutated_residues, self.__simultaenous_mutations, self.__refinement_range = self.parseArgs()
        self.__mutated_residues = [i.split("_") for i in self.__mutated_residues]
        self.__original_filename = self.__filename

    def parseArgs(self):
        """
        Parse arguments from command-line

        RETURNS
        -------
        input : string
                  list of report files to look for data
        mutated_residues : list of tuples
                      List of specified residues to be mutated and to which residue
        simultaenous_mutations: integer
                      Apply all mutations simultaneously
        refinement_range: float
                      Radius around the mutated residue to be refined with Prime
        """
        parser = ap.ArgumentParser(description='Script that perform the specified mutations \
            to a PDB file using Schrödinger suite.')
        optional = parser._action_groups.pop()
        required = parser.add_argument_group('required arguments')
        required.add_argument("-i", "--input", required=True, metavar="FILE",
                              type=str, help="path of the input PDB file")
        optional.add_argument("-M","--mutated_residues", metavar="LIST",type=str,
                              nargs='*',help="Residues that will be mutated (format: RES_NUM; GLU_123)")
        optional.add_argument("-S","--simultaenous_mutations", metavar="INTEGER",type=int,
                              help="Apply all the changes to the WT structure simultaneously", default=0)
        optional.add_argument("-R","--refinement_range", metavar="FLOAT",type=float,
                              help="Radius around the mutated residue to be refined with Prime", default=4.0)
        parser._action_groups.append(optional)
        args = parser.parse_args()

        return args.input, args.mutated_residues, args.simultaenous_mutations, args.refinement_range

    @property
    def filename(self):
        return self.__filename

    def get_residue_atoms(self,index):
        
        PDB = open(self.__filename)
        for line in PDB:
            if (line[0:4] == "ATOM"):
                if index == int(line[22:26].strip()):
                    return int(line[6:11].strip())

    def apply_mutation(self,index,final_residue):
        """
        Take a PDB file and mutate a particular residue of the system.

        PARAMETERS
        ----------
        index : integer
                  Index of the list of mutated residues
        final_residue : string
                  Residue name for the residue that wants to be mutated

        OUTPUT
        ------
        The PDB file with the mutated residues
        """

        if self.__simultaenous_mutations != 1:
            self.__filename = self.__original_filename

        # Convert PDB input file into mae file for proper use
        os.system("$SCHRODINGER/utilities/structconvert -ipdb %s -omae %s.mae"%(self.__filename,self.__filename[:-4]))
        Input_filename = "%s.mae"%(self.__filename[:-4])


        # Open the input mae file and create the output mae file where mutations will be added
        input_structure = structure.Structure.read(Input_filename)
        output_structure = structure.StructureWriter("mutated_%s"%(Input_filename))
        # output_structure.append(input_structure) # To append the input structure on the output file

        # Get the indices of the residues that want to be mutated and generate the list of the residues with their mutations
        Mutations = [('A',index,' ',final_residue)]

        # Generate the tool to create the mutated structures
        mutator = protein.Mutator(input_structure, Mutations)

        # Iterate over each mutation for each residue
        for mutation in mutator.generate():

            mutated_structure = mutation.struct
            residue_map       = mutation.residue_map

            res_str = ", ".join(str(res) for res in residue_map.values())
            print ('Residue affected by this mutation: %s' %(res_str))

            # Overwrite an output file of the mutated structure
            output_structure.append(mutated_structure)

        # os.system("$SCHRODINGER/utilities/prepwizard -epik_pH {} -propka_pH {} {} mutated_and_minimized_{}.pdb".format(self.__pH,self.__pH,"mutated_{}".format(Input_filename),self.filename[:-4]))

        os.system("$SCHRODINGER/utilities/structconvert -imae %s -opdb %s.pdb"%("mutated_%s"%(Input_filename),self.__filename[:-4]+"_"+final_residue))

        if self.__simultaenous_mutations != 1:
            self.__filename = self.__original_filename[:-4]+"_"+final_residue+".pdb"
        else:
            self.__filename = self.__filename[:-4]+"_"+final_residue+".pdb"

    def refine_mutation(self,index):
        """
        Take the mutated PDB file and refines the structure
        using Prime around the residue specified by the 
        index.

        PARAMETERS
        ----------
        index : integer
                  Index of the the residue around the one the refinement will be performed

        OUTPUT
        ------
        The refined mutated pdb file
        """

        os.system("$SCHRODINGER/utilities/structconvert -ipdb %s -omae %s.mae"%(self.__filename,self.__filename[:-4]))
        Input_filename = "%s.mae"%(self.__filename[:-4])

        st = structure.Structure.read(Input_filename)

        # Get the residue that was mutated
        ca = st.atom[self.get_residue_atoms(index)]
        mutated_residue = None
        for res in st.residue:
            ca_keys  = (ca.chain,  ca.resnum,  ca.inscode)
            res_keys = (res.chain, res.resnum, res.inscode)
            if ca_keys == res_keys:
                mutated_residue = res
                break

        # We want to use the reference to gather the residues to refine
        refine_residues = protein.get_residues_within(st,[mutated_residue],within = self.__refinement_range)

        # Create the refiner
        refiner = protein.Refiner(st, residues=refine_residues)

        # Run Prime minimization which returns the refined structure
        refined_struct = refiner.runPrimeMinimization('Refine')

    def minimization(self):
        """
        The method applies a minimization to a mutated pdb file
        """

        os.system("$SCHRODINGER/utilities/structconvert -ipdb mutated_{} -omae aux_{}.mae".format(self.__filename,self.__filename[:-4]))
        Input_filename = "aux_{}.mae".format(self.__filename[:-4])

        input_structure = structure.Structure.read(Input_filename)
        output_structure = structure.StructureWriter("minimized_{}".format(Input_filename))

        minimize.minimize_structure(input_structure)
        output_structure.append("min"+Input_filename)

        os.system("$SCHRODINGER/utilities/structconvert -imae {} -opdb minimized_{}.pdb".format("minimized_{}".format(Input_filename),self.__filename[:-4]))

    def main(self):
        """
        Main function
    
        It is called when this script is the main program called by the interpreter
        """

        for i, mutatable_residue in enumerate(self.__mutated_residues):
            self.apply_mutation(int(mutatable_residue[1]),mutatable_residue[0])
            if self.__simultaenous_mutations != 1:
                self.refine_mutation(int(mutatable_residue[1]))
                os.system("$SCHRODINGER/utilities/structconvert -imae refinement_job_%s/Refine-out.maegz -opdb %s%s%s.pdb"%(i+1,self.__original_filename[:-4],mutatable_residue[0],mutatable_residue[1]))
        if self.__simultaenous_mutations == 1:
            self.refine_mutation(int(self.__mutated_residues[-1][1]))
            os.system("$SCHRODINGER/utilities/structconvert -imae refinement_job_1/Refine-out.maegz -opdb %s.pdb"%(self.__filename[:-4]))


if __name__ == "__main__":
    """Call the main function"""
    PM = ProteinMutator()
    PM.main()

