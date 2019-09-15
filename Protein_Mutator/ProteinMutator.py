# -*- coding: utf-8 -*-

# Imports 
import os,sys
import glob
import argparse as ap
import multiprocessing as mp # Module to parallelize job.

# Schrödinger software import
from schrodinger import structure
from schrodinger.application.bioluminate import protein
from schrodinger.structutils import minimize

# Script information
__author__ = "Sergi Rodà Llordés"
__version__ ="1.0"
__maintainer__="Sergi Rodà Llordés"
__email__="sergi.rodallordes@bsc.es"

class ProteinMutator():

    def __init__(self,filename = "", mutated_residues = [0], pH = 7, simultaenous_mutations = 1):

        self.__filename = filename
        self.__mutated_residues = mutated_residues
        self.__pH = pH
        self.__simultaenous_mutations = simultaenous_mutations

    @property
    def filename(self):
        return self.__filename

    def get_residue_atoms(self,index):
        PDB = open(self.__filename)
        for line in PDB:
	    if "ATOM" in line:
            	if self.__mutated_residues[index] == int(line.split()[5]):
                	return int(line.split()[1])

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

        # Convert PDB input file into mae file for proper use
        os.system("$SCHRODINGER/utilities/structconvert -ipdb {} -omae {}.mae".format(self.__filename,self.filename[:-4]))
        Input_filename = "{}.mae".format(self.__filename[:-4])


        # Open the input mae file and create the output mae file where mutations will be added
        input_structure = structure.Structure.read(Input_filename)
        output_structure = structure.StructureWriter("mutated_{}".format(Input_filename))
        # output_structure.append(input_structure) # To append the input structure on the output file

        # Get the indices of the residues that want to be mutated and generate the list of the residues with their mutations
        Mutations = [('A',self.__mutated_residues[index],' ',final_residue)]

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

        os.system("$SCHRODINGER/utilities/structconvert -imae {} -opdb {}.pdb".format("mutated_{}".format(Input_filename),self.filename[:-4]+"_"+final_residue))

        if self.__simultaenous_mutations == 1:
            self.__filename = self.filename[:-4]+"_"+final_residue+".pdb"

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

        os.system("$SCHRODINGER/utilities/structconvert -ipdb {} -omae {}.mae".format(self.__filename,self.filename[:-4]))
        Input_filename = "{}.mae".format(self.__filename[:-4])

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
        refine_residues = protein.get_residues_within(st,[mutated_residue],within = 10.0)

        # Create the refiner
        refiner = protein.Refiner(st, residues=refine_residues)

        # Run Prime minimization which returns the refined structure
        refined_struct = refiner.runPrimeMinimization('Refine')

    def minimization(self):
        """
        The method applies a minimization to a mutated pdb file
        """

        os.system("$SCHRODINGER/utilities/structconvert -ipdb mutated_{} -omae aux_{}.mae".format(self.__filename,self.filename[:-4]))
        Input_filename = "aux_{}.mae".format(self.__filename[:-4])

        input_structure = structure.Structure.read(Input_filename)
        output_structure = structure.StructureWriter("minimized_{}".format(Input_filename))

        minimize.minimize_structure(input_structure)
        output_structure.append("min"+Input_filename)

        os.system("$SCHRODINGER/utilities/structconvert -imae {} -opdb minimized_{}.pdb".format("minimized_{}".format(Input_filename),self.filename[:-4]))

    def main(self):
        """
        Main function
    
        It is called when this script is the main program called by the interpreter
        """

        self.apply_mutation(0,'SER')
        self.apply_mutation(1,'HID')
        self.apply_mutation(2,'ASP')
        self.apply_mutation(3,'GLY')
        self.apply_mutation(4,'GLY')
        self.refine_mutation(1)
        # self.apply_mutation(2,'ASP')


if __name__ == "__main__":
    """Call the main function"""
    Trial = ProteinMutator("A.pdb",[153,334,345,336,355],1)
    Trial.main()

