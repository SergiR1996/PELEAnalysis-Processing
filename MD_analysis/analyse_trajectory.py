# -*- coding: utf-8 -*-

"""

                            MD ANALYSER (for molecular dynamics trajectories)

This script is for analysing the trajectories obtained by any software. It will use the MDtraj module for such
analysis. There will be some flag options which the user can pick on and the files and plots will be saved in a new
created directory

                                    Done by: Ruben Canadas Rodriguez & Sergi Rodà Llordés

"""

# Let's import packages
import matplotlib
matplotlib.use("tkagg")
import mdtraj as md
import glob, os, sys, argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


#plt.style.use(['seaborn-darkgrid'])


class OpenFiles:

    def __init__(self, xtc_file, top_file, trr_file=None, file=None):

        self.xtc_file = xtc_file
        self.top_file = top_file
        self.trr_file = trr_file
        self.file = file       #Any file of format: pdb,xtc,trr,dcd,binpos,Netcdf,hdf5..

    def load_xtc(self):

        """

        :param xtc_file: Trajectory in xtc file
        :param top_file: Topology file of such trajectory
        :return: trajectory loaded
        """

        trajectories = []

        for traj in self.xtc_file:
            trajectories.append(md.load_xtc(traj, self.top_file))

        return trajectories

    def load_trr(self):

        return md.load_trr(self.trr_file, self.top_file)



    def load_trajectory(self):

        return md.load(self.xtc_file, self.top_file)



    def information(self, file):


        with open("trajectory_info.txt", "w") as outfile:
            outfile.write("Number of frames: {}\nnumber of atoms: {}\nnumber of residues: {}\number of chains: {}\n\n".format(file.n_frames,
            file.n_atoms, file.n_residues, file.n_chains))


    def number_frames(self, file):



        return file.n_frames



class directory_manager:


    def __init__(self, directory):

        self.directory = directory

    def create_directory(self):

        if not os.path.exists(self.directory):

            os.mkdir(self.directory)

        else:

            pass




class TrajectoryProperties:

    """
    This class will have methods for computing properties of a trajectory: rmsd, distances ..
    """

    def __init__(self, traj, metric="angstrom"):

        self.traj = traj
        self.metric = metric

    def nanometer_to_angstrom(func):

        """
        This will be a decoration method
        :param func: Function to decorate
        :return:

        """

        def new_metric(self, pairs):
            return func(self,  pairs) * 10

        return new_metric


    # @nanometer_to_angstrom
    def traj_rmsd(self, reference, atom_indices):

        return md.rmsd(self.traj, reference, atom_indices = atom_indices)

    #@nanometer_to_angstrom
    def traj_rmsf(self, reference, atom_indices):

        self.traj.superpose(reference, atom_indices = atom_indices)
        return md.rmsf(self.traj, None, atom_indices = atom_indices)



    @nanometer_to_angstrom
    def compute_distance(self, atom_pairs):

        """


        :param atom_pairs: Each row gives the indices of two atoms involved in the interaction: np.ndarray, shape=(num_pairs, 2), dtype=int
        :return: distances : np.ndarray, shape=(n_frames, num_pairs), dtype=float (The distance, in each frame, between each pair of atoms)

        """

        return md.compute_distances(self.traj, atom_pairs)


    @nanometer_to_angstrom  # Decorator for transforming array in nanometer to array in angstroms
    def compute_displacements(self, atom_pairs):

        """

        :param atom_pairs: Each row gives the indices of two atoms: np.ndarray, shape[num_pairs, 2], dtype=int
        :return:  displacements : np.ndarray, shape=[n_frames, n_pairs, 3], dtype=float32

        """

        return md.displacements(self.traj, atom_pairs)


    @nanometer_to_angstrom #Decorator for transforming array in nanometer to array in angstroms
    def compute_contacts(self, residue_pairs):

        """

        :param residue_pairs:  An array containing pairs of indices (0-indexed) of residues to compute the contacts between
        :return: distances:  np.ndarray, shape=(n_frames, n_pairs); residues_pairs: np.ndarray, shape=(n_pairs, 2)

        """

        return md.compute_contacts(self.traj, residue_pairs)


    def compute_angles(self, angle_indices):

        """

        :param angle_indices: Each row gives the indices of three atoms which together make an angle (np.ndarray, shape=(num_angles, 3), dtype=int)
        :return: The angles are in radians (np.ndarray, shape=[n_frames, n_angles], dtype=float)

        """

        return md.compute_angles(self.traj, angle_indices)


    def compute_dihedrals(self, indices):

        """

        :param indices: Each row gives the indices of four atoms which together make a dihedral angle (np.ndarray, shape=(n_dihedrals, 4), dtype=int)
        :return: dihedrals : np.ndarray, shape=(n_frames, n_dihedrals), dtype=float. The output array gives,
        in each frame from the trajectory, each of the n_dihedrals torsion angles. The angles are measured in radians

        """

        return md.compute_dihedrals(self.traj, indices)



    def compute_sasa(self, radius=0.14, mode="residue"):


        """
        Compute the solvent accessible surface area of each atom or residue in each simulation frame

        :param radius: The radius of the probe, in nm
        :param mode: In mode == atom the extracted areas are resolved peratom In mode == residue,
        this is consolidated down to the per-residue SASA by summing over the atoms in each residue
        :return: The accessible surface area of each atom or residue in every frame. If mode == atom,
        the second dimension will index the atoms in the trajectory,
        whereas if mode == residue, the second dimension will index the residues. ( np.array, shape=(n_frames, n_features)

        """


        return md.shrake_rupley(self.traj, probe_radius = radius, mode=mode)




    def compute_radius_of_gyration(self):

        """
        Compute the radius of gyration for every frame.
        :return: Rg for every frame (ndarray)

        """

        return md.compute_rg(self.traj)



    def compute_inertia_tensor(self):

        """

        Compute the inertia tensor of a trajectory.
        :return: I_ab: np.ndarray, shape=(traj.n_frames, 3, 3), dtype=float64 (Inertia tensors for each frame)

        """

        return md.compute_inertia_tensor(self.traj)


class Plotter:



    def __init__(self, x_axis, y_axis, figure_name, title, x_label, y_label, style="ggplot", z_axis=None, plot=False, save=True):


        self.x_axis = x_axis
        self.y_axis = y_axis
        self.z_axis = z_axis
        self.plot = plot
        self.save = save
        self.path = "."
        self.dpis = 300
        self.title = title
        self.x_label = x_label
        self.y_label = y_label
        self.z_label=None
        self.cmap = "plasma"
        self.figure_name = figure_name

        plt.style.use(style)



    def scatter_plot(self):

        plt.plot(self.x_axis, self.y_axis)
        plt.title(self.title); plt.xlabel(self.x_label); plt.ylabel(self.y_label)

        if self.z_axis is not None:
            bar = plt.colorbar(); bar.set_label(self.z_label)
        if self.plot: plt.show()
        if self.save: plt.savefig(os.path.join(self.path, "md_{}_plot.png".format(self.figure_name)),dpi=self.dpis)
        plt.clf()


    def box_plot(self):


        sns.boxplot(self.y_axis, orient="v")
        plt.title(self.title);
        plt.xlabel(self.x_label);
        plt.ylabel("Trajectory")
        if self.plot:plt.show()
        if self.save: plt.savefig(os.path.join(self.path, "md_{}_boxplot.png".format(self.figure_name)), dpi=self.dpis)
        plt.clf()

    def density_plot(self):

        sns.distplot(self.y_axis)
        plt.title(self.title);
        plt.xlabel(self.y_label);
        plt.ylabel("Density (1/{})".format(self.y_label.split(" ")[1][1:-1]))
        if self.plot:plt.show()
        if self.save: plt.savefig(os.path.join(self.path, "md_{}_densityplot.png".format(self.figure_name)), dpi=self.dpis)
        plt.clf()

    def superpose_plots(self, Y_array):

        """

        :param Y_array: array of Y arrays for different trajectories (for instance when comparing the effect between
        mutants)
        :return:

        """

        #TODO: we have to load another trajectory so as to be compared

        for array in Y_array:
            plt.plot(self.x_axis, array)


        if self.plot: plt.show()
        if self.save: plt.savefig(self.path, dpi=self.dpis)
        plt.clf()






def parse_args():


    parser = argparse.ArgumentParser()

    parser.add_argument("traj", type=str, help="Trajectory", nargs = '*')
    parser.add_argument("top", type=str, help="Topology file (normally in .gro format)")
    parser.add_argument("-I", "--info", help="Save informfation of trajectory in file", action="store_true")
    parser.add_argument("-R","--rmsd",help="Compute the RMSD between one reference and the trajectory",action="store_true")
    parser.add_argument("-RF","--rmsf",help="Compute the RMSF between one reference and the trajectory",action="store_true")
    parser.add_argument("-LR","--localrmsd",help="Compute the local RMSD between one reference and the trajectory",type=str,nargs='*')
    parser.add_argument("-D","--distance", type=int, help="Two atoms (number) to compute it distance along the trajectory",nargs=2)
    parser.add_argument("-C", "--contact", type=int, help="Two residues (number of residues) to compute its contacts along the trajectory", nargs=2)
    parser.add_argument("-DIS", "--displacement", type=int, help="Atom pair for computing its displacements along the tractory", nargs=2)
    parser.add_argument("-G", "--gyration", help="Compute the gyration at each trajectory frame", action="store_true")
    parser.add_argument("-S", "--sasa",help="Compute the solvent accessible surface area of each atom or residue in each simulation frame (shrake rupley method)", action="store_true")
    parser.add_argument("-PS","--plot_style", type=str, help="Style of the plots (default=ggplot)", default="ggplot")
    parser.add_argument("-P", "--plot", help="Show plots", action="store_true")
    parser.add_argument("-SP", "--save_plot", help="Save plots in directory", action="store_true")
    parser.add_argument("-T","--time", help="Conversion factor from frames to time-scale", type=float, default=1)
    parser.add_argument("-AC", "--acid", help="Take into account for Acid-His distance", action="store_true")

    args = parser.parse_args()

    return args.traj, args.top, args.rmsd, args.rmsf, args.localrmsd, args.distance, args.contact, args.displacement, args.gyration , args.sasa, \
           args.plot_style, args.plot, args.save_plot, args.time, args.acid







def main():

    def average_property(properties):
        length, new_prop = [],[]

        for array in properties:
            length.append(len(array))

        minimum = min(length)

        for array in properties:
            new_prop.append(array[:minimum])

        return np.average(new_prop, axis=0)

    def execute_plots(x_axis,y_axis,xlabel,ylabel,title = "", figure_name="Plot"):
        plot_object = Plotter(x_axis, y_axis, x_label = xlabel, y_label = ylabel, title = title, figure_name = figure_name, plot=plot, save=save_plot)
        plot_object.scatter_plot()
        plot_object.box_plot()
        plot_object.density_plot()

    tra, top, rmsd, rmsf, local_rmsd, distance, contact, displacement, gyration, sasa, plot_style, plot, save_plot, time, acid = parse_args()

    xtc = OpenFiles(tra, top)
    number_of_frames = 10000000
    if ".xtc" in tra[0]:
        trajectory = xtc.load_xtc()
        for traj in trajectory:
            if xtc.number_frames(traj) < number_of_frames:
                number_of_frames = xtc.number_frames(traj)
    else:
        trajectory = xtc.load_trajectory()

    x_axis  = np.arange(0,number_of_frames,1)*time
    Residue_number = [i for i in range(len(trajectory[0].topology.select("name CA")))]

    RMSD, LRMSD, RMSF, Distances, Contacts = [],[],[],[],[]

    for traj in trajectory:

        prop = TrajectoryProperties(traj)

        if rmsd:

            RMSD.append(prop.traj_rmsd(traj,traj.topology.select("backbone")))

        if local_rmsd is not None:

            LRMSD.append(prop.traj_rmsd(traj,traj.topology.select("resid {}".format(" ".join(local_rmsd)))))

        if rmsf:

            RMSF.append(prop.traj_rmsf(traj,traj.topology.select("name CA")))

        if distance is not None:

            if acid:
                D1 = prop.compute_distance([distance])
                D2 = prop.compute_distance([[distance[0],distance[1]+1]])
                distances = []
                for i in range(len(D1)):
                    distances.append(min(D1[i],D2[i]))
                Distances.append(np.array(distances))
            
            else:
                Distances.append(prop.compute_distance([distance]))

        if contact is not None:

            Contacts.append(prop.compute_contacts([contact]))
            pl = Plotter(x_axis, contacts, x_label = "Time (ns)", y_label = "Distance ($\AA$)", figure_name="contact", plot=plot, save=save_plot)

    if rmsd:
        RMSD = average_property(RMSD)
        execute_plots(x_axis,RMSD,"Time (ns)","RMSD (nm)",title = "Global RMSD of the MD simulation",figure_name = "RMSD")

    if local_rmsd is not None:
        LRMSD = average_property(LRMSD)
        execute_plots(x_axis,LRMSD,"Time (ns)","RMSD (nm)",title = "Local RMSD of the MD simulation in residues {}".format(" ".join(local_rmsd)), figure_name = "LocalRMSD_{}".format("_".join(local_rmsd)))

    if rmsf:
        RMSF = average_property(RMSF)
        execute_plots(Residue_number,RMSF,"Residue number","RMSF (nm)",title = "RMSF of the MD simulation",figure_name = "RMSF")

    if distance is not None:
        Distances = average_property(Distances)
        execute_plots(x_axis,Distances,"Time (ns)","Distance ($\AA$)",title = "Distance of the MD simulation between atoms {}".format(" and ".join(str(index) for index in distance)),figure_name = "distance_{}".format("_".join(str(index) for index in distance)))

    if contact is not None:
        Contacts = average_property(Contacts)
        execute_plots(x_axis,Contacts,"Time (ns)","Distance ($\AA$)",title = "Distance of the MD simulation between residues {}".format(" and ".join(str(index) for index in contact)),figure_name = "contact_{}".format("_".join(str(index) for index in contact)))




if __name__=="__main__":
    main()








