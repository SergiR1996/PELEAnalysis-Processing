# -*- coding: utf-8 -*-

"""

                            MD ANALYSER (for molecular dynamics trajectories)

This script is for analysing the trajectories obtained by any software. It will use the MDtraj module for such
analysis. There will be some flag options which the user can pick on and the files and plots will be saved in a new
created directory

                                    Done by: Sergi Rodà Llordés

"""

# Global imports
import argparse as ap
import pickle

# Local imports 
from MDAnalysisTools import *

# Script information
__author__ = "Sergi Rodà"
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Sergi Rodà"
__email__ = "sergi.rodallordes@bsc.es"

def parseArgs():


    parser = ap.ArgumentParser()

    parser.add_argument("traj", type=str, help="Trajectory file/s", nargs = '*')
    parser.add_argument("top", type=str, help="Topology file")
    parser.add_argument("-I", "--info", help="Save informfation of trajectory in file", action="store_true")
    parser.add_argument("-R","--rmsd",help="Compute the RMSD between one reference and the trajectory",action="store_true")
    parser.add_argument("-RF","--rmsf",help="Compute the RMSF between one reference and the trajectory",action="store_true")
    parser.add_argument("-LR","--localrmsd",help="Compute the local RMSD between one reference and the trajectory",type=str,nargs='*')
    parser.add_argument("-D","--distance", type=int, help="Two atoms (number) to compute it distance along the trajectory",nargs=2)
    parser.add_argument("-CD","--catalytic_distance",help="Compute the distance against the catalytic residues",type=str,nargs='*')
    parser.add_argument("-HIE","--epsilon_protonated",help="Use the ND1 atom name for Ser-His distance", action="store_true")
    parser.add_argument("-C", "--contact", type=int, help="Two residues (number of residues) to compute its contacts along the trajectory", nargs=2)
    parser.add_argument("-DIS", "--displacement", type=int, help="Atom pair for computing its displacements along the tractory", nargs=2)
    parser.add_argument("-G", "--gyration", help="Compute the gyration at each trajectory frame", action="store_true")
    parser.add_argument("-S", "--sasa",help="Compute the solvent accessible surface area of each atom or residue in each simulation frame (shrake rupley method)", action="store_true")
    parser.add_argument("-PS","--plot_style", type=str, help="Style of the plots (default=ggplot)", default="ggplot")
    parser.add_argument("-P", "--plot", help="Show plots", action="store_true")
    parser.add_argument("-SP", "--save_plot", help="Save plots in directory", action="store_true")
    parser.add_argument("-T","--time", help="Conversion factor from frames to time-scale", type=float, default=1)
    parser.add_argument("-AC", "--acid", help="Take into account for Acid-His distance", action="store_true")
    parser.add_argument("-HB", "--hbond", help="Calculate the H-bonds along the trajectories", action="store_true")
    parser.add_argument("-PN", "--pickle_name", help="name of the stored pickle name", type=str, default="")
    parser.add_argument("-AN","--angle", type=int, help="Three atoms (number) to compute the angle between them along the trajectory",nargs=3)

    args = parser.parse_args()

    return args.traj, args.top, args.rmsd, args.rmsf, args.localrmsd, args.distance, args.catalytic_distance, args.epsilon_protonated, args.contact, args.displacement, args.gyration , args.sasa, \
           args.plot_style, args.plot, args.save_plot, args.time, args.acid, args.pickle_name, args.hbond, args.angle



def main():

    def average_property(properties):
        length, new_prop = [],[]

        for array in properties:
            length.append(len(array))

        minimum = min(length)

        for array in properties:
            new_prop.append(array[:minimum])

        return np.average(new_prop, axis=0)

    def save_pickle(pickle_name, property_name, property_y, property_x):
        if pickle_name != "":
            inf = open("{}_{}.pkl".format(property_name, pickle_name), "wb")
            pickle.dump(property_y, inf); pickle.dump(property_x, inf); inf.close()

    def execute_plots(x_axis,y_axis,xlabel,ylabel,title = "", figure_name="Plot"):
        plot_object = Plotter(x_axis, y_axis, x_label = xlabel, y_label = ylabel, title = title, figure_name = figure_name, plot=plot, save=save_plot)
        plot_object.scatter_plot()
        plot_object.box_plot()
        plot_object.density_plot()

    tra, top, rmsd, rmsf, local_rmsd, distance, catalytic_distance, epsilon_protonated, contact, displacement, gyration, sasa, plot_style, plot, save_plot, time, acid, pickle_name, hbond, angle = parseArgs()

    traj = OpenFiles(tra, top)
    number_of_frames = 10000000
    if ".xtc" in tra[0]:
        trajectory = traj.load_xtc()
        for t in trajectory:
            if traj.number_frames(t) < number_of_frames:
                number_of_frames = traj.number_frames(t)
    else:
        trajectory = traj.load_trajectory()
        number_of_frames = traj.number_frames(trajectory)

    x_axis  = np.arange(0,number_of_frames,1)*time
    Residue_number = [i for i in range(len(trajectory[0].topology.select("name CA")))]
    
    RMSD, LRMSD, RMSF, Distances, Contacts, Angles = [],[],[],[],[],[]

    for traj in trajectory:

        prop = TrajectoryProperties(traj)

        if hbond:
            HB.append(prop.compute_hydrogen_bonds(freq = 0.3))

        if rmsd:

            RMSD.append(prop.traj_rmsd(traj,traj.topology.select("backbone")))

        if local_rmsd is not None:

            #LRMSD.append(prop.traj_rmsd(traj,traj.topology.select("resSeq {} and protein".format(" ".join(local_rmsd)))))

            Res_indices=""

            for elem in local_rmsd:
                if elem==local_rmsd[0] and len(local_rmsd)==1:
                    Res_indices+="(resSeq {})".format(elem)
                elif elem==local_rmsd[0] and len(local_rmsd)==2:
                    Res_indices+="(resSeq {} ".format(elem)
                elif elem==local_rmsd[0]:
                    Res_indices+="(resSeq {}".format(elem)
                elif elem==local_rmsd[len(local_rmsd)-1]:
                    Res_indices+="or resSeq {})".format(elem)
                else:
                    Res_indices+=" or resSeq {} ".format(elem)

            LRMSD.append(prop.traj_rmsd(traj,traj.topology.select(Res_indices+" and protein")))

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
        
        if catalytic_distance is not None:

            if acid:
                try:
                    Asp_index_1=int(traj.topology.select("resSeq {} and name OD1 and protein".format(catalytic_distance[2])))
                    Asp_index_2=int(traj.topology.select("resSeq {} and name OD2 and protein".format(catalytic_distance[2])))
                except:
                    Asp_index_1=int(traj.topology.select("resSeq {} and name OE1 and protein".format(catalytic_distance[2])))
                    Asp_index_2=int(traj.topology.select("resSeq {} and name OE2 and protein".format(catalytic_distance[2])))
                try:
                    His_index=int(traj.topology.select("resSeq {} and name HD1 and protein".format(catalytic_distance[1])))
                except:
                    His_index=int(traj.topology.select("resSeq {} and name HE2 and protein".format(catalytic_distance[1])))
                D1 = prop.compute_distance([[Asp_index_1,His_index]])
                D2 = prop.compute_distance([[Asp_index_2,His_index]])
                distances = []
                for i in range(len(D1)):
                    distances.append(min(D1[i],D2[i]))
                Distances.append(np.array(distances))

            else:
                Ser_index=int(traj.topology.select("resSeq {} and name HG and protein".format(catalytic_distance[0])))
                if not epsilon_protonated:
                    His_index=int(traj.topology.select("resSeq {} and name NE2 and protein".format(catalytic_distance[1])))
                else:
                    His_index=int(traj.topology.select("resSeq {} and name ND1 and protein".format(catalytic_distance[1])))
                Distances.append(prop.compute_distance([[Ser_index,His_index]]))

        if angle is not None:

            Angles.append(prop.compute_angles([angle]))

        if contact is not None:

            Contacts.append(prop.compute_contacts([contact]))
            pl = Plotter(x_axis, contacts, x_label = "Time (ns)", y_label = "Distance ($\AA$)", figure_name="contact", plot=plot, save=save_plot)

    if hbond:
        HB_file = open("H_bond_results.txt","w")
        HB_dict = {}
        for elem in HB:
            for subelem in elem:
                if subelem not in HB_dict:
                    HB_dict[subelem] = 1
                else:
                    HB_dict[subelem] += 1
        for key,value in HB_dict.items():
            v = value/4
            HB_file.write("%s : %s\n"%(key,v))

    if rmsd:
        RMSD = average_property(RMSD)
        print("RMSD: "+str(RMSD.mean())+"(+-)"+str(RMSD.std())+"\n")
        save_pickle(pickle_name, "RMSD", RMSD, x_axis)
        execute_plots(x_axis,RMSD,"Time (ns)","RMSD (nm)",title = "Global RMSD of the MD simulation",figure_name = "RMSD")

    if local_rmsd is not None:
        LRMSD = average_property(LRMSD)
        print("Local RMSD: "+str(LRMSD.mean())+"(+-)"+str(LRMSD.std())+"\n")
        save_pickle(pickle_name, "LRMSD", LRMSD, x_axis)
        execute_plots(x_axis,LRMSD,"Time (ns)","RMSD (nm)",title = "Local RMSD of the MD simulation in residues {}".format(" ".join(local_rmsd)), figure_name = "LocalRMSD_{}".format("_".join(local_rmsd)))

    if rmsf:
        RMSF = average_property(RMSF)
        save_pickle(pickle_name, "RMSF", RMSF, Residue_number)
        execute_plots(Residue_number,RMSF,"Residue number","RMSF (nm)",title = "RMSF of the MD simulation",figure_name = "RMSF")

    if distance is not None:
        Distances = average_property(Distances)
        print("Distance: "+str(Distances.mean())+"(+-)"+str(Distances.std())+"\n")
        save_pickle(pickle_name, "Distance", Distances, x_axis)
        execute_plots(x_axis,Distances,"Time (ns)","Distance ($\AA$)",title = "Distance of the MD simulation between atoms {}".format(" and ".join(str(index) for index in distance)),figure_name = "distance_{}".format("_".join(str(index) for index in distance)))

    if catalytic_distance is not None:
        Distances = average_property(Distances)
        print("Distance: "+str(Distances.mean())+"(+-)"+str(Distances.std())+"\n")
        save_pickle(pickle_name, "Catalytic_distance", Distances, x_axis)
        if acid:
            execute_plots(x_axis,Distances,"Time (ns)","Distance ($\AA$)",title = "Acid-His catalytic distance of the MD simulation",figure_name = "distance_{}".format("_".join(str(index) for index in catalytic_distance[1:3])))
        else:
            execute_plots(x_axis,Distances,"Time (ns)","Distance ($\AA$)",title = "Ser-His catalytic distance of the MD simulation",figure_name = "distance_{}".format("_".join(str(index) for index in catalytic_distance[0:2])))

    if contact is not None:
        Contacts = average_property(Contacts)
        save_pickle(pickle_name, "Contact", Contacts, x_axis)
        execute_plots(x_axis,Contacts,"Time (ns)","Distance ($\AA$)",title = "Distance of the MD simulation between residues {}".format(" and ".join(str(index) for index in contact)),figure_name = "contact_{}".format("_".join(str(index) for index in contact)))

    if angle is not None:        
        Angles = average_property(Angles)
        print("Angle: "+str(Angles.mean())+"(+-)"+str(Angles.std())+"\n")
        save_pickle(pickle_name, "Angle", Angles, x_axis)
        execute_plots(x_axis,Angles,"Time (ns)","Angle (deg)",title = "Angle of the MD simulation between atoms {}".format(" and ".join(str(index) for index in angle)),figure_name = "angle_{}".format("_".join(str(index) for index in angle)))



if __name__=="__main__":
    main()
