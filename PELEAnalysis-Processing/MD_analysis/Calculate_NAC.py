# Global imports
import mdtraj as md
import numpy as np
import multiprocessing as mp
import pandas as pd
import argparse as ap

# Local imports
from MDAnalysisTools import *

# Script information
__author__ = "Sergi Rodà"
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Sergi Rodà"
__email__ = "sergi.roda@nostrumbiodiscovery.com"

def parseArgs():

    parser = ap.ArgumentParser()

    parser.add_argument("traj", type=str, help="Trajectory file/s", nargs = '*')
    parser.add_argument("top", type=str, help="Topology file")
    parser.add_argument("-NP","--n_processors", metavar="INTEGER",type=int,
                              help="number of processors to execute the code", default = 4)
    parser.add_argument("-o","--output", metavar="STRING", type=str,
                              help="filename for the output trajectory", default="analysis")
    parser.add_argument("-la","--lig_atom", metavar="STRING", type=str, nargs="*",
                              help="atom name of the ligand molecule", default="C1")
    parser.add_argument("-la2","--lig_atom2", metavar="STRING", type=str, nargs="*",
                              help="atom name of the ligand molecule", default="O1")
    parser.add_argument("-S","--serine", metavar="INTEGER", type=int,
                              help="residue number of the catalytic serine residue")
    parser.add_argument("-H","--histidine", metavar="INTEGER", type=int,
                              help="residue number of the catalytic histidine residue")
    parser.add_argument("-A","--acid", metavar="INTEGER", type=int,
                              help="residue number of the catalytic acid residue")
    parser.add_argument("-T","--threshold", metavar="FLOAT", type=float,
                              help="Threshold of the serine-substrate distance", default=3.5)
    parser.add_argument("-TH","--threshold_histidine", metavar="FLOAT", type=float,
                              help="Threshold of the histidine-substrate distance", default=3.5)
    parser.add_argument("-CD","--catalytic_distance", metavar="FLOAT", type=float,
                              help="Threshold of the catalytic distances", default=3.5)

    args = parser.parse_args()

    return args.traj, args.top, args.n_processors, args.output, args.lig_atom, args.lig_atom2, args.serine, args.histidine, args.acid, args.catalytic_distance, args.threshold_histidine, args.threshold

def calculateNAC(t):

    prop = TrajectoryProperties(t)

    if len(t.topology.select(f"resSeq {his_res} and name HD1"))==0:
        his_atom = "ND1"
        his_atom_h = "HE2"
    else:
        his_atom = "NE2"
        his_atom_h = "HD1"

    if len(t.top.select(f"resSeq {acid_res} and name OD1"))==0:
        acid_atom1 = "OE1"
        acid_atom2 = "OE2"
    else:
        acid_atom1 = "OD1"
        acid_atom2 = "OD2"

    if len(lig_atom) > 1:
        ser_sub_distances = []
        for la in lig_atom:
            NAC11 = int(t.topology.select(f"name {la}"))
            NAC12 = int(t.topology.select(f"resSeq {ser_res} and name OG"))
            ser_sub_distances.append(prop.compute_distance([[NAC11,NAC12]]))
        ser_sub_distances = np.concatenate(ser_sub_distances, axis=1)
        ser_sub_distances = np.min(ser_sub_distances, axis=1)
    else:
        NAC11 = int(t.topology.select(f"name {lig_atom[0]}"))
        NAC12 = int(t.topology.select(f"resSeq {ser_res} and name OG"))

    if len(lig_atom2) > 1:
        his_sub_distances = []
        for la2 in lig_atom2:
            NAC21 = int(t.topology.select(f"name {la2}"))
            NAC22 = int(t.topology.select(f"resSeq {his_res} and name {his_atom}"))
            his_sub_distances.append(prop.compute_distance([[NAC21,NAC22]]))
        his_sub_distances = np.concatenate(his_sub_distances, axis=1)
        his_sub_distances = np.min(his_sub_distances, axis=1)
    else:
        NAC21 = int(t.topology.select(f"name {lig_atom2[0]}"))
        NAC22 = int(t.topology.select(f"resSeq {his_res} and name {his_atom}"))

    NAC31 = int(t.topology.select(f"resSeq {ser_res} and name HG"))
    NAC32 = int(t.topology.select(f"resSeq {his_res} and name {his_atom}"))
    NAC41 = int(t.topology.select(f"resSeq {his_res} and name {his_atom_h}"))
    NAC421 = int(t.topology.select(f"resSeq {acid_res} and name {acid_atom1}"))
    NAC422 = int(t.topology.select(f"resSeq {acid_res} and name {acid_atom2}"))

    if (len(lig_atom2) <= 1) and (len(lig_atom2) <= 1):
        ser_sub_distances = np.array([el[0] for el in prop.compute_distance([[NAC11,NAC12]])])
        his_sub_distances = np.array([el[0] for el in prop.compute_distance([[NAC21,NAC22]])])
        ser_his_distances = np.array([el[0] for el in prop.compute_distance([[NAC31,NAC32]])])
        his_acid_distances = prop.compute_distance([[NAC41,NAC421], [NAC41,NAC422]])
        his_acid_distances = np.min(his_acid_distances, axis=1)
    elif (len(lig_atom2) <= 1) and (len(lig_atom2) > 1):
        raise Exception("If multiple substrate atoms are provided for the serine-substrate distance a same amount of atoms must be provided for histidine-substrate distance")
    elif (len(lig_atom2) > 1) and (len(lig_atom2) <= 1):
        raise Exception("If multiple substrate atoms are provided for the serine-substrate distance a same amount of atoms must be provided for histidine-substrate distance")
    else:
        ser_his_distances = np.array([el[0] for el in prop.compute_distance([[NAC31,NAC32]])])
        his_acid_distances = prop.compute_distance([[NAC41,NAC421], [NAC41,NAC422]])
        his_acid_distances = np.min(his_acid_distances, axis=1)

    Distances = np.stack((ser_sub_distances, his_sub_distances, ser_his_distances, his_acid_distances), axis=1)

    Distances_df = pd.DataFrame(Distances, columns = ["NAC1","NAC2","NAC3","NAC4"])

    Thresholds = {"NAC1": threshold,"NAC2": threshold_histidine,"NAC3": catalytic_distance, "NAC4": catalytic_distance}

    Slice_boolean =(Distances_df["NAC1"] <= Thresholds["NAC1"]) & (Distances_df["NAC2"] <= Thresholds["NAC2"]) & (Distances_df["NAC3"] <= Thresholds["NAC3"])

    NAC_t = t.slice(Slice_boolean)

    Dist_df = Distances_df.copy()
    for T in Thresholds:
        Dist_df = Dist_df[Dist_df[T] <= Thresholds[T]]

    return Dist_df.shape[0]/Distances_df.shape[0], NAC_t

def main():

    global lig_atom; global lig_atom2; global ser_res; global his_res; global acid_res; global catalytic_distance; global threshold_histidine; global threshold
    traj, top, n_processors, output, lig_atom, lig_atom2, ser_res, his_res, acid_res, catalytic_distance, threshold_histidine, threshold  = parseArgs()
    Traj = OpenFiles(traj, top)
    trajectories = Traj.load_xtc()

    NAC, NAC_t, NAC_all = [], [], []
    pool = mp.Pool(n_processors)
    NAC_all += pool.map(calculateNAC, trajectories)
    pool.terminate()

    for el,li in NAC_all:
        NAC.append(el)
        NAC_t.append(li)

    print(NAC)
    print(np.mean(NAC))
    NAC_t_all = md.join(NAC_t)
    NAC_t_all.save(f"NAC_{output}.xtc")


if __name__=="__main__":
    main()
