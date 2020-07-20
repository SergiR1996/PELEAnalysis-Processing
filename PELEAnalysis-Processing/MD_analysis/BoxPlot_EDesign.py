# Global imports
import pickle
import numpy as np
import os,sys
import argparse as ap
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns 
import pandas as pd
from random import randint
plt.rcParams.update({"font.size": 16, "mathtext.default":"regular"})

def parseArgs():
    """
    Parse arguments from command-line

    RETURNS
    -------
    args.input: list
            list of pickle files
    """

    parser = ap.ArgumentParser(description = "Script used to plot simultaneously several MD simulations with the pickle files")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    required.add_argument("-i", "--input", required = True, metavar = "FILE",
                          type = str, nargs = '*', help = "list of pickle files")
    optional.add_argument("-P", "--palette", metavar = "LIST",
                          type = str, nargs = '*', 
                          help = "list of colors for the palette in the plot",
                          default = [])   
    parser._action_groups.append(optional)
    args = parser.parse_args()

    return args.input, args.palette

def preparepickle(Pickle):
    """
    Take the data from the pickle file and add the name that will be added to the boxplot

    RETURNS
    -------
    MD: list
            list of metric values
    MD_names: list
            list of metric names
    """

    inf = open(Pickle, "rb")
    metric_MD = pickle.load(inf); inf.close()
    MD = metric_MD
    if "istance" in Pickle:
        MD = metric_MD[:,0]
        MD_names = ["/".join(Pickle.split("/")[-1].split("_")[2:]).split(".")[0]+Pickle.split("/")[1].split("_")[-1] for elem in MD] # Pickle.split("/")[1].split("_")[-1] is specific for EDesign
        if "Acid" in Pickle:
                MD_names = ["/".join(Pickle.split("/")[-1].split("_")[3:]).split(".")[0]+Pickle.split("/")[1].split("_")[-1] for elem in MD]
    elif "MSD" in Pickle:
        MD = MD[1:]
        MD_names = ["/".join(Pickle.split("/")[-1].split("_")[1:]).split(".")[0]+Pickle.split("/")[1].split("_")[-1] for elem in MD]
    
    return MD, MD_names
    
def Boxplot(Metric, Names, Proof, Proof_Num, Palette):
    """
    Plot the data in the form a boxplot

    OUTPUT
    ------

    It displays the boxplot of all the MD simulations and the metric
    """

    if Palette == []:
        # Create the color palette for the different MD simulations randomly
        for i in range(Proof_Num):
            pal_aux = "#"
            for j in range(6):
                    value = randint(0,15)
                    if value == 10:
                            value = "a"
                    elif value == 11:
                            value = "b"
                    elif value == 12:
                            value = "c"
                    elif value == 13:
                            value = "d"
                    elif value == 14:
                            value = "e"
                    elif value == 15:
                            value = "f"
                    pal_aux += str(value)
            Palette.append(pal_aux)
        print("'" + "' '".join(Palette) + "'")

    # Create the boxplot with the Seaborn library
    _, ax = plt.subplots()
    sns.boxplot(Names, Metric, palette=Palette, ax=ax, 
        linewidth=2.5, 
        flierprops={"marker":"o", "markerfacecolor":"darkgrey", "markeredgecolor":"k", "markersize":10}); 
    ax.set_xticklabels(ax.get_xticklabels(),rotation=90)

    # Add the Y label depending on the type of represented metric
    if "istance" in Proof:
        ax.set_ylabel("Serine-histidine distance ($\AA$)", labelpad=20)
        if "Acid" in Proof:
                ax.set_ylabel("Aspartic-histidine distance ($\AA$)", labelpad=20)
    elif "LRMSD" in Proof:
        ax.set_ylabel("Local RMSD of the catalytic triad (nm)", labelpad=20)
    elif "RMSD" in Proof and not "LRMSD" in Proof:
        ax.set_ylabel("Global RMSD (nm)", labelpad=20)

    # Display the plot, and the different options
    plt.setp(ax.artists, edgecolor="k")
    plt.setp(ax.lines, color="k")
    plt.show()

def main():
    """
    Main function

    It is called when this script is the main program called by the interpreter
    """

    Pickle_list, Palette = parseArgs()

    Names, Metric, Proof_Num = [],[],0

    for Pickle in Pickle_list:
        Aux = preparepickle(Pickle)
        Metric += list(Aux[0])
        Names += Aux[1]
        Proof_Num += 1
        if Pickle == Pickle_list[-1]:
                Proof = Pickle
    
    Boxplot(Metric, Names, Proof, Proof_Num, Palette)

if __name__ == "__main__":
    """Call the main function"""
    main()
