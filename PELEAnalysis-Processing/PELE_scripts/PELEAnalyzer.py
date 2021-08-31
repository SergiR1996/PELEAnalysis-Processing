# -*- coding: utf-8 -*-


# Global imports
from __future__ import unicode_literals
import os
import glob
import pickle
import argparse as ap
import pandas as pd
import numpy as n
import multiprocessing as mp

# Plot imports
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("tkagg")
import seaborn as sns

# Local imports
from PELEParseReports import *

# Script information
__author__ = "Sergi Rodà"
__license__ = "MIT"
__version__ = "1.0.2"
__maintainer__ = "Sergi Rodà"
__email__ = "sergi.rodallordes@bsc.es"


class PELEAnalyzer():

    def __init__(self):

        self.reports, self.output_path, self.threshold, self.column_number, self.window_size, self.catalytic_event,\
        self.to_drop, self.n_of_ester_groups, self.add_histidine, self.cysteine, self.threshold_histidine,\
        self.catalytic_distance, self.verbose, self.perform_plots, self.violin_plots,\
        self.num_steps, self.n_processors,self.analysis  = self.parseArgs()

    def parseArgs(self):
        """
        Parse arguments from command-line

        RETURNS
        -------
        reports : string
                  list of report files to look for data
        output_path : string
                      output directory where the csv file will be saved
        threshold: float
                      threshold for the desired metric
        column_number: integer
                      index of the column the desired metric
        window_size: integer
                      Number of steps to consider an entrance
        num_steps: ineger
                      Number of steps per PELE epoch
        n_processors: integer
                      Number of processors
        analysis: string
                      Analysis to perform
        """

        parser = ap.ArgumentParser(description='Script that performs different analysis of the \
            metrics of the reports file from a PELE simulation.')
        optional = parser._action_groups.pop()
        required = parser.add_argument_group('required arguments')
        required.add_argument("-i", "--input", required=True, metavar="FILE",
                              type=str, nargs='*', help="path to report files")
        optional.add_argument("-o", "--output", metavar="PATH", type=str,
                              help="output path to save figure", default="PELE_results")
        optional.add_argument("-TE","--threshold", metavar="THRESHOLD",type=float,
                              help="threshold for the desired metric", default = 4.0)
        optional.add_argument("-C","--column_number", metavar="INTEGER",type=int,
                              help="index of the column where the desired metrics resides", default=7)
        optional.add_argument("-W","--window_size", metavar="INTEGER",type=int,
                              help="number of steps between being out and entering the threshold to consider entrance", default = 2)
        optional.add_argument("-CE","--catalytic_event", metavar="LIST",type=str,
                              nargs='*',help="index of the column where catalytic distances reside (they must be 3)")
        optional.add_argument("-TD","--to_drop", metavar="LIST",type=str,
                              nargs='*',help="column names that want to be dropped", default=[])
        optional.add_argument("-NE","--n_of_ester_groups", metavar="INTEGER",type=int,
                              help="number of ester groups in the substrate used in the PELE simulation", default = 1)
        optional.add_argument("-AH","--add_histidine",
                              help="Add the catalytic histidine - ether O atom(s) of the substrate to the count of catalytic events",
                              action="store_true")
        optional.add_argument("-CS","--cysteine",
                              help="The code takes into account that the nucleophile residue is a cysteine for the catalytic events",
                              action="store_true")
        optional.add_argument("-TH","--threshold_histidine", metavar="THRESHOLD",type=float,
                              help="threshold for the distance between the N atom of the His and the ether O atom of the substrate to consider as catalytic", default = 6.5)
        optional.add_argument("-CD","--catalytic_distance", metavar="THRESHOLD",type=float,
                              help="threshold for the hydrogen bonds of the catalytic residues", default = 3.5)
        optional.add_argument("-V","--verbose",
                              help="Activate the verbose mode", action="store_true")
        optional.add_argument("-PP","--perform_plots",
                              help="Saves the plots of the calculated metrics/parameters", action="store_true")
        optional.add_argument("-VP","--violin_plots",
                              help="Perform violin plots instead of box plots", action="store_true")
        optional.add_argument("-NS","--num_steps", metavar="INTEGER",type=int,
                              help="number of steps per report", default = 40)
        optional.add_argument("-NP","--n_processors", metavar="INTEGER",type=int,
                              help="number of processors to execute the code", default = 4)
        optional.add_argument("-A","--analysis", metavar="STRING",type=str,
                              help="Type of analysis to perform", default="CATE")
        parser._action_groups.append(optional)
        args = parser.parse_args()

        reports = parseReports(args.input, parser)

        if args.catalytic_event is not None:
            args.catalytic_event = [int(i)-1 for i in args.catalytic_event]

        return reports, args.output, args.threshold, int(args.column_number)-1, args.window_size, \
        args.catalytic_event, args.to_drop, args.n_of_ester_groups, args.add_histidine, args.cysteine,\
        args.threshold_histidine, args.catalytic_distance, args.verbose, args.perform_plots,\
        args.violin_plots, args.num_steps, args.n_processors, args.analysis

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

    def Catalytic_events_and_means(self, report):
        """
        Take the PELE simulation report files and obtains the number of steps with 
        a distance (metric) smaller than a threshold and the number of times in the
        simulation the the value of the desired metric gets lower than the threshold

        RETURNS
        -------
        instances: integer
                      Number of inside steps
        entrance: integer
                      Number of entrances        
        """
        values_aux, cat_events, cat_trajectories = [], [], []

        rep = pd.read_csv(report,sep="    ")
        rep.dropna(axis=1,inplace=True)
        values_aux.append(rep.values.tolist())
        for i in range(self.n_of_ester_groups):
            if self.add_histidine:
                if i == 0:
                    if self.cysteine:
                        CATE = rep[(rep[rep.columns[self.catalytic_event[0]+i]] <= self.threshold) & \
                        (rep[rep.columns[self.catalytic_event[0]+i+1]] <= self.threshold_histidine) & \
                        (rep[rep.columns[self.catalytic_event[1]]] <= self.catalytic_distance) & \
                        (rep[rep.columns[self.catalytic_event[2]]] <= self.catalytic_distance)]
                    else:
                        CATE = rep[(rep[rep.columns[self.catalytic_event[0]+i]] <= self.threshold) & \
                        (rep[rep.columns[self.catalytic_event[0]+i+1]] <= self.threshold_histidine) & \
                        (rep[rep.columns[self.catalytic_event[1]]] <= self.catalytic_distance) & \
                        ((rep[rep.columns[self.catalytic_event[2]]] <= self.catalytic_distance) \
                        | (rep[rep.columns[self.catalytic_event[2]+1]] <= self.catalytic_distance))]
                else:
                    if self.cysteine:
                        CATE = rep[(rep[rep.columns[self.catalytic_event[0]+2^i]] <= self.threshold) & \
                        (rep[rep.columns[self.catalytic_event[0]+2**i+1]] <= self.threshold_histidine) & \
                        (rep[rep.columns[self.catalytic_event[1]]] <= self.catalytic_distance) & \
                        (rep[rep.columns[self.catalytic_event[2]]] <= self.catalytic_distance)]
                    else:
                        CATE = rep[(rep[rep.columns[self.catalytic_event[0]+2^i]] <= self.threshold) & \
                        (rep[rep.columns[self.catalytic_event[0]+2**i+1]] <= self.threshold_histidine) & \
                        (rep[rep.columns[self.catalytic_event[1]]] <= self.catalytic_distance) & \
                        ((rep[rep.columns[self.catalytic_event[2]]] <= self.catalytic_distance) \
                        | (rep[rep.columns[self.catalytic_event[2]+1]] <= self.catalytic_distance))]
            else:
                if self.cysteine:
                    CATE = rep[(rep[rep.columns[self.catalytic_event[0]+i]] <= self.threshold) & \
                    (rep[rep.columns[self.catalytic_event[1]]] <= self.catalytic_distance) & \
                    (rep[rep.columns[self.catalytic_event[2]]] <= self.catalytic_distance)]
                else:
                    CATE = rep[(rep[rep.columns[self.catalytic_event[0]+i]] <= self.threshold) & \
                    (rep[rep.columns[self.catalytic_event[1]]] <= self.catalytic_distance) & \
                    ((rep[rep.columns[self.catalytic_event[2]]] <= self.catalytic_distance) \
                    | (rep[rep.columns[self.catalytic_event[2]+1]] <= self.catalytic_distance))]
            CE = CATE.shape[0]
            cat_events.append(CE)
            if CE!=0:
                cat_trajectories.append(1)
                if self.verbose:
                    print("{} --> {}".format(report,list(CATE["numberOfAcceptedPeleSteps"])))
            else:
                cat_trajectories.append(0)

        return values_aux, cat_events, cat_trajectories

    def CalculateFreeEnergy(self, report):
        """
        Take the PELE simulation report files and returns the estimated difference 
        in the free energy of two differentiated states according to a/some metric/s.

        OUTPUT
        ------
        output_path.csv File with the difference in the binding energies and the 
        number of accepted PELE steps of each state.
        """

        G1, G2 = [],[]

        rep = pd.read_csv(report,sep="    ")
        rep.dropna(axis=1,inplace=True)
        for i_row in range(rep.shape[0]):
          if rep.loc[i_row][self.column_number]<=self.threshold:
#          if rep.loc[i_row][self.column_number]<=self.threshold and rep.loc[i_row][5]<=0.2:
            G1.append(rep.loc[i_row][4])
          else:
            G2.append(rep.loc[i_row][4])

        return G1, G2

    def EstimateEnantioselectivity(self, report):
        """
        Take the PELE simulation report files and filters the number of steps with 
        a distance (metric) smaller than a certain threshold

        RETURNS
        -------
        Rep: pandas DataFrame
                      Report file with the steps that have the metric below the threshold  
        """

        rep = pd.read_csv(report, sep="    ")
        rep.dropna(axis=1,inplace=True)
        Rep = rep[rep[rep.columns[self.column_number]] <= self.threshold]

        return [i for i in zip(Rep["Binding Energy"],Rep["dihedral"])]
          
    def Time_of_residence_and_number_of_entrances(self, report):
        """
        Take the PELE simulation report files and obtains the number of steps with 
        a distance (metric) smaller than a threshold and the number of times in the
        simulation the the value of the desired metric gets lower than the threshold

        RETURNS
        -------
        instances: integer
                      Number of inside steps
        entrance: integer
                      Number of entrances        
        """

        inside_bool, instances, entrance = False,[],0

        rep = pd.read_csv(report,sep="    ")
        rep.dropna(axis=1,inplace=True)
        step_in,i_aux,entrance_aux,instances_aux = 0,0,0,0
        for i_row in range(1,rep.shape[0]):
            if rep.loc[i_row][self.column_number] < self.threshold and inside_bool == False and i_row != (rep.shape[0]-1):
            # If the distance goes below the threshold and the previous step was above it, execute this
                if i_row != 1:
                    if (rep.loc[i_row][1] - i_aux) >= self.window_size and i_aux !=0: 
                        entrance_aux+=1
                    elif i_aux == 0 or entrance_aux == 0:
                        entrance_aux+=1
                    instances_aux+=1
                else:
                    instances_aux+=(rep.loc[i_row][1]-rep.loc[i_row-1][1])
                step_in = rep.loc[i_row][1]
                inside_bool = True
            elif rep.loc[i_row][self.column_number] < self.threshold and inside_bool == False and i_row == (rep.shape[0]-1):
            # Else if the distance goes below the threshold and this is the last accepted step, execute this
                instances_aux+=(self.num_steps-rep.loc[i_row][1])
                if (rep.loc[i_row][1] - i_aux) >= self.window_size and i_aux !=0:
                    entrance_aux+=1
                elif i_aux == 0 or entrance_aux == 0:
                    entrance_aux+=1
                instances.append(instances_aux)
            elif rep.loc[i_row][self.column_number] < self.threshold and inside_bool == True and i_row != (rep.shape[0]-1):
            # Else if the distance goes below the threshold but the previous step was already below it, execute this
                instances_aux+=(rep.loc[i_row][1]-rep.loc[i_row-1][1])
                step_in = rep.loc[i_row][1]
            elif rep.loc[i_row][self.column_number] > self.threshold and inside_bool == True:
            # Else if the distance goes above the threshold and the previous step was below it, execute this
                instances_aux+=(rep.loc[i_row][1]-rep.loc[i_row-1][1]-1)
                instances.append(instances_aux)
                instances_aux=0
                inside_bool = False
                i_aux = rep.loc[i_row][1]
            elif rep.loc[i_row][self.column_number] < self.threshold and inside_bool == True and i_row == (rep.shape[0]-1):
            # Else if the distance goes below the threshold but the previous step was already below it and it is the last step, execute this
                instances_aux+=(self.num_steps-step_in)
                instances.append(instances_aux)
            else:
                inside_bool = False
        entrance+=entrance_aux
        inside_bool = False
#        total_steps += rep.loc[rep.shape[0]-1][1]

        return instances, entrance

    def CATE_plot(self, values, column_names):
        """
        Function to perform the box/violin plots of the different metrics of the PELE simulation

        PARAMETERS
        ----------
        values: list of metrics averaged over every report
                List of lists

        OUTPUT
        ------
        Perform the different box/violin plots and stores the values and names as pickle files  
        """
        plot_TE_values, plot_TE_names, plot_IE_values, plot_IE_names\
        ,plot_SA_values, plot_SA_names, plot_values, plot_names = [], [], [], [], [], [], [], []
        
        for value in values:
            plot_TE_values += [value[3]]
            plot_TE_names += [column_names[0]]
            plot_IE_values += [value[4]]
            plot_IE_names += [column_names[1]]
            plot_SA_values += [value[5]]
            plot_SA_names += [column_names[2]]
            plot_values += value[6:]
            plot_names += column_names[3:]

        if self.violin_plots:
            sns.violinplot(plot_TE_names, plot_TE_values, flierprops={"marker":"o", "markerfacecolor":"darkgrey", "markeredgecolor":"k", "markersize":10});plt.savefig(self.output_path+"_TE");plt.close()
            sns.violinplot(plot_IE_names, plot_IE_values, flierprops={"marker":"o", "markerfacecolor":"darkgrey", "markeredgecolor":"k", "markersize":10});plt.savefig(self.output_path+"_IE");plt.close()
            sns.violinplot(plot_SA_names, plot_SA_values, flierprops={"marker":"o", "markerfacecolor":"darkgrey", "markeredgecolor":"k", "markersize":10});plt.savefig(self.output_path+"_SA");plt.close()
            sns.violinplot(plot_names, plot_values, flierprops={"marker":"o", "markerfacecolor":"darkgrey", "markeredgecolor":"k", "markersize":10});plt.savefig(self.output_path+"_DIST");plt.close()
        else:
            sns.boxplot(plot_TE_names, plot_TE_values, flierprops={"marker":"o", "markerfacecolor":"darkgrey", "markeredgecolor":"k", "markersize":10});plt.savefig(self.output_path+"_TE");plt.close()
            sns.boxplot(plot_IE_names, plot_IE_values, flierprops={"marker":"o", "markerfacecolor":"darkgrey", "markeredgecolor":"k", "markersize":10});plt.savefig(self.output_path+"_IE");plt.close()
            sns.boxplot(plot_SA_names, plot_SA_values, flierprops={"marker":"o", "markerfacecolor":"darkgrey", "markeredgecolor":"k", "markersize":10});plt.savefig(self.output_path+"_SA");plt.close()
            sns.boxplot(plot_names, plot_values, flierprops={"marker":"o", "markerfacecolor":"darkgrey", "markeredgecolor":"k", "markersize":10});plt.savefig(self.output_path+"_DIST");plt.close()

        inf_TE = open("{}_TE.pkl".format(self.output_path), "wb"); pickle.dump(plot_TE_values, inf_TE); pickle.dump(plot_TE_names, inf_TE); inf_TE.close()
        inf_IE = open("{}_IE.pkl".format(self.output_path), "wb"); pickle.dump(plot_IE_values, inf_IE); pickle.dump(plot_IE_names, inf_IE); inf_IE.close()
        inf_SA = open("{}_SA.pkl".format(self.output_path), "wb"); pickle.dump(plot_SA_values, inf_SA); pickle.dump(plot_SA_names, inf_SA); inf_SA.close()
        inf_DIST = open("{}_DIST.pkl".format(self.output_path), "wb"); pickle.dump(plot_values, inf_DIST); pickle.dump(plot_names, inf_DIST); inf_DIST.close()

    def CATE(self):
        """
        Function to calculate the average metrics and the catalytic
        events of the PELE simulation

        It is called when this script is the main program called by the interpreter
        """
        Results, results, values, cat_events, cat_trajectories = {}, [], [], [], []

        pool = mp.Pool(self.n_processors)
        results.append(pool.map(self.Catalytic_events_and_means,self.reports))
        pool.terminate()

        results = self.DecompressList(results)

        for elem in results:
            values.append(elem[0][0])
            cat_events.append(elem[1])
            cat_trajectories.append(elem[2])

        values = self.DecompressList(values)

        report = pd.read_csv(self.reports[0],sep="    ")
        report.dropna(axis=1,inplace=True)
        column_names = list(report.columns[3:])

        if self.perform_plots: self.CATE_plot(values,column_names)

        means = n.mean(values,axis=0)[3:]
        std = n.std(values,axis=0)[3:]
        means, std = [round(i,3) for i in means],[round(i,3) for i in std]
        for i in range(self.n_of_ester_groups):
            column_names.append("cat_events_{}".format(i+1));column_names.append("cat_n_trajectories_{}".format(i+1))

        for i,j in zip(n.sum(cat_events,axis=0),n.sum(cat_trajectories,axis=0)):
            means += [i,j]
            std += [0,0]

        for key,item,second_item in zip(column_names,means,std):
            Results[key] = (item, second_item)

        df = pd.DataFrame(Results,index=["Value (Mean)", "Standard_deviation"])
        df.drop(self.to_drop, axis = 1, inplace=True)
        df = df.round(3)
        df.to_csv(self.output_path+".csv")

    def CABE(self):
        """
        Function to calculate the mean interaction energy
        of the system according to a metric with a 
        particular threshold

        It is called when this script is the main program called by the interpreter
        """

        Results, G_values, G1_values, G2_values = {}, [], [], []

        pool = mp.Pool(self.n_processors)
        G_values.append(pool.map(self.CalculateFreeEnergy,self.reports))
        pool.terminate()

        G_values = self.DecompressList(G_values)

        for G_value_list in G_values:
            G1_values += G_value_list[0]
            G2_values += G_value_list[1]

        mean_G1 = n.mean(G1_values,axis=0); std_G1 = n.std(G1_values,axis=0)
        mean_G2 = n.mean(G2_values,axis=0); std_G2 = n.std(G2_values,axis=0)
        state_1 = len(G1_values) ; state_2 = len(G2_values)
        rel1 = 100*state_1/(state_1+state_2); rel2 = 100*state_2/(state_1+state_2)
        
        print(mean_G1-mean_G2,rel1,rel2)

        Results["mean,std, and frequency of state 1"] = (mean_G1, std_G1, state_1, rel1)
        Results["mean,std, and frequency of state 2"] = (mean_G2, std_G2, state_2, rel2)

        if self.perform_plots:
            G_names = ["State 1" for i in G1_values]+["State 2" for i in G2_values]
            G_val_grouped = G1_values + G2_values
            print(len(G_names),len(G_val_grouped))
            if self.violin_plots:
                sns.violinplot(G_names, G_val_grouped, flierprops={"marker":"o", "markerfacecolor":"darkgrey", "markeredgecolor":"k", "markersize":10})
            else:
                sns.boxplot(G_names, G_val_grouped, flierprops={"marker":"o", "markerfacecolor":"darkgrey", "markeredgecolor":"k", "markersize":10})
            plt.savefig(self.output_path+"_CABE");plt.close()
            inf_CABE = open("{}_CABE.pkl".format(self.output_path), "wb"); pickle.dump(G_val_grouped, inf_CABE); pickle.dump(G_names, inf_CABE); inf_CABE.close()


        df = pd.DataFrame(Results, index = ["Mean", "Standard deviation", "Frequency","Relative_abundance"])
        df.to_csv(self.output_path+".csv")

    def EE(self):
        """
        Function to calculate the number of ocurrences of pro-R and pro-S poses of the substrate
        in the active site

        It is called when this script is the main program called by the interpreter
        """
        List_of_reports, R, S, Ambiguous = [], [], [], []

        pool = mp.Pool(self.n_processors)
        List_of_reports.append(pool.map(self.EstimateEnantioselectivity,self.reports))
        pool.terminate()

        List_of_reports = self.DecompressList(List_of_reports)
        List_of_reports = self.DecompressList(List_of_reports)

        List_of_reports.sort(key=lambda x: x[0])

        Top_IE_values = List_of_reports[0:self.window_size]

        for energy, dihedral_value in Top_IE_values:
            if dihedral_value <= -40 and dihedral_value >= -140:
                R.append(dihedral_value)
            elif dihedral_value <= 140 and dihedral_value >= 40:
                S.append(dihedral_value)
            else:
                Ambiguous.append(dihedral_value)

        output_file = open("{}.txt".format(self.output_path),"wt")
        output_file.write("Ratio of R: {}\n".format(100*(len(R)/(len(R)+len(S)))))
        output_file.write("Ratio of S: {}\n".format(100*(len(S)/(len(R)+len(S)))))
        output_file.write("Mean dihedral value of R: {}\n".format(n.mean(R)))
        output_file.write("Mean dihedral value of S: {}\n".format(n.mean(S)))
        output_file.write("Mean dihedral value of ambigous: {}\n".format(n.mean(Ambiguous)))
        output_file.write("Number of ambigous poses: {}\n".format(len(Ambiguous)))

    def TRNE(self):
        """
        Function to calculate the time of residence, the inside steps, and the number of entrances

        It is called when this script is the main program called by the interpreter
        """

        results,instances,entrance = [],[],0

        pool = mp.Pool(self.n_processors)
        results.append(pool.map(self.Time_of_residence_and_number_of_entrances,self.reports))
        pool.terminate()

        results = self.DecompressList(results)

        for elem in results:
            if len(elem[0])!=0:
                instances.append(elem[0])
            entrance += elem[1]

        instances = self.DecompressList(instances)

        if self.perform_plots:
            if self.violin_plots:
                sns.violinplot(["Residence time" for i in instances], instances, flierprops={"marker":"o", "markerfacecolor":"darkgrey", "markeredgecolor":"k", "markersize":10})
            else:
                sns.boxplot(["Residence time" for i in instances], instances, flierprops={"marker":"o", "markerfacecolor":"darkgrey", "markeredgecolor":"k", "markersize":10})
            plt.savefig(self.output_path+"_TR");plt.close()
            inf_NI = open("{}_NI.pkl".format(self.output_path), "wb"); pickle.dump(instances, inf_NI); pickle.dump(["Residence time" for i in instances], inf_NI); inf_NI.close()

        output_file = open("{}.txt".format(self.output_path),"wt")
        output_file.write("Number of inside steps: {}\n".format(n.sum(instances)))
        output_file.write("Inside events: {}\n".format(len(instances)))
        output_file.write("Residence time: {}\n".format(n.mean(instances,axis=0)))
        output_file.write("Relative residence time: {}\n".format((100*n.mean(instances,axis=0))/(self.num_steps)))
        output_file.write("Number of entrances: {}\n".format(entrance))

if __name__ == "__main__":
    """Call the main function"""
    PELEanalyzer = PELEAnalyzer()
    if PELEanalyzer.analysis.upper() == "TRNE":
        if not os.path.exists("TRNE_"+PELEanalyzer.output_path):
            os.mkdir("TRNE_"+PELEanalyzer.output_path)
        PELEanalyzer.TRNE()
        os.system("mv {}*.* TRNE_{}/".format(PELEanalyzer.output_path,PELEanalyzer.output_path))
    if PELEanalyzer.analysis.upper() == "CATE":
        if not os.path.exists("CATE_"+PELEanalyzer.output_path):
            os.mkdir("CATE_"+PELEanalyzer.output_path)
        PELEanalyzer.CATE()
        os.system("mv {}*.* CATE_{}/".format(PELEanalyzer.output_path,PELEanalyzer.output_path))
    if PELEanalyzer.analysis.upper() == "CABE":
        if not os.path.exists("CABE_"+PELEanalyzer.output_path):
            os.mkdir("CABE_"+PELEanalyzer.output_path)
        PELEanalyzer.CABE()
        os.system("mv {}*.* CABE_{}/".format(PELEanalyzer.output_path,PELEanalyzer.output_path))
    if PELEanalyzer.analysis.upper() == "EE":
        if not os.path.exists("EE_"+PELEanalyzer.output_path):
            os.mkdir("EE_"+PELEanalyzer.output_path)
        PELEanalyzer.EE()
        os.system("mv {}*.* EE_{}/".format(PELEanalyzer.output_path,PELEanalyzer.output_path))
