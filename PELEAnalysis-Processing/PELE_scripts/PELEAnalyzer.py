# -*- coding: utf-8 -*-


# Global imports
from __future__ import unicode_literals
import os
import glob, json
import pickle
import argparse as ap
import pandas as pd
import numpy as n
import multiprocessing as mp

# Plot imports
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("agg")
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

        self.reports, self.output_path, self.window_size, self.to_drop,\
        self.catalytic_dict, self.verbose, self.perform_plots, self.violin_plots,\
        self.num_steps, self.n_processors, self.KT,self.separation  = self.parseArgs()

    def parseArgs(self):
        """
        Parse arguments from command-line

        RETURNS
        -------
        reports : string
                  list of report files to look for data
        output_path : string
                      output directory where the csv file will be saved
        catalytic_json: string
                      Dictionary with the key being the column name of the metrics
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
        optional.add_argument("-W","--window_size", metavar="INTEGER",type=int,
                              help="number of steps between being out and entering the threshold to consider entrance", default = 2)
        optional.add_argument("-TD","--to_drop", metavar="LIST",type=str,
                              nargs='*',help="column names that want to be dropped", default=[])
        optional.add_argument("-CJ","--catalytic_json", metavar="DICTIONARY [FILE]",type=str,
                                help="Dictionary with the key being the column name of the metrics that define a catalytic pose and the value the threshold", default=None)
        optional.add_argument("-V","--verbose",
                              help="Activate the verbose mode", action="store_true")
        optional.add_argument("-PP","--perform_plots",
                              help="Saves the plots of the calculated metrics/parameters", action="store_true")
        optional.add_argument("-VP","--violin_plots",
                              help="Perform violin plots instead of box plots", action="store_true")
        optional.add_argument("-NS","--num_steps", metavar="INTEGER",type=int,
                              help="number of steps per report", default = 0)
        optional.add_argument("-NP","--n_processors", metavar="INTEGER",type=int,
                              help="number of processors to execute the code", default = 4)
        optional.add_argument("-KT","--Boltzmann_constant", metavar="FLOAT",type=float,
                              help="Boltzmann constant and temperature value for the probabilities", default = 0.593)
        optional.add_argument("-S","--separation", metavar="STRING", type=str,
                              help="Type of delimiter between columns of the report file", default="    ")
        parser._action_groups.append(optional)
        args = parser.parse_args()

        reports = parseReports(args.input, parser)

        with open(args.catalytic_json) as json_data:
            catalytic_dict = json.load(json_data)

        return reports, args.output, args.window_size, \
        args.to_drop, catalytic_dict, args.verbose, args.perform_plots,\
        args.violin_plots, args.num_steps, args.n_processors, args.Boltzmann_constant, args.separation

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

    def Calculate_total_catalytic_events(self, set_of_thresholds, report, num_steps):
        """
        Helper function to calculate the catalytic events out of the total PELE steps
        whether are accepted or rejected

        PARAMETERS
        ----------
        set_of_thresholds: list of booleans
                List of booleans

        RETURNS
        -------
        Total_CE: The total number of catalytic events in the report file
                Integer/Float
        """

        Total_CE = 0

        for i,row in enumerate(set_of_thresholds):
            if i>0:
                if row==True and i != len(set_of_thresholds)-1 and set_of_thresholds[i-1]==True:
                    Total_CE += (report.loc[i][1]-report.loc[i-1][1])
                elif row==False and set_of_thresholds[i-1]==True:
                    Total_CE += (report.loc[i][1]-report.loc[i-1][1])
                elif row==True and i == len(set_of_thresholds)-1:
                    Total_CE += (num_steps-report.loc[i][1])
                else:
                    continue
            else:
                continue

        return  Total_CE

    def filter_dataframe(self, df):
        """
        Filters the DataFrame based on dynamic conditions provided as a dictionary.
        
        Args:
            df (pd.DataFrame): The DataFrame to filter.
            conditions (dict): A dictionary where each key is a column name and the value is either:
                - A value to filter for equality, or
                - A tuple (operator, value) for more complex comparisons. 
                Supported operators include '==', '!=', '>', '<', '>=', '<='.
        
        Returns:
            pd.DataFrame: The filtered DataFrame.
        """
        query_list = []

        filtered_boolean = df["#Task"]==1
        for col, cond in conditions.items():
            filtered_boolean = filtered_boolean & (df[col] <= cond)
            # Check if the condition is a tuple with an operator and a value.
            operator = "<="
            value = cond
            value_str = str(value)
            
            # Use backticks around the column name in case it contains special characters.
            query_list.append(f"`{col}` {operator} {value_str}")
        
        # Join all condition strings with 'and'
        query_str = " and ".join(query_list)
        
        # Use the query method to filter the DataFrame
        filtered_df = df.query(query_str)
        return filtered_df, filtered_boolean

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

        values_aux, cat_events, cat_trajectories, total_catalytic_events, total_num_steps = [], [], [],  [], 0

        rep = pd.read_csv(report,sep=self.separation)
        logfilename = "/".join(report.split("/")[:-1])+"/logFile_"+report.split("/")[-1].split("metric")[0][-3:].replace("_","").replace("t","")+".txt"
        if os.path.isfile(logfilename):
            logfile = open(logfilename, "rt")
            for line in logfile:
                if "New Step" in line:
                    num_steps = int(line.split()[-1])
        if self.num_steps != 0:
            num_steps = self.num_steps
        total_num_steps+=num_steps
        rep.dropna(axis=1,inplace=True)
        values_aux.append(rep.values.tolist())

        catalytic_series, catalytic_boolean = self.filter_dataframe(rep)
        Total_CE = self.Calculate_total_catalytic_events(catalytic_boolean, rep, num_steps)
        CE = catalytic_series.shape[0]
        cat_events.append(CE)
        total_catalytic_events.append(Total_CE)
        if CE!=0:
            cat_trajectories.append(1)
            if self.verbose:
                print(f"{report} --> {list(catalytic_series['numberOfAcceptedPeleSteps'])}")
        else:
            cat_trajectories.append(0)
        catalytic_indexes = catalytic_series.index.values.tolist()
        non_catalytic_series = rep.loc[~rep.index.isin(catalytic_indexes)]

        return values_aux, cat_events, cat_trajectories, rep.shape[0], total_catalytic_events, total_num_steps, list(rep["BindingEnergy"]), list(catalytic_series["BindingEnergy"]), list(non_catalytic_series["BindingEnergy"]), list(rep["currentEnergy"]), list(catalytic_series["currentEnergy"]), list(non_catalytic_series["currentEnergy"])

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

        inf_TE = open(f"{self.output_path}_TE.pkl", "wb"); pickle.dump(plot_TE_values, inf_TE); pickle.dump(plot_TE_names, inf_TE); inf_TE.close()
        inf_IE = open(f"{self.output_path}_IE.pkl", "wb"); pickle.dump(plot_IE_values, inf_IE); pickle.dump(plot_IE_names, inf_IE); inf_IE.close()
        inf_SA = open(f"{self.output_path}_SA.pkl", "wb"); pickle.dump(plot_SA_values, inf_SA); pickle.dump(plot_SA_names, inf_SA); inf_SA.close()
        inf_DIST = open(f"{self.output_path}_DIST.pkl", "wb"); pickle.dump(plot_values, inf_DIST); pickle.dump(plot_names, inf_DIST); inf_DIST.close()

    def CATE(self):
        """
        Function to calculate the average metrics and the catalytic
        events of the PELE simulation

        It is called when this script is the main program called by the interpreter
        """
        Results, results, values, cat_events, cat_trajectories, total_accepted_steps, total_cat_events, total_num_steps, BindE_total, BindE_cat, BindE_ncat, TotE_total, TotE_cat, TotE_ncat = {}, [], [], [], [], [], [], 0, [], [], [], [], [], []

        logfilename = "/".join(self.reports[0].split("/")[:-1])+"/logFile_"+self.reports[0].split("/")[-1].split("metric")[0][-3:].replace("_","").replace("t","")+".txt"
        
        if not os.path.isfile(logfilename) and self.num_steps == 0:
            raise ValueError("Logfiles are missing, use the num_of_steps flag instead")

        pool = mp.Pool(self.n_processors)
        results.append(pool.map(self.Catalytic_events_and_means,self.reports))
        pool.terminate()

        results = self.DecompressList(results)

        for elem in results:
            values.append(elem[0][0])
            cat_events.append(elem[1])
            cat_trajectories.append(elem[2])
            total_accepted_steps.append(elem[3])
            total_cat_events.append(elem[4])
            total_num_steps += elem[5]
            BindE_total.append(elem[6])
            BindE_cat.append(elem[7])
            BindE_ncat.append(elem[8])
            TotE_total.append(elem[9])
            TotE_cat.append(elem[10])
            TotE_ncat.append(elem[11])

        values = self.DecompressList(values)
        BindE_total = n.array(self.DecompressList(BindE_total))
        BindE_cat = n.array(self.DecompressList(BindE_cat))
        BindE_ncat = n.array(self.DecompressList(BindE_ncat))
        TotE_total = n.array(self.DecompressList(TotE_total))
        TotE_cat = n.array(self.DecompressList(TotE_cat))
        TotE_ncat = n.array(self.DecompressList(TotE_ncat))

        report = pd.read_csv(self.reports[0],sep=self.separation)
        report.dropna(axis=1,inplace=True)
        column_names = list(report.columns[3:])

        if self.perform_plots: self.CATE_plot(values,column_names)

        # Calculate the Ebc
        energy_minimum = n.min(TotE_total)
        relative_energy = TotE_total-energy_minimum
        Z = n.sum(n.exp(-relative_energy/self.KT))

        relative_energy_cat = TotE_cat-energy_minimum
        probability_cat = n.exp(-relative_energy_cat/self.KT)/Z
        Ebc = n.sum(probability_cat*BindE_cat)

        relative_energy = TotE_ncat-energy_minimum
        probability_ncat = n.exp(-relative_energy/self.KT)/Z
        Ebnc = n.sum(probability_ncat*BindE_ncat)

        dEbc_dEbnc = Ebc-Ebnc

        means = n.mean(values,axis=0)[3:]
        std = n.std(values,axis=0)[3:]
        means, std = [round(i,3) for i in means],[round(i,3) for i in std]
        column_names.append("cat_events");column_names.append("cat_n_trajectories")
        column_names.append("cat_events (%)");column_names.append("total_cat_events")
        column_names.append("total_cat_events (%)")
        column_names.append("Catalytic Free Binding Energy");column_names.append("Non-Catalytic Free Binding Energy")
        column_names.append("Difference in Free Binding Energy")
        means += [n.sum(cat_events),n.sum(cat_trajectories),100*n.sum(cat_events)/n.sum(total_accepted_steps),n.sum(total_cat_events),100*n.sum(total_cat_events)/(total_num_steps),Ebc,Ebnc,dEbc_dEbnc]
        std += ["-","-","-","-","-","-","-","-"]

        for key,item,second_item in zip(column_names,means,std):
            Results[key] = (item, second_item)

        df = pd.DataFrame(Results,index=["Value (Mean)", "Standard_deviation"])
        df.drop(self.to_drop, axis = 1, inplace=True)
        df = df.round(3)
        df.to_csv(self.output_path+".csv")

        output_file = open(f"{self.output_path}_catalytic_events.txt", "wt")
        output_file.write(f"\nNumber of accepted catalytic events in all groups: {n.sum(cat_events)}\n")
        output_file.write(f"Relative frequency of accepted catalytic events in all groups: {round(100*n.sum(cat_events)/n.sum(total_accepted_steps),3)} %\n")
        output_file.write(f"Number of total catalytic events in all groups: {n.sum(total_cat_events)}\n")
        output_file.write(f"Relative frequency of total catalytic events in all groups: {round(100*n.sum(total_cat_events) / (total_num_steps), 3)} %\n")

if __name__ == "__main__":
    """Call the main function"""
    PELEanalyzer = PELEAnalyzer()
    if not os.path.exists(PELEanalyzer.output_path):
        os.mkdir(PELEanalyzer.output_path)
        PELEanalyzer.CATE()
        os.system(f"mv {PELEanalyzer.output_path}*.* {PELEanalyzer.output_path}/")