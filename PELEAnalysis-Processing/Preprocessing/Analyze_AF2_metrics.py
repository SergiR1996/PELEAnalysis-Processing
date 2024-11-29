import pickle
import numpy as np
import pandas as pd
import sys, glob, os
import seaborn as sns
import operator as op
import matplotlib.pyplot as plt

def plot_results(information, label, Output):
        sorted_keys, sorted_vals = zip(*sorted(information.items(), key=op.itemgetter(0)))
        SK = []
        for ke, el in zip(sorted_keys, sorted_vals):
            #print(ke, el)
            SK += [ke]*len(el)
        SV = [it for sl in sorted_vals for it in sl]
        data = {"Variant": SK, label: SV}
        df = pd.DataFrame(data)
        order = df.groupby('Variant')[label].median().sort_values(ascending=False).index

        plt.figure(figsize=(19, 10.5))
        sns.boxplot(x="Variant",y=label, data=df, order=order)
        plt.xticks(rotation = 90)
        plt.ylabel(label, fontsize=22)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.tight_layout()
        plt.savefig("{}.png".format(Output),dpi=300)

def generate_scoredict(confidences, binderlen) -> None:
        '''
        Collect the confidence values, slicing them to the binder and target regions
        then add the parsed scores to the score_dict
        '''

        plddt_array = confidences['plddt']
        plddt = np.mean(plddt_array)

        plddt_binder = np.mean(plddt_array[:binderlen])
        plddt_target = np.mean(plddt_array[binderlen:])

        pae = confidences['predicted_aligned_error']
 
        pae_interaction1 = np.mean(pae[:binderlen,binderlen:])
        pae_interaction2 = np.mean(pae[binderlen:,:binderlen])
        pae_binder = np.mean(pae[:binderlen,:binderlen])
        pae_target = np.mean(pae[binderlen:,binderlen:])
        pae_interaction_total = (pae_interaction1+pae_interaction2)/2

        score_dict = {
                "plddt_total" : plddt,
                "plddt_binder" : plddt_binder,
                "plddt_target" : plddt_target,
                "pae_binder" : pae_binder,
                "pae_target" : pae_target,
                "pae_interaction_binder": pae_interaction1,
                "pae_interaction_target": pae_interaction2,
                "pae_interaction" : pae_interaction_total,
                "iptm": float(confidences['iptm'])
        }

        return score_dict

ene_files = [i for i in os.listdir("preds") if os.path.isdir(os.path.join("preds",i))]

pLDDT_total, pLDDT_antigen, pLDDT_antibody, pae_antigen, pae_antibody, pae_interaction_antigen, pae_interaction_antibody, pae_interaction, iptm = {}, {}, {}, {}, {}, {}, {}, {}, {}

for ene_file in ene_files:
     pLDDT_total[ene_file], pLDDT_antigen[ene_file], pLDDT_antibody[ene_file], pae_antigen[ene_file], pae_antibody[ene_file], pae_interaction_antigen[ene_file], pae_interaction_antibody[ene_file], pae_interaction[ene_file], iptm[ene_file] = [], [], [], [], [], [], [], [], []
     for j in glob.glob(f"preds/{ene_file}/result*pkl"):
        pickle_file = open(j,"rb")
        pickle_data = pickle.load(pickle_file)
        Scores_dict = generate_scoredict(pickle_data, 389) # 231 is the index where the antibody finishes and the antigen starts

        pLDDT_total[ene_file].append(Scores_dict["plddt_total"])
        pLDDT_antigen[ene_file].append(Scores_dict["plddt_target"])
        pLDDT_antibody[ene_file].append(Scores_dict["plddt_binder"])
        pae_antigen[ene_file].append(Scores_dict["pae_target"])
        pae_antibody[ene_file].append(Scores_dict["pae_binder"])
        pae_interaction_antigen[ene_file].append(Scores_dict["pae_interaction_target"])
        pae_interaction_antibody[ene_file].append(Scores_dict["pae_interaction_binder"])
        pae_interaction[ene_file].append(Scores_dict["pae_interaction"])
        iptm[ene_file].append(Scores_dict["iptm"])

plot_results(pLDDT_total,"pLDDT",f"pLDDT_total_{sys.argv[1]}")
plot_results(pLDDT_antigen,"pLDDT",f"pLDDT_antigen_{sys.argv[1]}")
plot_results(pLDDT_antibody,"pLDDT",f"pLDDT_antibody_{sys.argv[1]}")
plot_results(pae_antigen,"pAE",f"pae_antigen_{sys.argv[1]}")
plot_results(pae_antibody,"pAE",f"pae_antibody_{sys.argv[1]}")
plot_results(pae_interaction_antigen,"pAE_interface",f"pae_interface_antigen_{sys.argv[1]}")
plot_results(pae_interaction_antibody,"pAE_interface",f"pae_interface_antibody_{sys.argv[1]}")
plot_results(pae_interaction,"pAE_interface",f"pae_interface_{sys.argv[1]}")
plot_results(iptm,"ipTM",f"ipTM_{sys.argv[1]}")

df = pd.DataFrame(pLDDT_total.items(), columns=['Variant', 'pLDDT value'])
df_antigen = pd.DataFrame(pLDDT_antigen.items(), columns=['Variant', 'pLDDT value'])
df_antibody = pd.DataFrame(pLDDT_antibody.items(), columns=['Variant', 'pLDDT value'])
df_pae_antigen = pd.DataFrame(pae_antigen.items(), columns=['Variant', 'pAE value'])
df_pae_antibody = pd.DataFrame(pae_antibody.items(), columns=['Variant', 'pAE value'])
df_pae_interaction_antigen = pd.DataFrame(pae_interaction_antigen.items(), columns=['Variant', 'pAE_interface value'])
df_pae_interaction_antibody = pd.DataFrame(pae_interaction_antibody.items(), columns=['Variant', 'pAE_interface value'])
df_pae_interaction = pd.DataFrame(pae_interaction.items(), columns=['Variant', 'pAE_interface value'])
df_iptm = pd.DataFrame(iptm.items(), columns=['Variant', 'ipTM value'])

df.to_csv(f"{sys.argv[1]}_total.csv")
df_antigen.to_csv(f"{sys.argv[1]}_antigen.csv")
df_antibody.to_csv(f"{sys.argv[1]}_antibody.csv")
df_pae_antigen.to_csv(f"{sys.argv[1]}_pae_antigen.csv")
df_pae_antibody.to_csv(f"{sys.argv[1]}_pae_antibody.csv")
df_pae_interaction_antigen.to_csv(f"{sys.argv[1]}_pae_interaction_antigen.csv")
df_pae_interaction_antibody.to_csv(f"{sys.argv[1]}_pae_interaction_antibody.csv")
df_pae_interaction.to_csv(f"{sys.argv[1]}_pae_interaction.csv")
df_iptm.to_csv(f"{sys.argv[1]}_iptm.csv")
