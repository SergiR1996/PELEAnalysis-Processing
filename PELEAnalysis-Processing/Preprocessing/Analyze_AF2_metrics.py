import pickle
import numpy as np
import math, itertools
import pandas as pd
import mdtraj as md
import sys, glob, os
import seaborn as sns
import operator as op
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
import matplotlib.ticker as ticker
from Bio import PDB
import multiprocessing as mp

FAVORED_REGIONS = [
    (-180, -40, -150, 150),  # ϕ: -180° to -40° and ψ: -150° to 150°
    (-120, -40, -150, 150)   # Another set of favored regions (can be expanded for more residues)
]

ALLOWED_REGIONS = [
    (-180, 180, -150, 150),  # A wider range for allowed regions
    (-120, 120, -150, 150)   # Another allowed region
]

FAVORED_REGIONS_GLY = [
    (-180, -60, -180, 180)  # ϕ: -180° to -60° and ψ: -150° to 150°
    ]

ALLOWED_REGIONS_GLY = [
    (-180, 180, -180, 180)  # A wider range for allowed regions
]

FAVORED_REGIONS_PRO = [
    (-75, -60, -150, 150)  # ϕ: -180° to -60° and ψ: -150° to 150°
    ]

ALLOWED_REGIONS_PRO = [
    (-90, -50, -150, 150)  # A wider range for allowed regions
]

FAVORED_REGIONS_PPRO = [
    (-60, -50, -150, 150)  # ϕ: -180° to -60° and ψ: -150° to 150°
    ]

ALLOWED_REGIONS_PPRO = [
    (-90, -40, -140, 140)  # A wider range for allowed regions
]

npz_file_path = "gaussian_density.npz"
npz_file = np.load(npz_file_path)
reference_map = npz_file["general"]
reference_map_gly = npz_file["gly"]
reference_map_pro = npz_file["pro"]
reference_map_ppro = npz_file["prepro"]
binderlen = int(sys.argv[3])

def is_in_allowed_region(phi, psi, regions):
    for region in regions:
        if region[0] <= phi <= region[1] and region[2] <= psi <= region[3]:
            return True
    return False    

def make_ramachandran_plots(phi_angles, psi_angles, outname, cmap, reference_map):
    fig, ax = plt.subplots(figsize=(19, 10.5))
    percentile_1 = np.percentile(reference_map, 60)
    percentile_2 = np.percentile(reference_map, 90)
    ax.imshow(np.rot90(reference_map), interpolation="bilinear", cmap=cmap, norm=mplcolors.BoundaryNorm(boundaries=[0, percentile_1, percentile_2, 1], ncolors=cmap.N), origin="upper", extent=(-180, 180, -180, 180))
    ax.scatter(phi_angles, psi_angles, color='blue', s=20, edgecolors="black")
    ax.set_title(f"Ramachandran Plot for {outname}", fontsize=24)
    ax.set_xlabel("Phi (ϕ) Angle (degrees)", fontsize=21)
    ax.set_ylabel("Psi (ψ) Angle (degrees)", fontsize=21)
    ax.set_xlim(-180, 180); ax.set_ylim(-180, 180)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(45))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(45))
    ax.plot([-180, 180], [0, 0], "--", linewidth=0.5, color="black")
    ax.plot([0, 0], [-180, 180], "--", linewidth=0.5, color="black")
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tight_layout()
    plt.savefig(f"{outname}_Ramachandran_{sys.argv[2]}.png", dpi=300)

# Function to compute Ramachandran plot data
def ramachandran_angles(pdb):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb)
    phi_angles, psi_angles, phi_angles_gly, psi_angles_gly, phi_angles_pro, psi_angles_pro, phi_angles_ppro, psi_angles_ppro, score = [], [], [], [], [], [], [], [], 0
    for model in structure:
        for chain in model:
            polypeptides = PDB.PPBuilder().build_peptides(chain)
            res_names = [res.get_resname() for res in chain.get_residues()]
            phi_psi_list = polypeptides[0].get_phi_psi_list()
            for m, (phi, psi) in enumerate(phi_psi_list[1:-1]):
                if res_names[m] in ['CYS', 'ASP', 'SER', 'GLN', 'LYS', 'ILE', 'THR', 'PHE', 'ASN', 'HIS', 'LEU', 'ARG', 'TRP', 'ALA', 'VAL', 'GLU', 'TYR', 'MET'] and res_names[m+1]!="PRO":
                        phi_angles.append(math.degrees(phi)); psi_angles.append(math.degrees(psi))
                        if not is_in_allowed_region(phi, psi, FAVORED_REGIONS):
                                score += 0.5
                                if not is_in_allowed_region(phi, psi, ALLOWED_REGIONS):
                                        score += 0.5
                elif res_names[m] in ['CYS', 'ASP', 'SER', 'GLN', 'LYS', 'ILE', 'THR', 'PHE', 'ASN', 'HIS', 'LEU', 'ARG', 'TRP', 'ALA', 'VAL', 'GLU', 'TYR', 'MET'] and res_names[m+1]=="PRO":
                        phi_angles_ppro.append(math.degrees(phi)); psi_angles_ppro.append(math.degrees(psi))
                        if not is_in_allowed_region(phi, psi, FAVORED_REGIONS_PPRO):
                                score += 0.5
                                if not is_in_allowed_region(phi, psi, ALLOWED_REGIONS_PPRO):
                                        score += 0.5
                elif res_names[m]=="GLY":
                        phi_angles_gly.append(math.degrees(phi)); psi_angles_gly.append(math.degrees(psi))
                        if not is_in_allowed_region(phi, psi, FAVORED_REGIONS_GLY):
                                score += 0.5
                                if not is_in_allowed_region(phi, psi, FAVORED_REGIONS_GLY):
                                        score += 0.5
                elif res_names[m]=="PRO":
                        phi_angles_pro.append(math.degrees(phi)); psi_angles_pro.append(math.degrees(psi))
                        if not is_in_allowed_region(phi, psi, FAVORED_REGIONS_PRO):
                                score += 0.5
                                if not is_in_allowed_region(phi, psi, ALLOWED_REGIONS_PRO):
                                        score += 0.5

    make_ramachandran_plots(phi_angles, psi_angles, pdb.split('.')[0], mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']), reference_map)
    make_ramachandran_plots(phi_angles_ppro, psi_angles_ppro, pdb.split('.')[0] + " Pre-Pro", mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']), reference_map_ppro)
    make_ramachandran_plots(phi_angles_gly, psi_angles_gly, pdb.split('.')[0] + " Gly", mplcolors.ListedColormap(['#FFFFFF', '#FFE8C5', '#FFCC7F']), reference_map_gly)
    make_ramachandran_plots(phi_angles_pro, psi_angles_pro, pdb.split('.')[0] + " Pro", mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']), reference_map_pro)

    return score

# Function to calculate steric clashes (atoms that are too close)
def calculate_clashes(pdb, cutoff=0.3):
    struct = md.load_pdb(pdb)
    ca_atoms = struct.topology.select('name CA')
    ca_atom_pairs = np.array(list(itertools.product(ca_atoms, ca_atoms)))
    distances_aux = md.compute_distances(trj, atom_pairs=ca_atom_pairs)
    distances = distances_aux[distances_aux != 0]
    clash_count = distances[distances < cutoff].shape[0]
                
    return clash_count

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

def generate_scoredict(pickle_filename) -> None:
        '''
        Collect the confidence values, slicing them to the binder and target regions
        then add the parsed scores to the score_dict
        '''
        pickle_file = open(pickle_filename,"rb")
        confidences = pickle.load(pickle_file)

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

        return score_dict["plddt_total"], score_dict["plddt_binder"], score_dict["plddt_target"], score_dict["pae_binder"], score_dict["pae_target"], score_dict["pae_interaction_binder"], score_dict["pae_interaction_target"], score_dict["pae_interaction"], score_dict["iptm"]

ene_files = [i for i in os.listdir(sys.argv[1]) if os.path.isdir(os.path.join(sys.argv[1],i))]

ramachandran_data, clash_data, pLDDT_total, pLDDT_antigen, pLDDT_antibody, pae_antigen, pae_antibody, pae_interaction_antigen, pae_interaction_antibody, pae_interaction, iptm = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}

for ene_file in ene_files:
     ramachandran_data[ene_file], clash_data[ene_file], pLDDT_total[ene_file], pLDDT_antigen[ene_file], pLDDT_antibody[ene_file], pae_antigen[ene_file], pae_antibody[ene_file], pae_interaction_antigen[ene_file], pae_interaction_antibody[ene_file], pae_interaction[ene_file], iptm[ene_file] = [], [], [], [], [], [], [], [], [], [], []
     PDBs = [i for i in glob.glob(f"{sys.argv[1]}/{ene_file}/ranked*.pdb")]
     pool = mp.Pool(int(sys.argv[4]))
     pLDDT_total[ene_file], pLDDT_antigen[ene_file], pLDDT_antibody[ene_file], pae_antigen[ene_file], pae_antibody[ene_file], pae_interaction_antigen[ene_file], pae_interaction_antibody[ene_file], pae_interaction[ene_file], iptm[ene_file] = zip(*pool.map(generate_scoredict,glob.glob(f"{sys.argv[1]}/{ene_file}/result*pkl")))
     ramachandran_data[ene_file] = pool.map(ramachandran_angles,PDBs)
     clash_data[ene_file] = pool.map(calculate_clashes,PDBs)

plot_results(pLDDT_total,"pLDDT",f"pLDDT_total_{sys.argv[2]}")
plot_results(pLDDT_antigen,"pLDDT",f"pLDDT_antigen_{sys.argv[2]}")
plot_results(pLDDT_antibody,"pLDDT",f"pLDDT_antibody_{sys.argv[2]}")
plot_results(pae_antigen,"pAE",f"pae_antigen_{sys.argv[2]}")
plot_results(pae_antibody,"pAE",f"pae_antibody_{sys.argv[2]}")
plot_results(pae_interaction_antigen,"pAE_interface",f"pae_interface_antigen_{sys.argv[2]}")
plot_results(pae_interaction_antibody,"pAE_interface",f"pae_interface_antibody_{sys.argv[2]}")
plot_results(pae_interaction,"pAE_interface",f"pae_interface_{sys.argv[2]}")
plot_results(iptm,"ipTM",f"ipTM_{sys.argv[2]}")

ramachandran_df = pd.DataFrame(ramachandran_data.items(), columns=['Variant', 'Ramachandran score'])
clash_df = pd.DataFrame(clash_data.items(), columns=['Variant', 'Clash score'])
df = pd.DataFrame(pLDDT_total.items(), columns=['Variant', 'pLDDT value'])
df_antigen = pd.DataFrame(pLDDT_antigen.items(), columns=['Variant', 'pLDDT value'])
df_antibody = pd.DataFrame(pLDDT_antibody.items(), columns=['Variant', 'pLDDT value'])
df_pae_antigen = pd.DataFrame(pae_antigen.items(), columns=['Variant', 'pAE value'])
df_pae_antibody = pd.DataFrame(pae_antibody.items(), columns=['Variant', 'pAE value'])
df_pae_interaction_antigen = pd.DataFrame(pae_interaction_antigen.items(), columns=['Variant', 'pAE_interface value'])
df_pae_interaction_antibody = pd.DataFrame(pae_interaction_antibody.items(), columns=['Variant', 'pAE_interface value'])
df_pae_interaction = pd.DataFrame(pae_interaction.items(), columns=['Variant', 'pAE_interface value'])
df_iptm = pd.DataFrame(iptm.items(), columns=['Variant', 'ipTM value'])

ramachandran_df.to_csv(f"{sys.argv[2]}_ramachandran.csv")
clash_df.to_csv(f"{sys.argv[2]}_clash.csv")
df.to_csv(f"{sys.argv[2]}_total.csv")
df_antigen.to_csv(f"{sys.argv[2]}_antigen.csv")
df_antibody.to_csv(f"{sys.argv[2]}_antibody.csv")
df_pae_antigen.to_csv(f"{sys.argv[2]}_pae_antigen.csv")
df_pae_antibody.to_csv(f"{sys.argv[2]}_pae_antibody.csv")
df_pae_interaction_antigen.to_csv(f"{sys.argv[2]}_pae_interaction_antigen.csv")
df_pae_interaction_antibody.to_csv(f"{sys.argv[2]}_pae_interaction_antibody.csv")
df_pae_interaction.to_csv(f"{sys.argv[2]}_pae_interaction.csv")
df_iptm.to_csv(f"{sys.argv[2]}_iptm.csv")
