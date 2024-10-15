import mdtraj as md
import itertools
import numpy as np
import matplotlib.pyplot as plt
import igraph as ig
import argparse as ap

# Script information
__author__ = "Sergi Rodà"
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Sergi Rodà"
__email__ = "sergi.rodallordes@bsc.es"

def parseArgs():
    """
        Parse arguments from command-line

        RETURNS
        -------
        traj : string
                  Filename of the trajectory file(s)
        top : string
                  Filename of the topology file
        distance: float
                  Distance threshold
        SPM_significance: float
                  SPM significance threshold
    """
   
    parser = ap.ArgumentParser(description='Script that returns a the most important edges between Ca atoms and the commands to display the Shortest Path Map (SPM) on PyMOL')
    required = parser.add_argument_group('required arguments')
    optional = parser._action_groups.pop()

    required.add_argument("traj", type=str, help="Trajectory file/s", nargs = '*')
    required.add_argument("top", type=str, help="Topology file")
    optional.add_argument("-D","--distance", help="Distance threshold to filter the connection of nodes (in nm)", type=float, default=0.6)
    optional.add_argument("-S","--SPM_significance", help="Shortest Path Map (SPM) significance threshold to filter the most important connection of nodes", type=float, default=0.3)
    parser._action_groups.append(optional)

    args = parser.parse_args()

    return args.traj, args.top, args.distance, args.SPM_significance

def PyMOL_SMP(central_edges):
    """
        Print the commands to have the Shortest Path Map representation in PyMOL

        OUTPUT
        ------
        It prints the commands for the PyMOL console
    """

    # Display the filtered central edges and their normalized weights
    print("Central edges with high frequency scores:")
    node_weights = {}
    for n, (u, v, weight) in enumerate(central_edges):
        #print(f"Edge ({u}, {v}): normalized weight = {weight:.4f}")
        # Get the coordinates of the Cα atoms (you need to know hw to fetch the Cα atoms for each residue)
        # Replace this part with how you fetch Cα coordinates for your residu
        if u not in node_weights:
            node_weights[u] = []
        if v not in node_weights:
            node_weights[v] = []
        node_weights[u].append(weight); node_weights[v].append(weight)
        print(f"sele ca_u, polymer and name CA and resi {u + 1}; sele ca_v, polymer and name CA and resi {v + 1}")

        # Draw a line (or cylinder) between the two residues with weight mapped to thickness/color
        print(f"distance dist_nodes_{n}, ca_u, ca_v; set dash_radius, {weight}, dist_nodes_{n}")

    for node in node_weights:
        print(f"sele node_ca, polymer and name CA and resi {node + 1}; set sphere_scale, {np.mean(node_weights[node])}, node_ca; show spheres, node_ca")

    print("set dash_gap, 0; set dash_color, black")

def load_traj(traj_file, top_file):
    trajectories = []
    if len(traj_file) != 0:
        for traj in traj_file:
            trajectories.append(md.load(traj, top = top_file))
    else:
        trajectories.append(md.load(traj_file, top = top_file))
    return trajectories

def main():
    # Load the trajectory and topology
    trj_str,top_str, distance_threshold, SPM_significance_threshold = parseArgs()
    trjs = load_traj(trj_str, top_str)#trj = md.load(trj_str, top = top_str)
    trj = md.join(trjs)
    ca_atoms = trjs[0].topology.select('name CA') #ca_atoms = trj.topology.select('name CA')
    num_atoms = len(ca_atoms)

    # Align trajectory to the first frame using C-alpha atoms (for example)
    trj.superpose(trj, 0, atom_indices=ca_atoms)

    # Reshape to a 2D array (frames, atom_coordinates) using C-alpha atoms
    positions = trj.xyz[:, ca_atoms, :].reshape(trj.n_frames, -1) # traj.xyz.reshape(traj.n_frames, -1)

    # Calculate the mean over all frames
    mean_positions = np.mean(positions, axis=0)

    # Subtract the mean positions to get fluctuations
    fluctuations = positions - mean_positions

    # Calculate the covariance matrix
    covariance_matrix = np.cov(fluctuations.T)

    # Get the standard deviations (square root of the diagonal elements)
    stddev = np.sqrt(np.diag(covariance_matrix))

    # Create the correlation matrix by normalizing the covariance matrix
    correlation_matrix = covariance_matrix / np.outer(stddev, stddev)

    # Plot the correlation matrix
    plt.imshow(correlation_matrix, cmap='coolwarm', vmin=-1, vmax=1)
    plt.colorbar(label='Correlation')
    plt.title('Positional Correlation Matrix')
    plt.savefig(f"{trj_str[0].split('.')[0]}_CorrMat.png", dpi=300)

    ## Alternative option to get the correlation matrix
    # Perform PCA using MDTraj
    #pca = md.compute_pca(traj, n_components=5)
    # Access the PCA results
    #print(pca[0])  # First principal component

    # Calculate the distance matrix
    ca_atom_pairs = np.array(list(itertools.product(ca_atoms, ca_atoms)))
    distances = md.compute_distances(trj, atom_pairs=ca_atom_pairs)  # = md.compute_distances(trj, atom_pairs=np.arange(trj.top.n_atoms))
    #contact_map = md.compute_contacts(traj, contacts='all', threshold=0.5)[0]
    # Compute the mean distance across all frames
    mean_distances = np.mean(distances, axis=0)
    n_shape = mean_distances.size; m_shape = int(np.sqrt(n_shape))
    distance_matrix = mean_distances.reshape(m_shape, m_shape)
    #distance_matrix = np.zeros((num_atoms, num_atoms))

    # Plot the mean distance matrix
    plt.figure()
    plt.imshow(distance_matrix, cmap='viridis', interpolation='none')
    plt.colorbar(label='Mean Distance (nm)')
    plt.title('Inter-Residue Mean Distance Matrix')
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    plt.savefig(f"{trj_str[0].split('.')[0]}_DistMat.png", dpi=300)

    ### Make the SPM ###

    # Initialize graph
    g = ig.Graph()
    g.add_vertices(num_atoms)

    # Add edges based on distance and correlation
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):  # Avoid self-loops and duplicate edges
            if distance_matrix[i][j] < distance_threshold:
                correlation_value = correlation_matrix[i][j]
                if correlation_value != 0:  # Avoid log(0)
                    edge_distance = -np.log(abs(correlation_value))
                    g.add_edge(i, j, weight=edge_distance)

    # Use Dijkstra's algorithm to find the shortest path lengths from a given node
    #source_node = 0  # Example: start from node 0
    shortest_paths = g.shortest_paths(weights="weight") #g.shortest_paths(source=source_node, weights="weight")
    #print(shortest_paths)

    # Display the shortest paths from the source node
    #for target_node, path_length in enumerate(shortest_paths[0]):
    #    print(f"Shortest path from node {source_node} to node {target_node}: {path_length:.4f}")

    # Matrix to store frequency of edge usage
    frequency_matrix = np.zeros((num_atoms, num_atoms))

    for source_node in range(num_atoms):
        for target_node in range(num_atoms):
            if source_node != target_node:
                # Get the actual path (sequence of nodes) for the shortest path
                path = g.get_shortest_paths(source_node, to=target_node, weights="weight", output="vpath")[0]

                # Increment the frequency for each edge in the path
                for i in range(len(path) - 1):
                    u = path[i]
                    v = path[i + 1]
                    frequency_matrix[u][v] += 1
                    frequency_matrix[v][u] += 1 # Since it's an undirected graph

    # Find the maximum frequency score (fmax)
    fmax = np.max(frequency_matrix)

    # Normalize the frequency scores: fij/fmax
    for edge in g.es:
        u, v = edge.tuple
        fij = frequency_matrix[u][v]
        normalized_weight = fij / fmax if fmax > 0 else 0
        edge["normalized_weight"] = normalized_weight

    # Optionally, filter out edges with low centrality (e.g., keep only edges with normalized weight > 0.1)
    central_edges = [(e.source, e.target, e["normalized_weight"]) for e in g.es if e["normalized_weight"] > SPM_significance_threshold]

    fig, ax = plt.subplots()
    ig.plot(
        g,
        target=ax,
        layout='circle',
        vertex_color='steelblue',
        vertex_label=range(g.vcount()),#    edge_width=g.es['weight'],
        edge_label=g.es["normalized_weight"],
        edge_color='#666',
        edge_align_label=True,
        edge_background='white'
    )
    plt.savefig(f"{trj_str[0].split('.')[0]}_SPM.png", dpi=300)

    PyMOL_SMP(central_edges)

if __name__ == "__main__":
    """Call the main function"""
    main()