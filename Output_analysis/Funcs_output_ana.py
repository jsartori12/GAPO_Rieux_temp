#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 11:44:12 2024

@author: joao
"""
#### Python default libs
import time
from tqdm import tqdm
import descriptor_converter
#### Biopython 
from Bio import Align
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices
from Bio.Phylo.TreeConstruction import DistanceCalculator

#### Data manipulation
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.manifold import TSNE
import seaborn as sns
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler

#from sklearn import StandardScaler

def Read_data(path, condition):
    # Open the file in read mode
    with open(path, 'r') as file:
        # Read the content of the file
        content = file.read()
    
        # Split the content by commas and newlines
        rows = [row.split(',') for row in content.split('\n')]
    
        # Create a Pandas DataFrame
        df = pd.DataFrame(rows, columns=['Sequence', 'dG', 'Cycle'])
        df['dG'] = df['dG'].astype(float)
        df['Cycle'] = df['Cycle'].astype(str)
    
    # Print the DataFrame
    #print(df)
    df.drop(df.index[-1], inplace=True)
    df.reset_index(drop=True, inplace=True)
    df["Condition"] = condition
    df['Sequence_name'] = 'Sequence_' + df.index.astype(str) + '_' + str(condition)
    df.set_index('Sequence_name', inplace = True)
    return df

def List_min_dgs(dfs, labels_list, apt_function, out,plot_width=10, plot_height=6, wt_value=None, min_or_max = min):
    lowests_cond = []
    min_dG_rows = []
    
    if min_or_max == min:
        for df in dfs:
            # Find the lowest dG values for each Cycle
            lowests_cond.append([min(df[df["Cycle"] == str(i)][apt_function]) for i in range(len(set(df["Cycle"])))])
            
            # Get the rows with the minimum dG values for each Cycle
            min_dG_indices = df.groupby('Cycle')[apt_function].idxmin()
            min_dG_rows.append(df.loc[min_dG_indices])
    
        # Set the size of the plot
        plt.figure(figsize=(plot_width, plot_height))
        
        # Plot the values for each condition
        for i, lows in enumerate(lowests_cond):
            plt.plot(lows, linestyle='-', label=labels_list[i])
        
        # Optionally add a red line for the WT value
        if wt_value is not None:
            # Get the x-axis limits
            x_min, x_max = plt.xlim()
            plt.axhline(y=wt_value, color='red', linestyle='--', xmin=x_min, xmax=x_max, label='WT Value')
    if min_or_max == max:
        for df in dfs:
            # Find the lowest dG values for each Cycle
            lowests_cond.append([max(df[df["Cycle"] == str(i)][apt_function]) for i in range(len(set(df["Cycle"])))])
            
            # Get the rows with the minimum dG values for each Cycle
            min_dG_indices = df.groupby('Cycle')[apt_function].idxmin()
            min_dG_rows.append(df.loc[min_dG_indices])

        # Set the size of the plot
        plt.figure(figsize=(plot_width, plot_height))
        
        # Plot the values for each condition
        for i, lows in enumerate(lowests_cond):
            plt.plot(lows, linestyle='-', label=labels_list[i])
        
        # Optionally add a red line for the WT value
        if wt_value is not None:
            # Get the x-axis limits
            x_min, x_max = plt.xlim()
            plt.axhline(y=wt_value, color='red', linestyle='--', xmin=x_min, xmax=x_max, label='WT Value')
    # Add labels and title
    plt.xlabel('Cycle')
    plt.ylabel(f'Lowest {apt_function}')
    plt.title('Genetic Algorithm Conditions')
    # Add a legend
    plt.legend()
    
    # Save and show the plot
    plt.savefig(f'best_per_cycle_{apt_function}_{out}.png', dpi=300)
    plt.show()

    return min_dG_rows

# def List_min_dgs(df1, df2, df3, labels_list, plot_width=10, plot_height=6, wt_value=None):
#     lowests_cond1 = [min(df1[df1["Cycle"] == str(i)]["dG"]) for i in range(len(set(df1["Cycle"])))]
#     lowests_cond2 = [min(df2[df2["Cycle"] == str(i)]["dG"]) for i in range(len(set(df2["Cycle"])))]
#     lowests_cond3 = [min(df3[df3["Cycle"] == str(i)]["dG"]) for i in range(len(set(df3["Cycle"])))]
    
#     min_dG_indices_1 = df1.groupby('Cycle')['dG'].idxmin()
#     min_dG_rows_1 = df1.loc[min_dG_indices_1]
    
#     min_dG_indices_2 = df2.groupby('Cycle')['dG'].idxmin()
#     min_dG_rows_2 = df2.loc[min_dG_indices_2]
    
#     min_dG_indices_3 = df3.groupby('Cycle')['dG'].idxmin()
#     min_dG_rows_3 = df3.loc[min_dG_indices_3]
    
#     # Set the size of the plot
#     plt.figure(figsize=(plot_width, plot_height))
    
#     # Plot the values
#     plot1, = plt.plot(lowests_cond1, linestyle='-', label=labels_list[0])
#     plot2, = plt.plot(lowests_cond2, linestyle='-', label=labels_list[1])
#     plot3, = plt.plot(lowests_cond3, linestyle='-', label=labels_list[2])
    
#     # Optionally add a red line for the WT value
#     if wt_value is not None:
#         # Get the x-axis limits
#         x_min, x_max = plt.xlim()
#         plt.axhline(y=wt_value, color='red', linestyle='--', xmin=x_min, xmax=x_max, label='WT Value')


#     # Add labels and title
#     plt.xlabel('Cycle')
#     plt.ylabel('Lowest dG')
#     plt.title('Genetic Algorithm Conditions')
#     # Add a legend
#     plt.legend()
    
#     # Save and show the plot
#     plt.savefig('best_per_cycle.png', dpi=300)
#     plt.show()

#     return min_dG_rows_1, min_dG_rows_2, min_dG_rows_3

# def Plot_bests(df1, df2, df3, natural, out):
    
#     lowests_cond1 = [min(df1[df1["Cycle"] == str(i)]["dG"]) for i in range(len(set(df1["Cycle"])))]
#     lowests_cond2 = [min(df2[df2["Cycle"] == str(i)]["dG"]) for i in range(len(set(df2["Cycle"])))]
#     lowests_cond3 = [min(df3[df3["Cycle"] == str(i)]["dG"]) for i in range(len(set(df3["Cycle"])))]
#     min_dG_indices_1 = df1.groupby('Cycle')['dG'].idxmin()
#     min_dG_rows_1 = df1.loc[min_dG_indices_1]
#     plot1, = plt.plot(lowests_cond1, linestyle='-', label = "Replica 1")
#     plot2, = plt.plot(lowests_cond2, linestyle='-', label = "Replica 2")
#     plot3, = plt.plot(lowests_cond3, linestyle='-', label = "Replica 3")
#     plot4, = plt.plot(natural, linestyle=':', color = "Red", label = "Original")
#     plt.gcf().set_size_inches(10, 6)
#     # Add labels and title
#     #plt.xlabel('Cycle')
#    # plt.ylabel('Lowest dG')
#     #plt.title('Genetic Algorithm Conditions')
#     #plt.legend()
#     plt.savefig(f'{out}.png', dpi=300, transparent=True)
        
#     fig, ax = plt.subplots()
#     plt.ylim(-20, 0)
#     ax.set_xlim(0, len(lowests_cond1) - 1)
#     line1, = ax.plot([], [], linestyle='-', label="Replica 1")
#     line2, = ax.plot([], [], linestyle='-', label="Replica 2")
#     line3, = ax.plot([], [], linestyle='-', label="Replica 2")
#     line4, = ax.plot([], [], linestyle=':', color="Red", label="Original")
#     ax.legend()
    
#     def update(frame):
#         if frame == len(lowests_cond1):
#             anim.event_source.stop()
#         line1.set_data(range(frame), lowests_cond1[:frame])
#         line2.set_data(range(frame), lowests_cond2[:frame])
#         line3.set_data(range(frame), lowests_cond3[:frame])
#         line4.set_data(range(frame), natural[:frame])
#         return line1, line2, line3, 
    


#     anim = FuncAnimation(fig, update, frames=len(lowests_cond1) + 1, interval=200)
    
#     anim.save('animation_pbee.gif', writer='pillow',dpi=300)
#     plt.show()

#     #plt.xlim(0, 151)  # Adjust the limits according to your needs
#     # Add a legend
#     return min_dG_rows_1


def Plot_bests(dfs, natural, apt_function,out):
    lowests_cond = []
    
    # Extract lowest dG values for each DataFrame
    for df in dfs:
        lowests_cond.append([min(df[df["Cycle"] == str(i)][apt_function]) for i in range(len(set(df["Cycle"])))])

    # Plot and save static plot with all lines
    plt.figure(figsize=(10, 6))
    
    for i, lows in enumerate(lowests_cond):
        plt.plot(lows, linestyle='-', label=f"Replica {i + 1}")
    
    plt.plot(natural, linestyle=':', color="Red", label="Original")
    plt.ylim(-20, 0)
    plt.xlim(0, len(lowests_cond[0]) - 1)
    plt.xlabel('Cycle')
    plt.ylabel(f'Lowest {apt_function}')
    plt.title('Genetic Algorithm Conditions')
    plt.legend()
    plt.savefig(f'{out}.png', dpi=300, transparent=True)
    plt.close()

    # Set up the figure and axis for the animation
    fig, ax = plt.subplots()
    ax.set_xlim(0, len(lowests_cond[0]) - 1)
    ax.set_ylim(-20, 0)
    
    # Static red line for natural condition
    ax.plot(natural, linestyle=':', color="Red", label="Original")

    # Initialize empty line objects for each condition
    lines = [ax.plot([], [], linestyle='-', label=f"Replica {i + 1}")[0] for i in range(len(lowests_cond))]
    
    ax.legend()
    
    def update(frame):
        for line, lows in zip(lines, lowests_cond):
            line.set_data(range(frame), lows[:frame])
        return lines

    # Create and save the animation
    anim = FuncAnimation(fig, update, frames=len(lowests_cond[0]) + 1, interval=200, blit=True)
    anim.save(f'{out}_animation.gif', writer='pillow', dpi=300)
    plt.show()

# def Plot_bests(df1, df2, df3, natural, out):
    
#     lowests_cond1 = [min(df1[df1["Cycle"] == str(i)]["dG"]) for i in range(len(set(df1["Cycle"])))]
#     lowests_cond2 = [min(df2[df2["Cycle"] == str(i)]["dG"]) for i in range(len(set(df2["Cycle"])))]
#     lowests_cond3 = [min(df3[df3["Cycle"] == str(i)]["dG"]) for i in range(len(set(df3["Cycle"])))]
#     min_dG_indices_1 = df1.groupby('Cycle')['dG'].idxmin()
#     min_dG_rows_1 = df1.loc[min_dG_indices_1]

#     # Plot and save static plot with all lines
#     plt.figure(figsize=(10, 6))
#     plt.plot(lowests_cond1, linestyle='-', label="Replica 1")
#     plt.plot(lowests_cond2, linestyle='-', label="Replica 2")
#     plt.plot(lowests_cond3, linestyle='-', label="Replica 3")
#     plt.plot(natural, linestyle=':', color="Red", label="Original")
#     plt.ylim(-20, 0)
#     plt.xlim(0, len(lowests_cond1) - 1)
#     plt.xlabel('Cycle')
#     plt.ylabel('Lowest dG')
#     plt.title('Genetic Algorithm Conditions')
#     plt.legend()
#     plt.savefig(f'{out}.png', dpi=300, transparent=True)
#     plt.close()

#     # Set up the figure and axis for the animation
#     fig, ax = plt.subplots()
#     ax.set_xlim(0, len(lowests_cond1) - 1)
#     ax.set_ylim(-20, 0)
    
#     # Static red line
#     ax.plot(natural, linestyle=':', color="Red", label="Original")

#     line1, = ax.plot([], [], linestyle='-', label="Replica 1")
#     line2, = ax.plot([], [], linestyle='-', label="Replica 2")
#     line3, = ax.plot([], [], linestyle='-', label="Replica 3")
#     ax.legend()
    
#     def update(frame):
#         line1.set_data(range(frame), lowests_cond1[:frame])
#         line2.set_data(range(frame), lowests_cond2[:frame])
#         line3.set_data(range(frame), lowests_cond3[:frame])
#         return line1, line2, line3

#     anim = FuncAnimation(fig, update, frames=len(lowests_cond1) + 1, interval=200, blit=True)
#     anim.save('animation_pbee.gif', writer='pillow', dpi=300)
#     plt.show()

def Get_best_ind(df):
    
    # Find the index of the row with the minimum value in the 'Score' column
    min_index = df['dG'].idxmin()

    # Use the index to select the corresponding row
    min_row = df.loc[min_index]
    
    return min_row

    

def Distance_matrix_init_final(df, out):
    #init_vs_final = pd.concat([df[df["Cycle"] == str(0)], df[df["Cycle"] == str(150)]])
    init_vs_final = df[df["Cycle"] == str(150)]
    init_vs_final = init_vs_final.reset_index()
    
    
    # Create an alignment file in FASTA format for synthetic sequences
    alignment_file = (f"{out}.fasta")
    with open(alignment_file, "w") as f:
        for i in init_vs_final.index:
            f.write(f">{init_vs_final['index'].iloc[i]}\n{init_vs_final['Sequence'].iloc[i]}\n")
    
    
    alignment = AlignIO.read(alignment_file, "fasta")
    
    # Create a distance calculator
    calculator = DistanceCalculator('blosum62')
    dm = calculator.get_distance(alignment)
    
    df1 = pd.DataFrame(dm.matrix)
    return df1

def Create_align_fasta(df, out):
    # Create an alignment file in FASTA format for synthetic sequences
    alignment_file = (f"{out}.fasta")
    with open(alignment_file, "w") as f:
        for i in df.index:
            f.write(f">{i}#{df['dG'][i]}\n{df['Sequence'][i]}\n")
    
    return alignment_file

def Diss_dist(alg_obj):
    
    ############## calculate dissimilatiry matrix
    alignment = AlignIO.read(alg_obj, "fasta")
    
    # Create a distance calculator
    calculator = DistanceCalculator('blosum62')
    dm = calculator.get_distance(alignment)
    
    df1 = pd.DataFrame(dm.matrix)
    return df1

def algn_pairwise(seq1, seq2):
    """
    Perform pairwise sequence alignment and return a list of alignment objects.

    Parameters:
    - seq1 (str): The first input sequence for alignment.
    - seq2 (str): The second input sequence for alignment.

    Returns:
    list: A list containing pairwise alignment objects.
    """
    #### https://biopython.org/docs/1.77/api/Bio.Align.html
    # Define your sequences
    ref_seq = seq1
    test_seq = seq2
    
    # Create a pairwise aligner with specified substitution matrix (BLOSUM62)
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

    # Perform pairwise sequence alignment
    alignment = aligner.align(ref_seq, test_seq)

    # Convert the alignment object into a list of alignments
    alignments = list(alignment[0])

    return alignments

def calculate_indentity(seq1, seq2):
    """
    Compute the total number of identical matches between two alignments.

    Parameters:
    - align1 (str): The first sequence alignment.
    - align2 (str): The second sequence alignment.

    Returns:
    int: The total number of identical matches between the two alignments.
    """
    
    # Initialize the counter for identical matches
    compute_total = 0
    
    align1, align2 = algn_pairwise(seq1, seq2)
    
    
    # Iterate over each pair of aligned residues
    for j, (res1, res2) in enumerate(zip(align1, align2)):
        # If neither residue is a gap ('-') and the residues are identical, update the count
        if res1 != "-" and res2 != "-" and res1 == res2:
            compute_total += 1
    
    sequence_length = len(seq1)
    
    # Calculate sequence identity query vs reference
    identity = (compute_total / sequence_length) * 100
    
    return round(identity, 2)

def calculate_identity2(seq1, seq2):
    """
    Calculate the identity score between two sequences of the same length using NumPy.

    Parameters:
    - seq1: First sequence (string or list/array of characters).
    - seq2: Second sequence (string or list/array of characters).

    Returns:
    - float: Identity score (proportion of positions with identical residues).
    """
    # Convert sequences to NumPy arrays if they are strings
    seq1_array = np.array(list(seq1))
    seq2_array = np.array(list(seq2))

    # Ensure the sequences are of the same length
    if seq1_array.shape[0] != seq2_array.shape[0]:
        raise ValueError("Sequences must be of the same length")

    # Calculate the number of identical positions
    identical_positions = np.sum(seq1_array == seq2_array)

    # Calculate the identity score as the fraction of identical positions
    identity_score = identical_positions / len(seq1_array)

    return identity_score

# def calculate_identity_matrix(df, sequence_column='Sequence'):
#     """
#     Calculate identity matrix for a group of sequences in a DataFrame.

#     Parameters:
#     - df: DataFrame containing the sequences.
#     - sequence_column: Name of the column containing the sequences.

#     Returns:
#     - DataFrame: Identity matrix.
#     """
#     sequences = df[sequence_column]

#     # Initialize an empty matrix
#     identity_matrix = pd.DataFrame(index=sequences.index, columns=sequences.index, dtype=float)

#     # Calculate identity scores and fill the matrix
#     with tqdm(total=len(sequences)**2) as pbar:
#         for i, seq1 in enumerate(sequences):
#             for j, seq2 in enumerate(sequences):
#                 #print(seq1, seq2)
    
#                 if i == j:
#                     # Diagonal elements (self-comparison) are set to 1
#                     identity_matrix.iloc[i, j] = 0
#                 else:
#                     identity = calculate_indentity(seq1, seq2)
#                     # Fill the matrix with the identity score
#                     identity_matrix.iloc[i, j] = identity
#                 pbar.update(1)  # Update the progress bar
#     return identity_matrix
def calculate_identity_matrix_distance(df,out, sequence_column='Sequence'):
    """
    Calculate identity matrix for a group of sequences in a DataFrame using NumPy.

    Parameters:
    - df: DataFrame containing the sequences.
    - sequence_column: Name of the column containing the sequences.

    Returns:
    - np.ndarray: Identity matrix as a NumPy array.
    """
    sequences = df[sequence_column].values  # Extract sequences as a NumPy array
    num_sequences = len(sequences)

    # Initialize an empty identity matrix with zeros
    identity_matrix = np.zeros((num_sequences, num_sequences), dtype=float)

    # Calculate identity scores and fill the matrix
    with tqdm(total=num_sequences**2) as pbar:
        for i, seq1 in enumerate(sequences):
            for j, seq2 in enumerate(sequences):
                if i == j:
                    # Diagonal elements (self-comparison) are set to 1
                    identity_matrix[i, j] = 0  # Adjust if needed to set to 1 for identity
                else:
                    # Compute identity between two sequences
                    identity = calculate_identity2(seq1, seq2)
                    # Fill the matrix with the identity score
                    identity_matrix[i, j] = 1 - identity
                pbar.update(1)  # Update the progress bar
     # Convert NumPy array to Pandas DataFrame
    identity_df = pd.DataFrame(identity_matrix)

    # Save the DataFrame as a CSV file
    identity_df.to_csv(f"{out}_distance_matrix.csv", index=False)
    
    return identity_matrix

def calculate_identity_matrix_distance2(df, out, sequence_column='Sequence'):
    """
    Calculate identity matrix for a group of sequences in a DataFrame using NumPy, optimized by calculating only the lower triangle.
    
    Parameters:
    - df: DataFrame containing the sequences.
    - sequence_column: Name of the column containing the sequences.

    Returns:
    - np.ndarray: Identity matrix as a NumPy array.
    """
    sequences = df[sequence_column].values  # Extract sequences as a NumPy array
    num_sequences = len(sequences)

    # Initialize an empty identity matrix with zeros
    identity_matrix = np.zeros((num_sequences, num_sequences), dtype=float)

    # Calculate identity scores for only one triangle
    with tqdm(total=num_sequences * (num_sequences + 1) // 2) as pbar:
        for i in range(num_sequences):
            for j in range(i, num_sequences):  # Only compute for j >= i (lower triangle)
                if i == j:
                    # Diagonal elements (self-comparison) are set to 0
                    identity_matrix[i, j] = 0  # Set to 1 if you want identity for self-comparison
                else:
                    # Compute identity between two sequences
                    identity = calculate_identity2(sequences[i], sequences[j])
                    # Fill the lower triangle and mirror to the upper triangle
                    identity_matrix[i, j] = 1 - identity
                    identity_matrix[j, i] = 1 - identity  # Symmetric copy
                pbar.update(1)  # Update the progress bar

    # Convert NumPy array to Pandas DataFrame
    identity_df = pd.DataFrame(identity_matrix)

    # Save the DataFrame as a CSV file
    identity_df.to_csv(f"{out}_distance_matrix.csv", index=False)

    return identity_matrix

def red_tsne(data,dimension,prp, label_list, out):
    X = np.array(data)
    tsn = TSNE(n_components= dimension,perplexity= prp,random_state= 42,n_iter= 5000,n_jobs= -1)
    X_embedded = tsn.fit_transform(X)
    X_embedded.shape
    dt_tsne_all = pd.DataFrame(X_embedded)
    dt_tsne_all["Conditions"] = label_list
    dt_tsne_all.to_csv(f"{out}_tsne.csv")
    return dt_tsne_all

def pca_dataset(dataset, dimensions, label_list):

    pca_features= PCA(n_components=dimensions)
    principalComponents_features = pca_features.fit_transform(dataset)

    pca_features_all = PCA(n_components= min(np.shape(dataset)))
    principalComponents_features_all = pca_features_all.fit_transform(dataset)

    x=0
    explained_variance = [(x+pca_features.explained_variance_[i]/np.sum(pca_features_all.explained_variance_))*100 for i in range(dimensions)]
    
    df_pca = pd.DataFrame(principalComponents_features)
    
    df_pca["Conditions"] = label_list
    
    return df_pca, explained_variance


# def Ploting_tSNE(df):
    
#     fig = plt.figure(figsize=(8, 8))
#     sns.scatterplot(df, x=df.iloc[:,0], y=df.iloc[:,1], hue = df.iloc[:,2],palette='rocket', alpha=1.0)    
#     #plt.legend = sns.legend(title=df.columns[2], bbox_to_anchor=(1, 1), loc='upper left')
#     fig.suptitle('t-SNE in Identity matrix data')
#     fig.savefig("tSNE.png", dpi=600, transparent=True, bbox_inches='tight', pad_inches=0)

# def Ploting_tSNE(df):
#     fig = plt.figure(figsize=(10, 10))
#     ax = fig.add_subplot(1, 1, 1)
#     sns.scatterplot(df, x=df.iloc[:,0], y=df.iloc[:,1], hue=df.iloc[:,2], palette='rocket', alpha=1.0)
#     #sns.scatterplot(df, x=df.iloc[:,0], y=df.iloc[:,1], palette='rocket', alpha=1.0, ax=ax)
#     #ax.set_title('t-SNE in Identity matrix data')
    
#     # Set the plot boundaries to include a bit of margin
#     # ax.set_xlim(df.iloc[:,0].min() - 0.5, df.iloc[:,0].max() + 0.5)
#     # ax.set_ylim(df.iloc[:,1].min() - 0.5, df.iloc[:,1].max() + 0.5)
#     ax.get_legend().remove()
#     # Save the plot with transparent background, but keep the frame
#     fig.savefig("tSNE_new.png", dpi=600, transparent=True, bbox_inches='tight')
#     plt.close(fig)  # close the figure to release memory

# def Ploting_tSNE(df):
#     fig = plt.figure(figsize=(10, 10))
#     ax = fig.add_subplot(1, 1, 1)


# def Ploting_tSNE_panel(df, out, column_name):
#     unique_values = df[column_name].unique()
#     num_unique = len(unique_values)
    
#     # Create subplots for Population and dG, arranged in a grid
#     fig, axs = plt.subplots(2, num_unique, figsize=(5 * num_unique, 10))
    
#     for i, val in enumerate(unique_values):
#         # Filter the DataFrame based on the unique value
#         df_filtered = df[df[column_name] == val]
        
#         # Plot for Population
#         palette_pop = sns.dark_palette("#69d", reverse=True, as_cmap=False)
#         sns.scatterplot(
#             data=df_filtered,
#             x=df_filtered.iloc[:, 0], 
#             y=df_filtered.iloc[:, 1], 
#             hue=df_filtered["Population"],
#             palette=palette_pop, 
#             alpha=1.0, 
#             ax=axs[0, i]
#         )
#         axs[0, i].set_title(f"Population - {column_name}: {val}")
#         axs[0, i].get_legend().remove()
        
#         # Plot for dG
#         palette_dg = sns.color_palette("plasma_r", as_cmap=True)
#         sns.scatterplot(
#             data=df_filtered,
#             x=df_filtered.iloc[:, 0], 
#             y=df_filtered.iloc[:, 1], 
#             hue=df_filtered["dG"],
#             palette=palette_dg, 
#             alpha=1.0, 
#             ax=axs[1, i]
#         )
#         axs[1, i].set_title(f"dG - {column_name}: {val}")
#         axs[1, i].get_legend().remove()

#     # Adjust the layout to avoid overlapping
#     plt.tight_layout()
    
#     # Save the plot
#     fig.savefig(f"{out}_panel.png", dpi=600, bbox_inches='tight')
#     fig.savefig(f"{out}_panel_transp.png", dpi=600, transparent=True, bbox_inches='tight')
    
#     # Close the figure to release memory
#     plt.close(fig)

# def Ploting_tSNE(df, out, color_method):
#     fig = plt.figure(figsize=(10, 10))
#     ax = fig.add_subplot(1, 1, 1)
    
#     # Define a mapping from the index values (1, 2, 3) to specific markers
#     #marker_mapping = {1: 'o', 2: 's', 3: 'D'}  # 'o' for circle, 's' for square, 'D' for diamond
    
#     # Add a new column to the DataFrame with mapped markers
#     #df['marker'] = df.iloc[:, 3].map(marker_mapping)
#     if color_method == "population":
#         palette = sns.dark_palette("#69d", reverse=True, as_cmap=False)
#         #palette = sns.color_palette("dark:salmon_r")
#         sns.scatterplot(
#             data=df,
#             x=df.iloc[:, 0], 
#             y=df.iloc[:, 1], 
#             hue=df["Population"],
#             palette=palette, 
#             alpha=1.0, 
#             ax=ax
#         )
        
#         ax.get_legend().remove()
        
#         # Save the plot with transparent background, but keep the frame
#         fig.savefig(f"{out}_new.png", dpi=600, bbox_inches='tight')
#         fig.savefig(f"{out}_new_transp.png", dpi=600, transparent=True, bbox_inches='tight')
#         plt.close(fig)  # close the figure to release memory
#     if color_method == "dG":
#         palette = sns.color_palette("plasma_r", as_cmap=True)
        

#         # Plot using the style argument based on the index column (the 4th column in the DataFrame)
#         sns.scatterplot(
#             data=df,
#             x=df.iloc[:, 0], 
#             y=df.iloc[:, 1], 
#             hue=df["dG"],
#             palette=palette, 
#             alpha=1.0, 
#             ax=ax
#         )
        
#         ax.get_legend().remove()
        
#         # Save the plot with transparent background, but keep the frame
#         fig.savefig(f"{out}_new.png", dpi=600, bbox_inches='tight')
#         fig.savefig(f"{out}_new_transp.png", dpi=600, transparent=True, bbox_inches='tight')
#         plt.close(fig)  # close the figure to release memory

def plot_panel(df):
    # Set the style and figure
    sns.set(style="whitegrid")
    fig, axs = plt.subplots(1, 2, figsize=(14, 6))
    
    # First plot: Colored by Population
    scatter1 = axs[0].scatter(df.iloc[:, 0], df.iloc[:, 1], c=df['Population'], cmap='viridis', s=50, edgecolor='k')
    axs[0].set_title('Scatter Plot Colored by Population')
    axs[0].set_xlabel('X Dimension')
    axs[0].set_ylabel('Y Dimension')
    fig.colorbar(scatter1, ax=axs[0], label='Population')

    # Second plot: Colored by dG
    scatter2 = axs[1].scatter(df.iloc[:, 0], df.iloc[:, 1], c=df['dG'], cmap='coolwarm', s=50, edgecolor='k')
    axs[1].set_title('Scatter Plot Colored by dG')
    axs[1].set_xlabel('X Dimension')
    axs[1].set_ylabel('Y Dimension')
    fig.colorbar(scatter2, ax=axs[1], label='dG')
    
    # Adjust the layout
    plt.tight_layout()
    plt.show()

def plot_panel_3d(df, out,df_lowests,dG_min=None, dG_max=None):
    
    ##### do the same for the lowests
    # Set the style and figure
    sns.set(style="whitegrid")
    fig = plt.figure(figsize=(14, 6))
    
    # First plot: Colored by Population (3D)
    ax1 = fig.add_subplot(121, projection='3d')
    scatter1 = ax1.scatter(df_lowests.iloc[:, 0], df_lowests.iloc[:, 1], df_lowests['dG'], c=df_lowests['Population'], cmap='viridis', s=50, edgecolor='k')
    ax1.set_title('3D Scatter Plot Colored by Population')
    ax1.set_xlabel('X Dimension')
    ax1.set_ylabel('Y Dimension')
    ax1.set_zlabel('')
    if dG_min is not None and dG_max is not None:
        ax1.set_zlim(dG_min, dG_max)
    fig.colorbar(scatter1, ax=ax1, label='Population')

    # Second plot: Colored by dG (3D)
    ax2 = fig.add_subplot(122, projection='3d')
    scatter2 = ax2.scatter(df_lowests.iloc[:, 0], df_lowests.iloc[:, 1], df_lowests['dG'], c=df_lowests['dG'], cmap='coolwarm', s=50, edgecolor='k')
    ax2.set_title('3D Scatter Plot Colored by dG')
    #ax2.set_xlabel('X Dimension')
    #ax2.set_ylabel('Y Dimension')
    ax2.set_zlabel('')
    # Optionally set dG limits
    if dG_min is not None and dG_max is not None:
        ax2.set_zlim(dG_min, dG_max)
    fig.colorbar(scatter2, ax=ax2, label='')
    
    # Adjust the layout
    plt.tight_layout()
    plt.show()
    # Save the plot
    fig.savefig(f"{out}_lowests_panel.png", dpi=600, bbox_inches='tight')
    fig.savefig(f"{out}_lowests_panel_transp.png", dpi=600, transparent=True, bbox_inches='tight')

    #######
    # Set the style and figure
    sns.set(style="whitegrid")
    fig = plt.figure(figsize=(14, 6))
    
    # First plot: Colored by Population (3D)
    ax1 = fig.add_subplot(121, projection='3d')
    scatter1 = ax1.scatter(df.iloc[:, 0], df.iloc[:, 1], df['dG'], c=df['Population'], cmap='viridis', s=50, edgecolor='k')
    ax1.set_title('3D Scatter Plot Colored by Population')
    ax1.set_xlabel('X Dimension')
    ax1.set_ylabel('Y Dimension')
    ax1.set_zlabel('')
    if dG_min is not None and dG_max is not None:
        ax1.set_zlim(dG_min, dG_max)
    fig.colorbar(scatter1, ax=ax1, label='Population')

    # Second plot: Colored by dG (3D)
    ax2 = fig.add_subplot(122, projection='3d')
    scatter2 = ax2.scatter(df.iloc[:, 0], df.iloc[:, 1], df['dG'], c=df['dG'], cmap='coolwarm', s=50, edgecolor='k')
    ax2.set_title('3D Scatter Plot Colored by dG')
    #ax2.set_xlabel('X Dimension')
    #ax2.set_ylabel('Y Dimension')
    ax2.set_zlabel('')
    # Optionally set dG limits
    if dG_min is not None and dG_max is not None:
        ax2.set_zlim(dG_min, dG_max)
    fig.colorbar(scatter2, ax=ax2, label='')
    
    # Adjust the layout
    plt.tight_layout()
    plt.show()
    # Save the plot
    fig.savefig(f"{out}_panel.png", dpi=600, bbox_inches='tight')
    fig.savefig(f"{out}_panel_transp.png", dpi=600, transparent=True, bbox_inches='tight')





def Ploting_tSNE_panel(df, out, column_name, plot_lowest_dG=False):
    unique_values = df[column_name].unique()
    num_unique = len(unique_values)
    
    # Create subplots for Population and dG, arranged in a grid
    fig, axs = plt.subplots(2, num_unique, figsize=(5 * num_unique, 10))
    
    for i, val in enumerate(unique_values):
        # Filter the DataFrame based on the unique value
        df_filtered = df[df[column_name] == val]
        
        if plot_lowest_dG:
            # Keep only the rows with the lowest dG for each unique population
            min_dg_rows = df_filtered.loc[df_filtered.groupby('Population')['dG'].idxmin()]
            df_filtered = min_dg_rows

        # Plot for Population
        palette_pop = sns.color_palette("crest", as_cmap=True)
        sns.scatterplot(
            data=df_filtered,
            x=df_filtered.iloc[:, 0], 
            y=df_filtered.iloc[:, 1], 
            hue=df_filtered["Population"],
            palette=palette_pop, 
            alpha=1.0, 
            ax=axs[0, i]
        )
        axs[0, i].set_title(f"Population - {column_name}: {val}")
        axs[0, i].get_legend().remove()
        # Remove x and y labels
        axs[0, i].set_xlabel(None)
        axs[0, i].set_ylabel(None)
        
        # Add labels for each point if plot_lowest_dG is True
        if plot_lowest_dG:
            for j, row in df_filtered.iterrows():
                axs[0, i].annotate(
                    f"{row['Population']}", 
                    (row.iloc[0], row.iloc[1]), 
                    textcoords="offset points", 
                    xytext=(5, 5), 
                    ha='center'
                )
        
        # Plot for dG
        palette_dg = sns.color_palette("plasma_r", as_cmap=True)
        sns.scatterplot(
            data=df_filtered,
            x=df_filtered.iloc[:, 0], 
            y=df_filtered.iloc[:, 1], 
            hue=df_filtered["dG"],
            palette=palette_dg, 
            alpha=1.0, 
            ax=axs[1, i]
        )
        axs[1, i].set_title(f"dG - {column_name}: {val}")
        axs[1, i].get_legend().remove()
        # Remove x and y labels
        axs[1, i].set_xlabel(None)
        axs[1, i].set_ylabel(None)
        
        # Add labels for each point if plot_lowest_dG is True
        if plot_lowest_dG:
            for j, row in df_filtered.iterrows():
                axs[1, i].annotate(
                    f"{row['dG']:.2f}", 
                    (row.iloc[0], row.iloc[1]), 
                    textcoords="offset points", 
                    xytext=(5, 5), 
                    ha='center'
                )

    # Adjust the layout to avoid overlapping
    plt.tight_layout()
    
    # Save the plot
    fig.savefig(f"{out}_panel.png", dpi=600, bbox_inches='tight')
    fig.savefig(f"{out}_panel_transp.png", dpi=600, transparent=True, bbox_inches='tight')
    
    # Close the figure to release memory
    plt.close(fig)

def Ploting_tSNE_for_all(df, out, color_method):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1)
    
    # Define a mapping from the index values (1, 2, 3) to specific markers
    #marker_mapping = {1: 'o', 2: 's', 3: 'D'}  # 'o' for circle, 's' for square, 'D' for diamond
    
    # Add a new column to the DataFrame with mapped markers
    #df['marker'] = df.iloc[:, 3].map(marker_mapping)
    if color_method == "population":
        palette = sns.dark_palette("#69d", reverse=True, as_cmap=False)
        #palette = sns.color_palette("dark:salmon_r")
        # Plot using the style argument based on the index column (the 4th column in the DataFrame)
        sns.scatterplot(
            data=df,
            x=df.iloc[:, 0], 
            y=df.iloc[:, 1], 
            hue=df["Population"],
            style=df.iloc[:, 3],  # Use the index column directly for style (shapes)
            
            palette=palette, 
            alpha=1.0, 
            ax=ax
        )
        
        ax.get_legend().remove()
        
        # Save the plot with transparent background, but keep the frame
        fig.savefig(f"{out}_new.png", dpi=600, bbox_inches='tight')
        fig.savefig(f"{out}_new_transp.png", dpi=600, transparent=True, bbox_inches='tight')
        plt.close(fig)  # close the figure to release memory
    if color_method == "dG":
        #palette = sns.color_palette("rocket", as_cmap=True)
        #palette = sns.color_palette("icefire", as_cmap=True)
        palette = sns.color_palette("plasma_r", as_cmap=True)
        
        #scaler = MinMaxScaler()

# Reshape the data to a 2D array (required by sklearn)
        #df["normalized_column"] = scaler.fit_transform(df[["dG"]])
        #print(df["normalized_column"])
        # Plot using the style argument based on the index column (the 4th column in the DataFrame)
        sns.scatterplot(
            data=df,
            x=df.iloc[:, 0], 
            y=df.iloc[:, 1], 
            hue=df["dG"],
            style=df.iloc[:, 3],  # Use the index column directly for style (shapes)
            
            palette=palette, 
            alpha=1.0, 
            ax=ax
        )
        
        ax.get_legend().remove()
        
        # Save the plot with transparent background, but keep the frame
        fig.savefig(f"{out}_new.png", dpi=600, bbox_inches='tight')
        fig.savefig(f"{out}_new_transp.png", dpi=600, transparent=True, bbox_inches='tight')
        plt.close(fig)  # close the figure to release memory


def vhse_vector(df):
    df_descriptors = pd.DataFrame()

    for j in df.index:
        data_cleaned_vhse = descriptor_converter.get_descriptors(df["Sequence"][j], descriptor='vhse')
        CDRs_positions = sorted(df["Position"].unique().tolist())
        vhse_cdrs = data_cleaned_vhse[data_cleaned_vhse.index.isin(CDRs_positions)]
        new_list = []
        for i in vhse_cdrs.index:
            new_list.append(vhse_cdrs.loc[i].values.tolist())
        flat_list = [item for sublist in new_list for item in sublist]
        df_temp = pd.DataFrame([flat_list])
       
        # Concatenate the descriptor for each sequence
        df_descriptors = pd.concat([df_descriptors, df_temp], axis=0, ignore_index=True)

    # Reset index and add additional columns
    df_descriptors = df_descriptors.reset_index(drop=True)
    df_descriptors["LogFC_exp"] = df["LogFC_exp"]
    df_descriptors["LogFC_bnd"] = df["LogFC_bnd"]
    df_descriptors["LogFC_bind_adjusted"] = df["LogFC_bind_adjusted"]
    df_descriptors["dG"] = df["dG"]

    return df_descriptors

def vhse_zcales_converter(df, CDR_positions, descriptor):
    df_descriptors = pd.DataFrame()
    
    if descriptor == "vhse":
        for j in df.index:
            data_cleaned_vhse = descriptor_converter.get_descriptors(df["Sequence"][j], descriptor='vhse')
            #CDRs_positions = sorted(df["Position"].unique().tolist())
            vhse_cdrs = data_cleaned_vhse[data_cleaned_vhse.index.isin(CDR_positions)]
            new_list = []
            for i in vhse_cdrs.index:
                new_list.append(vhse_cdrs.loc[i].values.tolist())
            flat_list = [item for sublist in new_list for item in sublist]
            df_temp = pd.DataFrame([flat_list])
       
        # Concatenate the descriptor for each sequence
            df_descriptors = pd.concat([df_descriptors, df_temp], axis=0, ignore_index=True)
    if descriptor == "zcales":
        for j in df.index:
            data_cleaned_vhse = descriptor_converter.get_descriptors(df["Sequence"][j], descriptor='zscales')
            #CDRs_positions = sorted(df["Position"].unique().tolist())
            vhse_cdrs = data_cleaned_vhse[data_cleaned_vhse.index.isin(CDR_positions)]
            new_list = []
            for i in vhse_cdrs.index:
                new_list.append(vhse_cdrs.loc[i].values.tolist())
            flat_list = [item for sublist in new_list for item in sublist]
            df_temp = pd.DataFrame([flat_list])
           
            # Concatenate the descriptor for each sequence
            df_descriptors = pd.concat([df_descriptors, df_temp], axis=0, ignore_index=True)


    # df_descriptors["Sequence"] = df["Sequence"]
    # df_descriptors["dG"] = df["dG"]
    # df_descriptors["Cycle"] = df["Cycle"]
    # df_descriptors["Condition"] = df["Condition"]
    
    #df_descriptors[label] = df[label]

    return df_descriptors



