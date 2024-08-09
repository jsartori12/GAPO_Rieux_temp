#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 11:42:53 2023

@author: joao
"""

from Bio import Align
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import pandas as pd
import time
from Bio.Align import substitution_matrices

# Load your sequences from a file
# sequences = SeqIO.parse("uniprotkb_family_AB_hydrolase_superfami_2023_11_16.fasta", "fasta")

# # Create a list of SeqRecord objects
# seq_records = [record for record in sequences]

def pre_filter_sequences(input_string):
    """
    Pre-filter sequences based on valid protein characters.

    This function checks if a given input string contains only characters
    that are valid in a protein sequence. It uses a predefined set of valid
    protein characters for the check.

    Parameters:
    - input_string (str): The input string representing a protein sequence.

    Returns:
    - bool: True if the input string contains only valid protein characters,
      False otherwise.
    """

    # Define a list of valid protein characters
    valid_protein_characters = list("ACDEFGHIKLMNPQRSTVWY")

    # Check if all characters in the input string are valid protein characters
    return all(char in valid_protein_characters for char in input_string)


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
    #alignments = list(alignment[0])
    alignments = alignment[0]

    return alignments

def calculate_identity_cover(seq_ref, query_seq):
    """
    Calculate sequence identity and query coverage between two sequences.

    Parameters:
    - seq_ref (str): Reference sequence for alignment.
    - query_seq (str): Query sequence for alignment.

    Returns:
    - identity (float): Sequence identity percentage.
    - coverage (float): Query coverage percentage.
    """
    # Aligning sequences
    alignment = algn_pairwise(seq1 = seq_ref, seq2 = query_seq)
    
    # Extract alignment information
    matches = compute_identical_matches(align1 = alignment[0], align2 = alignment[1])
    sequence_length = len(seq_ref)
    
    # Calculate sequence identity query vs reference
    identity = (matches / sequence_length) * 100
    
    # Calculate query cover
    query_sequence = alignment[1]
    aligned_query = alignment[1]
    coverage = sum(1 for q in aligned_query if q != "-") / len(query_sequence) * 100

    return round(identity, 2), round(coverage, 2)

def filter_sequences(sequences_record, id_cutoff, cover_cutoff):
    """
    Filter sequences based on identity and query coverage thresholds.

    Parameters:
    - sequences_record (list): List of Bio.SeqRecord objects.
    - id_cutoff (float): Identity threshold for filtering.
    - cover_cutoff (float): Query coverage threshold for filtering.

    Returns:
    - df_all_seqs (pd.DataFrame): Filtered sequences with Seq_ID, Sequence, Identity, and Query_Cover columns.
    """
    start_time = time.time()  # Record the start time
    
    # Create dataframe with Seq_ID | Sequence | Identity | Query Cover
    df_all_seqs = pd.DataFrame(columns=range(1, 5), index=range(1, 2))
    
    # Define new column names
    new_column_names = ["Seq_ID", "Sequence", "Identity", "Query_Cover"]

    # Rename columns
    df_all_seqs.rename(columns=dict(zip(df_all_seqs.columns, new_column_names)), inplace=True)
    df_dummy = df_all_seqs.copy()

    seq_ref = str(sequences_record[0].seq)

    for i in range(0, len(sequences_record)):
        query_seq = str(sequences_record[i].seq)
        if pre_filter_sequences(query_seq):
            # Calculate identity and coverage for the current sequence
            #print(query_seq)
            #print(sequences_record[i].id)
            id_seq, cover_seq = calculate_identity_cover(seq_ref=seq_ref, query_seq=query_seq)
            
            # Check if the sequence passes the filtering criteria
            if id_seq > id_cutoff and id_seq < 100 and cover_seq > cover_cutoff:
                # Create a temporary DataFrame for the current sequence
                df_temp = df_dummy.copy()
                df_temp["Seq_ID"][1] = sequences_record[i].id
                df_temp["Sequence"][1] = query_seq
                df_temp["Identity"][1] = id_seq
                df_temp["Query_Cover"][1] = cover_seq
                
                # Concatenate the temporary DataFrame to the main DataFrame
                df_all_seqs = pd.concat([df_all_seqs, df_temp], ignore_index=True)
            load = i / len(sequences_record) * 100
        
        
        print(f"{round(load,2)}% done... \U0001f600")
        
    df_all_seqs["Seq_ID"][0] = "PETase"
    df_all_seqs["Sequence"][0] = 'KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL'

    end_time = time.time()  # Record the end time
    execution_time = end_time - start_time  # Calculate the execution time
    
    df_all_seqs.to_csv("Results.csv")
    
    # Write the results to a text file
    with open("output.txt", "w") as f:
        f.write(f"Execution Time: {round(execution_time/60, 1)} minutes\n")
                
    return df_all_seqs

def compute_identical_matches(align1, align2):
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
    
    # Iterate over each pair of aligned residues
    for j, (res1, res2) in enumerate(zip(align1, align2)):
        # If neither residue is a gap ('-') and the residues are identical, update the count
        if res1 != "-" and res2 != "-" and res1 == res2:
            compute_total += 1
    
    return compute_total

def generate_synthetic_sequence(df_filtered, fixed_indexes):
    """
    Generates a synthetic sequence by aligning each sequence in the DataFrame with a reference sequence,
    while keeping specified positions in the reference sequence unchanged.

    Parameters:
    - df_filtered (pd.DataFrame): DataFrame containing a column named "Sequence" with sequences to be aligned.
    - fixed_indexes (list of int): List of indexes in the reference sequence that should remain unchanged.

    Returns:
    - pd.DataFrame: DataFrame with synthetic sequences.

    Note: Requires the align_pairwise function to be defined elsewhere.
    """
    # Select the reference sequence from the first row of the DataFrame
    reference_sequence = df_filtered["Sequence"].iloc[0]

    # Iterate over each row in the DataFrame starting from the second row
    for i in range(1, len(df_filtered)):
        # Align the reference sequence with the current sequence in the DataFrame
        aligned = algn_pairwise(reference_sequence, df_filtered["Sequence"].iloc[i])

        # Extract aligned sequences from the alignment result
        refseq = aligned[0]
        aligned_sequence = aligned[1]

        # Initialize a list with characters from the reference sequence
        synthetic_sequence = list(refseq)

        # Iterate over each pair of aligned residues
        for j, (res1, res2) in enumerate(zip(refseq, aligned_sequence)):
            # Skip positions that are in the fixed_indexes list
            if j in fixed_indexes:
                continue
            
            # If neither residue is a gap ('-'), update the synthetic sequence
            if res1 != "-" and res2 != "-":
                synthetic_sequence[j] = aligned_sequence[j]

        # Remove gaps from the synthetic sequence
        synthetic_sequence_nogap = "".join([char for char in synthetic_sequence if char != "-"])

        # Update the "Sequence" column in the DataFrame with the synthetic sequence
        df_filtered.at[i, "Sequence"] = synthetic_sequence_nogap

    # Return the DataFrame with synthetic sequences
    return df_filtered

def calculate_total_distance(df):
    sampled_df = df
    total_distance = sampled_df.sum().sum()
    
    return total_distance    

def sample_sequences(df, sample_size, n_samples):
    
    samples_list = []
    distances_list = []
    
    for i in range(0, n_samples):
        
        
        samples_list_temp = []
        
        sampled_columns = df.sample(n=sample_size, axis=1)
        sampled_columns_names = sampled_columns.columns.tolist()
        
        samples_list_temp.append(sampled_columns_names)
        dist_temp = calculate_total_distance(df = sampled_columns )
        distances_list.append(dist_temp)
        samples_list.append(samples_list_temp)
    
    index_higher_dist = distances_list.index(max(distances_list))
    return samples_list[index_higher_dist], distances_list[index_higher_dist]
def create_multifasta(df):
    alignment_file = "seqs.fasta"
    with open(alignment_file, "w") as f:
        for i in range(0, len(df)):
            f.write(f">{df['ID'].iloc[i]}\n{df['Sequence'].iloc[i]}\n")
