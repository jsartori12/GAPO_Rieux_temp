#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 14:20:19 2024

@author: joao
"""


import PETases_sequences_alignment
from Bio import SeqIO
import pandas as pd
from Bio import Align
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import pandas as pd
import time
from Bio.Align import substitution_matrices
import PETases_sequences_alignment
from PETases_sequences_alignment import filter_sequences
from PETases_sequences_alignment import generate_synthetic_sequence
from Bio.Phylo.TreeConstruction import DistanceCalculator



# Replace 'your_file.fasta' with the path to your FASTA file
fasta_file = "uniprotkb_lysozyme_AND_reviewed_true_2024_08_09.fasta"

# Initialize lists to store sequence data
ids = []
sequences = []
descriptions = []

# Read the FASTA file and extract data
for record in SeqIO.parse(fasta_file, "fasta"):
    ids.append(record.id)
    sequences.append(str(record.seq))
    descriptions.append(record.description)

# Create a DataFrame
df = pd.DataFrame({
    "ID": ids,
    "Sequence": sequences,
    "Description": descriptions
})



new_id = "reference"
new_sequence = "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
new_description = "reference sequence from PDB:2lzt"

fixed_residues = [5, 29, 47, 51, 55, 56, 57, 58, 59, 60, 61, 62, 63, 70, 71, 72, 73, 74, 75, 79, 93, 97, 102, 106, 107, 108, 114, 126]


# Create a DataFrame for the new sequence
new_row = pd.DataFrame({
    "ID": [new_id],
    "Sequence": [new_sequence],
    "Description": [new_description]
})

# Concatenate the new row with the original DataFrame, placing it at the top
df = pd.concat([new_row, df]).reset_index(drop=True)

synthetic_sequences = PETases_sequences_alignment.generate_synthetic_sequence(df, fixed_residues)

# Create an alignment file in FASTA format for synthetic sequences
alignment_file = "aligned_sequences.fasta"
with open(alignment_file, "w") as f:
    for i in synthetic_sequences.index:
        f.write(f">{synthetic_sequences['ID'].iloc[i]}\n{synthetic_sequences['Sequence'].iloc[i]}\n")

alignment = AlignIO.read(alignment_file, "fasta")

# Create a distance calculator
calculator = DistanceCalculator('blosum62')
dm = calculator.get_distance(alignment)


df = pd.DataFrame(dm.matrix)

# Define new column names
new_column_names = dm.names
# Rename columns
df.rename(columns=dict(zip(df.columns, new_column_names)), inplace=True)
df.rename(index=dict(zip(df.index, new_column_names)), inplace=True)

df.to_csv("distance_matrix_synthethic_seqs.csv")

dist_matrix = pd.read_csv("distance_matrix_synthethic_seqs.csv")
new_column_names = dist_matrix.iloc[:,0]
dist_matrix.rename(index=dict(zip(dist_matrix.index, new_column_names)), inplace=True)
dist_matrix = dist_matrix.iloc[:,1:len(dist_matrix.columns)]


highpop, highdist = PETases_sequences_alignment.sample_sequences(df = dist_matrix, sample_size = 200, n_samples = 10000)   

sample_natural = synthetic_sequences[synthetic_sequences["ID"].isin(highpop[0])] 

PETases_sequences_alignment.create_multifasta(df = sample_natural)
