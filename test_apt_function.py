#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 00:09:57 2024

@author: joao
"""

#### ESM stuff
import torch
import esm
import random
import math
import argparse

# Load ESM-2 model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results

def apt_esm(seq, temperature=1):
    data = [
        ("protein1", seq)
    ]

    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

    # Predict token probabilities
    with torch.no_grad():
        token_probs = model(batch_tokens, repr_layers=[33])["logits"]

    # Apply temperature
    token_probs /= temperature

    softmax = torch.nn.Softmax(dim=-1)
    probabilities = softmax(token_probs)

    # Calculate the average probability of the tokens at each position, excluding start and end tokens
    average_probs = []
    for i in range(1, batch_lens[0] - 1):  # Exclude start and end tokens
        token_idx = batch_tokens[0, i].item()
        token_prob = probabilities[0, i, token_idx].item()
        average_probs.append(token_prob)

    avg_probability = sum(average_probs) / len(average_probs)
    
    return avg_probability

parser = argparse.ArgumentParser(description="Calculate average token probabilities using ESM model.")
parser.add_argument("--seq", type=str, help="Protein sequence")
parser.add_argument("--temperature", type=float, default=1.0, help="Temperature for scaling probabilities")
args = parser.parse_args()


avg_probability = apt_esm(args.seq, args.temperature)
print(f"Average Probability: {avg_probability}")

# Optionally save the average probability to a file
with open("avg_probability.txt", "w") as f:
    f.write(str(avg_probability))