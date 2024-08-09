#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 03:52:23 2023

@author: lucas
"""

import os

from genetic_algorithm_rosetta import genetic_algo, genetic_algo_sequence

from apt_function import *
import apt_function

from pyrosetta import *
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation

import numpy as np
from numpy.random import uniform
from random import sample
import random
import pandas as pd


############ Running GAPO - Structure

# Initialize PyRosetta
pyrosetta.init()

# Create a scoring function using the "ref2015_cart.wts" weight set
scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")

# Creates pose object from input PDB
starting_pose = pose_from_pdb('2lzt_sem_hetatm_relax_3_times.pdb')
# Relax the starting pose by packing and relaxing it iteratively for 3 times
scorefxn(starting_pose)

# Define a list of single-letter amino acid codes
gene_values = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
starting_pose_seq = [x for x in starting_pose.sequence()]
len(starting_pose_seq)


# Residues to mutate during optimization
fixed_residues = [6, 30, 48, 52, 56, 57, 58, 59, 60, 61, 62, 63, 64, 71, 72, 73, 74, 75, 76, 80, 94, 98, 103, 107, 108, 109, 115, 127]


Natural_Sequences = apt_function.Read_sequence("Natural_lysos.fasta")
init_population = Natural_Sequences

# Initiates GA object
GA = genetic_algo(pose=starting_pose, opt_direction='down',initial_population = init_population, gene_values=gene_values, gene_type='discrete',
              vector_size=len(starting_pose_seq), pop_size=len(init_population), mutation_rate=0.9, segment_fluctuation=0,
              apt_function=apt_rosetta, selection_method='tournament', threads=False,
              convergence_threshold=0, n_cycles=100, tournament_cycles=int(np.round(len(init_population)/4)), tournament_size=4, benchmark=False, 
              lista_fixed=fixed_residues, crossing_over_type='mask', file_name="apt_rosetta_random.txt", cpus  = 25, mutation_type = "random")

# Run the Genetic Algorithm
GA.execute()








