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


############ Running GAPO - Sequence


gene_values = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']


lysozime = "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"

fixed_residues = [6, 30, 48, 52, 56, 57, 58, 59, 60, 61, 62, 63, 64, 71, 72, 73, 74, 75, 76, 80, 94, 98, 103, 107, 108, 109, 115, 127]

Natural_Sequences = apt_function.Read_sequence("Natural_lysos.fasta")
init_population = Natural_Sequences

# Initiates GA object
GA = genetic_algo_sequence(opt_direction='up',initial_population = init_population, gene_values=gene_values, gene_type='discrete',
             vector_size=len(init_population[0]), pop_size=len(init_population), mutation_rate=0.9, segment_fluctuation=0,
             apt_function=apt_esm, selection_method='tournament', threads=False,
             convergence_threshold=0, n_cycles=100, lista_fixed = fixed_residues, tournament_cycles=int(np.round(len(init_population)/4)), tournament_size=4, benchmark=False,  
             crossing_over_type='mask', file_name="apt_esm_random.txt", mutation_type = "random")

# Run the Genetic Algorithm
GA.execute()

