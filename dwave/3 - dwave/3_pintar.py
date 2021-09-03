# ref: https://github.com/prince-ph0en1x/QuASeR/blob/819fcfa85a3a9486ee27a95808f245e51ab1d5de/QA_DeNovoAsb/denovo_009.py

import numpy as np
import math
import matplotlib.pyplot as plt
import networkx as nx
import sys
import time
from collections import defaultdict

import dimod
import neal
import minorminer
import dwave_networkx as dnx
from dwave.cloud import Client
from dwave.embedding import embed_qubo, unembed_sampleset
from dwave.embedding.utils import edgelist_to_adjacency
from dwave.system.samplers import DWaveSampler
from dwave.embedding.chain_breaks import majority_vote

# ----------------------- Preparing the data ------------------------

"""
Overlap between pair-wise reads
"""
def align(read1, read2, mm):
	l1 = len(read1)
	l2 = len(read2)
	for shift in range(l1-l2,l1):
		mmr = 0
		r2i = 0
		for r1i in range(shift,l1):
			if read1[r1i] != read2[r2i]:
				mmr += 1
			r2i += 1
			if mmr > mm:
				break
		if mmr <= mm:
			return l2-shift
	return 0

"""
Convert set of reads to adjacency matrix of pair-wise overlap for TSP
"""
def reads_to_tspAdjM(reads, max_mismatch = 0):
	n_reads = len(reads)
	O_matrix = np.zeros((n_reads,n_reads)) # edge directivity = (row id, col id)
	for r1 in range(0,n_reads):
		for r2 in range(0,n_reads):
			if r1!=r2:
				O_matrix[r1][r2] = align(reads[r1],reads[r2],max_mismatch)
	O_matrix = O_matrix / np.linalg.norm(O_matrix)
	return O_matrix

"""
For experimenting with embedding the QUBO model in Chimera graph
"""
def embed_qubo_chimera(Q_matrix, plt_show = True, filename=''):
	connectivity_structure = dnx.chimera_graph(3,3) # try to minimize
	G = nx.from_numpy_matrix(Q_matrix)
	max_chain_length = 0
	while(max_chain_length == 0):
		embedded_graph = minorminer.find_embedding(G.edges(), connectivity_structure)
		for _, chain in embedded_graph.items():
		    if len(chain) > max_chain_length:
		        max_chain_length = len(chain)
	# print("max_chain_length",max_chain_length) # try to minimize
	
	dnx.draw_chimera_embedding(connectivity_structure, embedded_graph)
	if plt_show:
		plt.show()
	if filename != '':
		plt.savefig('../../thesis/figures/experiments/{}'.format(filename),
					bbox_inches='tight')

# ----------------------- Experiment ------------------------

"""
Solve de novo assembly on D-Wave annealer
"""
def plot_chimera_deNovo(reads, plt_show=False, filename=''):
	tspAdjM = reads_to_tspAdjM(reads)
	embed_qubo_chimera(tspAdjM, plt_show=plt_show, filename=filename)


"""
EXPERIMENT 3
"""

reads = ['ATGGCGTGCA', 'GCGTGCAATG', 'TGCAATGGCG', 'AATGGCGTGC']
plot_chimera_deNovo(reads, filename='exp2_4reads.png')

reads = ['AAA' for _ in range(13)]
plot_chimera_deNovo(reads, plt_show=False, filename='exp2_{}reads.png'.format(len(reads)))

