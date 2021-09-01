# ref: https://github.com/prince-ph0en1x/QuASeR/blob/819fcfa85a3a9486ee27a95808f245e51ab1d5de/QA_DeNovoAsb/denovo_009.py

import numpy as np
import math
import matplotlib.pyplot as plt
import networkx as nx
import sys
import time
import itertools
import random
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
Convert adjacency matrix of pair-wise overlap for TSP to QUBO matrix of TSP
"""
def tspAdjM_to_quboAdjM(tspAdjM, p0, p1, p2):
	n_reads = len(tspAdjM)
	# Initialize
	Q_matrix = np.zeros((n_reads**2,n_reads**2)) # Qubit index semantics: {c(0)t(0) |..| c(i)-t(j) | c(i)t(j+1) |..| c(i)t(n-1) | c(i+1)t(0) |..| c(n-1)t(n-1)}
	# Assignment reward (self-bias)
	p0 = -1.6
	for ct in range(0,n_reads**2):
		Q_matrix[ct][ct] += p0
	# Multi-location penalty
	p1 = -p0 # fixed emperically by trail-and-error
	for c in range(0,n_reads):
		for t1 in range(0,n_reads):
			for t2 in range(0,n_reads):
				if t1!=t2:
					Q_matrix[c*n_reads+t1][c*n_reads+t2] += p1
	# Visit repetation penalty
	p2 = p1
	for t in range(0,n_reads):
		for c1 in range(0,n_reads):
			for c2 in range(0,n_reads):
				if c1!=c2:
					Q_matrix[c1*n_reads+t][c2*n_reads+t] += p2
	# Path cost
	# kron of tspAdjM and a shifted diagonal matrix
	for ci in range(0,n_reads):
		for cj in range(0,n_reads):
			for ti in range(0,n_reads):
				tj = (ti+1)%n_reads
				Q_matrix[ci*n_reads+ti][cj*n_reads+tj] += -tspAdjM[ci][cj]
	return Q_matrix

"""
Convert QUBO matrix of TSP to QUBO dictionary of weighted adjacency list
"""
def quboAdjM_to_quboDict(Q_matrix):
	n_reads = int(math.sqrt(len(Q_matrix)))
	Q = {}
	for i in range(0,n_reads**2):
		ni = 'n'+str(int(i/n_reads))+'t'+str(int(i%n_reads))
		for j in range(0,n_reads**2):
			nj = 'n'+str(int(j/n_reads))+'t'+str(int(j%n_reads))
			if Q_matrix[i][j] != 0:
				Q[(ni,nj)] = Q_matrix[i][j]
	return Q

"""
Given a solution vector returned by DWave, rebuilds the associated quboDict.

Example: In a 2 nodes problem, the vector:
[0, 1, 1, 0]
Will be transform into:
{'n0t0': 0, 'n0t1': 1, 'n1t0': 1, 'n1t1': 0}
"""
def rebuild_quboDict_from_vector(v, num_reads):
	result = {}
	for k, val in enumerate(v):
		n = int(k / num_reads)
		t = k % num_reads
		result['n{}t{}'.format(n, t)] = val
	return result

# ----------------------- Plotting ------------------------

"""
For experimenting with embedding the QUBO model in Chimera graph
"""
def embed_qubo_chimera(Q_matrix, plotIt = False):
	connectivity_structure = dnx.chimera_graph(3,3) # try to minimize
	G = nx.from_numpy_matrix(Q_matrix)
	max_chain_length = 0
	while(max_chain_length == 0):
		embedded_graph = minorminer.find_embedding(G.edges(), connectivity_structure)
		for _, chain in embedded_graph.items():
		    if len(chain) > max_chain_length:
		        max_chain_length = len(chain)
	print("max_chain_length", max_chain_length) # try to minimize
	if plotIt:
		dnx.draw_chimera_embedding(connectivity_structure, embedded_graph)
		plt.show()

"""
For experimenting with embedding the QUBO model in Chimera graph
"""
def embed_qubo_pegasus(Q_matrix, plotIt = False):
	connectivity_structure = dnx.pegasus_graph(15) # try to minimize
	G = nx.from_numpy_matrix(Q_matrix)
	max_chain_length = 0
	while(max_chain_length == 0):
		embedded_graph = minorminer.find_embedding(G.edges(), connectivity_structure)
		for _, chain in embedded_graph.items():
		    if len(chain) > max_chain_length:
		        max_chain_length = len(chain)
	print("max_chain_length", max_chain_length) # try to minimize
	if plotIt:
		dnx.draw_chimera_embedding(connectivity_structure, embedded_graph)
		plt.show()

"""
	 --------------------------------- TEST MANAGEMENT ---------------------------------
"""

def chop_solution(sol, read_length, fixed_overlap):
	reads = []
	init_index = 0
	step = read_length - fixed_overlap
	reads.append( sol[:init_index + read_length] )
	while True:
		init_index += step
		reads.append( sol[init_index : init_index + read_length] )

		if init_index + read_length >= len(sol):
			break

	return reads

"""
Variation from: https://docs.python.org/2/library/itertools.html#recipes
Returns a random cartesian product
"""
def random_product(iter, length):
    return ''.join(random.choice(iter) for _ in range(length))

def create_test(num_reads=4, read_length=150, fixed_overlap=50):
	length = (read_length - fixed_overlap) *(num_reads-1) + read_length
	sol = random_product('ACGT', length)
	reads = chop_solution(sol, read_length, fixed_overlap)
	return (sol, reads)

# ----------------------- Running the experiment ------------------------

"""
Solve de novo assembly on D-Wave annealer
"""
def try_embedding(n_reads, solver, print_embedding=True):
	start_total = time.time()

	# Create test
	start_test = time.time()
	_, reads = create_test(num_reads=n_reads)
	test_creation_time = time.time() - start_test

	# Prepare the data
	tspAdjM = reads_to_tspAdjM(reads)
	quboAdjM = tspAdjM_to_quboAdjM(tspAdjM, -1.6, 1.6, 1.6)
	Q = quboAdjM_to_quboDict(quboAdjM)

	# Try to embed Q into a valid graph for the pegasus topology
	start_embedding = time.time()
	edgelist = solver.edges
	adjdict = edgelist_to_adjacency(edgelist)
	embed = minorminer.find_embedding(Q, edgelist)
	Q_embeded = embed_qubo(Q, embed, adjdict)
	embedding_time = time.time() - start_embedding

	total_time = time.time() - start_total

	print('Successful: {} reads, {} nodes graph, {:.6f}s test creation time, {:.6f}s embedding time, {:.6f}s total time'.format(
		len(reads), len(reads)*len(reads), test_creation_time, embedding_time, total_time))

	#print('Embedding of graph with {} reads ({} nodes graph) successful'.format(
	#			len(reads), len(reads)*len(reads)))
	

"""
EXPERIMENT 4
"""
def test_limits():

	n_reads = 2
	while True:
		# Create the solver (connecting to D-Wave) and the Sampler
		config_file='../dwave.conf'
		client = Client.from_config(config_file, profile='ocete')
		solver = client.get_solver('Advantage_system1.1') # Use Advantage to test Pegasus

		try_embedding(n_reads, solver, print_embedding=False)
		n_reads += 1

	client.close()

test_limits()