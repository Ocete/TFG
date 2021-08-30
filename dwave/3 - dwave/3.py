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

# ----------------------- Running the experiment ------------------------

"""
Solve a QUBO model using dimod exact solver
"""
def solve_qubo_simulated(Q, all=False, print_it=False, save_it=False, num_reads=10000):
	#solver = dimod.SimulatedAnnealingSampler()
	solver = neal.SimulatedAnnealingSampler()
	
	start = time.time()
	response = solver.sample_qubo(Q, num_reads=num_reads)
	qpu_time = time.time() - start

	# Count the number of ocurrences
	ocurrences = defaultdict(int)
	for sample, energy in response.data(['sample', 'energy']):
		frozen_sample  = frozenset(sample.items())
		ocurrences[frozen_sample] += 1

	# Print the results
	if print_it:
		printed = defaultdict(bool)
		for sample, energy in response.data(['sample', 'energy']):
			frozen_sample  = frozenset(sample.items())
			if not printed[frozen_sample]:
				printed[frozen_sample] = True
				print('({}, {}) --> {:.4f}'.format(sample, ocurrences[frozen_sample], energy))

	if save_it:
		target = open('results.txt', 'w')
		target.write('{\n')
		printed_2 = defaultdict(bool)

		for sample, energy in response.data(['sample', 'energy']):
			frozen_sample  = frozenset(sample.items())
			if not printed_2[frozen_sample]:
				printed_2[frozen_sample] = True
				target.write('{}: ({}, {:.4f}),\n'.format(frozen_sample, ocurrences[frozen_sample], energy))

		target.write('}')

	return qpu_time

def compute_mean_time(repeat=10):
	mean = 0.0
	for _ in range(repeat):
		mean += deNovo_on_DWave()
	print('Mean time in {} executions: {}'.format(repeat, mean/repeat))


"""
Solve a QUBO  model using D-Wave solver
"""
def solve_qubo_dwave(Q, num_reads=100):
	# Create the solver (connecting to D-Wave) and the Sampler
	config_file='../dwave.conf'
	client = Client.from_config(config_file, profile='ocete')
	solver = client.get_solver() # Available QPUs: DW_2000Q_2_1 (2038 qubits), DW_2000Q_5 (2030 qubits)
	dwsampler = DWaveSampler(config_file=config_file)

	# We need to embed Q into a valid graph for the D-Wave architecture
	edgelist = solver.edges
	adjdict = edgelist_to_adjacency(edgelist)
	embed = minorminer.find_embedding(Q, edgelist)
	Q_embeded = embed_qubo(Q, embed, adjdict)

	# Obtain the response from the solver. This is the actual D-Wave execution!
	start = time.time()
	response_qpt = dwsampler.sample_qubo(Q_embeded, num_reads=num_reads)
	qpu_time = time.time() - start
	client.close()

	# Transform the response from the embeded graph to our original architecture
	bqm = dimod.BinaryQuadraticModel.from_qubo(Q)
	unembedded = unembed_sampleset(response_qpt, embed, bqm, chain_break_method=majority_vote)

	# Sort the solutions from lowest energy and format them to quboDict format
	unformated_solutions_list = sorted(unembedded.record, key=lambda x: +x[1])
	solutions_list = []
	for sol, energy, num_appereances in unformated_solutions_list:
		solutions_list.append(
			tuple([rebuild_quboDict_from_vector(sol, len(reads)), energy, num_appereances]))

	return solutions_list, qpu_time

"""
Write the solutions in the specific wanted format
"""
def write_solutions(solutions):
	target = open('results.txt', 'w')
	target.write('{\n')

	for sample, energy, ocurrences in solutions:
		frozen_sample  = frozenset(sample.items())
		target.write('{}: ({}, {:.4f}),\n'.format(
			frozen_sample, ocurrences, energy))

	target.write('}')

"""
Solve de novo assembly on D-Wave annealer
"""
def deNovo_on_DWave_QA(reads, print_solutions=False, write_solutions=False):
	tspAdjM = reads_to_tspAdjM(reads)
	#print('tspAdjM: {}'.format(tspAdjM))

	quboAdjM = tspAdjM_to_quboAdjM(tspAdjM, -1.5, 1.5, 1.5) # self-bias, multi-location, repetition
	#print('quboAdjM: {}'.format(quboAdjM))

	quboDict = quboAdjM_to_quboDict(quboAdjM)
	#print('quboDict: {}'.format(quboDict))

	solutions, qpu_time = solve_qubo_dwave(quboDict)
	#print('solutions_list: {}'.format(solutions_list))

	# Print some solutions if requested
	if print_solutions:
		print("Minimum Energy Configurations from D-Wave\t===>")
		for sol in solutions[:30]:
			print(sol)

	if write_solutions:
		write_solutions(solutions)
		
	return qpu_time

"""
EXPERIMENT 3
"""

reads = ['ATGGCGTGCA', 'GCGTGCAATG', 'TGCAATGGCG', 'AATGGCGTGC']
deNovo_on_DWave_QA(reads, print_solutions=True)

# compute_mean_time()
