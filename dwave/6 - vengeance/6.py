# ref: https://github.com/prince-ph0en1x/QuASeR/blob/819fcfa85a3a9486ee27a95808f245e51ab1d5de/QA_DeNovoAsb/denovo_009.py

import numpy as np
import math
import matplotlib.pyplot as plt
import networkx as nx
import sys
import time
from collections import defaultdict
import random

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
	for ct in range(0,n_reads**2):
		Q_matrix[ct][ct] += p0
	# Multi-location penalty
	for c in range(0,n_reads):
		for t1 in range(0,n_reads):
			for t2 in range(0,n_reads):
				if t1!=t2:
					Q_matrix[c*n_reads+t1][c*n_reads+t2] += p1
	# Visit repetation penalty
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

"""
Given a list of reads, produces the QUBO matrix
"""
def from_reads_to_QUBO(reads, self_bias=-1.6, repetition=1.6, multilocation=1.6):
	tspAdjM = reads_to_tspAdjM(reads)
	quboAdjM = tspAdjM_to_quboAdjM(tspAdjM, self_bias, repetition, multilocation)
	return quboAdjM_to_quboDict(quboAdjM)

"""
	 --------------------------------- TEST CREATION ---------------------------------
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

"""
	 --------------------------------- SOLUTIONS RECONSTRUCTION ---------------------------------
"""

def get_num_reads(cycle):
	return int(math.sqrt(len(list(cycle.values()))))

# The argument is a dict: {'n0t0': 1, ...}
def is_valid_cycle(cycle):
	num_reads = get_num_reads(cycle)
	if sum(cycle.values()) != num_reads:
		return False

	for n in range(num_reads):
		occ = 0
		for t in range(num_reads):
			occ += cycle['n{}t{}'.format(n, t)]
		if occ != 1:
			return False

	for t in range(num_reads):
		occ = 0
		for n in range(num_reads):
			occ += cycle['n{}t{}'.format(n, t)]
		if occ != 1:
			return False

	return True

# Create a string with the obtained order. Example: '1320'
def transform_cycle_into_str(cycle):
	num_reads = get_num_reads(cycle)
	if not is_valid_cycle(cycle):
		return ''

	# Create a string with the obtained order. Example: '1320'
	result = ''
	for t in range(num_reads):
		for n in range(num_reads):
			if cycle['n{}t{}'.format(n, t)] == 1:
				result += '{}-'.format(str(n))

	# Shift the string so the 0 is on the first position
	double_result = result + result
	i = 0
	while True:
		if double_result[i-1:i+2] == '-0-':
			return double_result[i:i + len(result)]
		i += 1

def get_max_chain_length(embedded_graph):
	return max([ len(x) for _,x in embedded_graph.items()])

"""
	----------------------- Running the experiment ------------------------
"""

"""
Write the solutions in the specific wanted format
"""
def write_solutions(solutions, file_name='results.txt'):
	target = open(file_name, 'w')
	target.write('{\n')

	for sample, energy, ocurrences in solutions:
		frozen_sample  = frozenset(sample.items())
		target.write('{}: ({}, {:.4f}),\n'.format(
			frozen_sample, ocurrences, energy))

	target.write('}')

token = 'DEV-e9ec3f66ba24ec5d5e1739a84289bd7a911823c3'
endpoint = 'https://cloud.dwavesys.com/sapi'

"""
Solve a QUBO  model using D-Wave solver
"""
def solve_qubo_dwave(Q, n_genome_reads, num_reads=100, label=''):
	# Create the solver (connecting to D-Wave) and the Sampler
	config_file='../dwave_adv.conf'
	client = Client.from_config(config_file, profile='ocete')
	solver = client.get_solver('Advantage_system1.1')
	dwsampler = DWaveSampler(solver={'qpu': True, 'topology__type': 'pegasus'},
							 token=token, endpoint=endpoint)

	# We need to embed Q into a valid graph for the D-Wave architecture
	adjdict = edgelist_to_adjacency(solver.edges)
	embed = minorminer.find_embedding(Q, solver.edges)
	Q_embeded = embed_qubo(Q, embed, adjdict)

	# Obtain the response from the solver. This is the actual D-Wave execution!
	start = time.time()
	response_qpt = dwsampler.sample_qubo(Q_embeded, num_reads=num_reads, label=label)
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
			[rebuild_quboDict_from_vector(sol, n_genome_reads), energy, num_appereances])

	return solutions_list, qpu_time, get_max_chain_length(embed)

"""
	----------------------- EVALUATE RESULTS ------------------------
"""

# Given a list of responses [sample, energy, occurences], count the number of valid cycles
def count_valid_responses(responses):
	valid = 0
	for sample, _, occ in responses:
		if is_valid_cycle(sample):
			valid += occ
	return valid

# Returns the number of occurrences that match the actual solution.
# Sol is passed stringified.
def count_matches_responses_solution(responses, sol_str):
	matches = 0
	for sample, _, occ in responses:
		if transform_cycle_into_str(sample) == sol_str:
			matches += occ
	return matches

# Returns '0-1-2-[...]-(num_reads-1)-'
def solution_stringified(num_reads):
	result = ''
	for i in range(num_reads):
		result += '{}-'.format(str(i))

	return result

def compute_energy_best_sol(Q, self_bias=-1.6):
	n_reads = len(Q[0])
	energy = self_bias*n_reads
	for i in range(n_reads):
		energy -= Q[i, (i+1)%n_reads]
	return energy

"""
Solve de novo assembly on D-Wave annealer
"""
def solve_denovo_assembly(reads, self_bias, multilocation, repetition,
			use_QA=False, print_solutions=False, 
			write_solutions=False, filename='results.txt', label=''):
	Q = from_reads_to_QUBO(reads, self_bias, multilocation, repetition)

	max_chain_length = 0
	solutions, qpu_time, max_chain_length = solve_qubo_dwave(Q, num_reads=10000, n_genome_reads=len(reads), label=label)

	# Print some solutions if requested
	if print_solutions:
		print("Minimum Energy Sampled\t===>")
		for sol in solutions[:30]:
			print(sol)

	# Write results if requested
	if write_solutions:
		write_solutions(solutions)
		
	return qpu_time, solutions, max_chain_length

def stress_test(_from=3, _to=15, use_QA=False):

	print('Params\t| Valid cycles\t| Solution reached\t| E_d\t| Best energy\t| Max chain length\t|Time (s)')

	params = [1.6, 3, 5, 10, 20, 50]
	num_reads = 10

	for p in params:
		_, reads = create_test(num_reads=num_reads)
		total_time, responses, max_chain_length = \
			solve_denovo_assembly(reads, -p, p, p, label='Param: {}'.format(p), print_solutions=False, use_QA=use_QA)

		valid_cycles = count_valid_responses(responses)
		min_energy = responses[0][1]

		best_response_stringify = transform_cycle_into_str(responses[0][0])
		sol_stringify = solution_stringified(num_reads)
		best_sol_reached = best_response_stringify == sol_stringify
		n_sol_reached = count_matches_responses_solution(responses, sol_stringify)

		Q = reads_to_tspAdjM(reads)
		best_energy = compute_energy_best_sol(Q)
		e_delta = best_energy - min_energy

		print('({}, {}, {})\t| {}\t\t| {}\t\t\t| {:.6f}\t| {:.6f}\t| {}\t\t| {:.6f}'.format(
				str(-p), p, p, valid_cycles, n_sol_reached, e_delta, best_energy, max_chain_length, total_time))

"""
EXPERIMENT 5
"""

# Call for SQA
# stress_test()

stress_test(10, 11, use_QA=True)
