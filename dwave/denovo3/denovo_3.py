# ref: https://github.com/prince-ph0en1x/QuASeR/blob/819fcfa85a3a9486ee27a95808f245e51ab1d5de/QA_DeNovoAsb/denovo_009.py

"""
Explore embedding
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import networkx as nx
import sys
import itertools
import random
import time

import dimod
import minorminer
# import dwave_networkx as dnx
from dwave.cloud import Client
from dwave.embedding import embed_qubo, embed_ising, unembed_sampleset
from dwave.embedding.utils import edgelist_to_adjacency
from dwave.system.samplers import DWaveSampler
from dwave.embedding.chain_breaks import majority_vote

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

def create_test(length, read_length, fixed_overlap):
	sol = random_product('ACGT', length)
	reads = chop_solution(sol, read_length, fixed_overlap)
	return (sol, reads)

def print_test_results(test_number, expected_sol, sol, print_all=False):
	sucess = expected_sol == sol
	if not sucess:
		print('Test {} results: Failure \nExpected: {} \nObtained: {}'.format(
				test_number, expected_sol, sol))
	elif print_all:
		print('Test {} results: Sucess \nExpected: {} \nObtained: {}'.format(
				test_number, expected_sol, sol))
	else:
		print('Test {} results: Sucess'.format(test_number))

"""
	 --------------------------------- BUILDING THE GEN BACK FROM THE RESPONSE ---------------------------------
"""

"""
Retrieves from the response the node in instante t,
given the number of reads N.
"""
def get_node_in_instant(t, response, N):
	for i in range(N):
		if response['n{}t{}'.format(i, t)] == 1:
			return i
	print('Error building gen in instant {}'.format(t))


def	stitch(str1, str2):
	l1 = len(str1)
	for i in range(len(str2), 0, -1):
		if str1[ l1-i :] == str2[:i]:
			return str1 + str2[i:]

	print('Error in stitching: {} and {}'.format(str1, str2))
	return str1


"""
Given the response and the reads,
returns the obtained response.
"""
def rebuild_gen_from_response(response, reads):
	if not response:
		print('Error rebuild_gen_from_response: response is empty: {} --endline--'.format(
			response))
		return ''

	N = len(reads)
	reads_order = []
	for i in range(N):
		t = get_node_in_instant(i, response, N)
		reads_order.append(reads[i])

	result = reads_order[0]
	for str2 in reads_order[1:]:
		result = stitch(result, str2)
	return result
	
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
	--------------------------------- DWAVE STUFF ---------------------------------
"""

"""
Read *.qubo file to form the QUBO model
"""
def quboFile_to_quboDict(filename):
	f = open(filename, "r")
	qubo_header = f.readline().split()
	Q = {}
	for i in range(0,int(qubo_header[4])):
		x = f.readline().split()
		Q[(x[0],x[1])] = float(x[2])
	for i in range(0,int(qubo_header[5])):
		x = f.readline().split()
		Q[(x[0],x[1])] = float(x[2])
	f.close()
	return Q

"""
Overlap between pair-wise reads
"""
def align(read1,read2,mm):
	l1 = len(read1)
	l2 = len(read2)
	for shift in range(l1-l2,l1):
		mmr = 0
		r2i = 0
		for r1i in range(shift,l1):
			if read1[r1i]!=read2[r2i]:
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
	for i in range(0, n_reads**2):
		ni = 'n'+str(int(i/n_reads))+'t'+str(int(i%n_reads))
		for j in range(0,n_reads**2):
			nj = 'n'+str(int(j/n_reads))+'t'+str(int(j%n_reads))
			if Q_matrix[i][j] != 0:
				Q[(ni,nj)] = Q_matrix[i][j]
	return Q

"""
Solve a QUBO model using dimod exact solver
"""
def solve_qubo_exact(Q, print_solutions=False, all=False):
	solver = dimod.ExactSolver()
	response = solver.sample_qubo(Q)
	minE = min(response.data(['sample', 'energy']), key=lambda x: x[1])

	sol = {}
	for sample, energy in response.data(['sample', 'energy']):
		if energy == minE[1]:
			sol = sample

			if print_solutions:
				print('{} --> {}'.format(sample, energy))
		elif print_solutions and all:
			print('{} --> {}'.format(sample, energy))

	return sol

"""
Solve an Ising model using dimod exact solver
"""
"""
def solve_ising_exact(hii,Jij, plotIt = False):
	solver = dimod.ExactSolver()
	response = solver.sample_ising(hii,Jij)
	print("Minimum Energy Configurations\t===>")
	minE = min(response.data(['sample', 'energy']), key=lambda x: x[1])
	for sample, energy in response.data(['sample', 'energy']):
		if energy == minE[1]:
			print(sample,energy)
	if plotIt:
		y = []
		for sample, energy in response.data(['sample', 'energy']): y.append(energy)
		plt.plot(y)
		plt.xlabel('Solution landscape')
		plt.ylabel('Energy')
		plt.show()
"""

"""
Solve an Ising model using D-Wave solver
"""
"""
def solve_ising_dwave(hii,Jij):
	config_file='../dwave.conf'
	client = Client.from_config(config_file, profile='ocete')
	solver = client.get_solver() # Available QPUs: DW_2000Q_2_1 (2038 qubits), DW_2000Q_5 (2030 qubits)
	dwsampler = DWaveSampler(config_file=config_file)

	edgelist = solver.edges
	adjdict = edgelist_to_adjacency(edgelist)
	embed = minorminer.find_embedding(Jij.keys(),edgelist)
	[h_qpu, j_qpu] = embed_ising(hii, Jij, embed, adjdict)

	response_qpt = dwsampler.sample_ising(h_qpu, j_qpu, num_reads=solver.max_num_reads())
	client.close()

	bqm = dimod.BinaryQuadraticModel.from_ising(hii, Jij)
	unembedded = unembed_sampleset(response_qpt, embed, bqm, chain_break_method=majority_vote)
	print("Maximum Sampled Configurations from D-Wave\t===>")
	solnsMaxSample = sorted(unembedded.record,key=lambda x: -x[2])
	for i in range(0,10):
		print(solnsMaxSample[i])
	print("Minimum Energy Configurations from D-Wave\t===>")
	solnsMinEnergy = sorted(unembedded.record,key=lambda x: +x[1])
	for i in range(0,10):
		print(solnsMinEnergy[i])
"""

"""
Solve a QUBO  model using D-Wave solver
"""
def solve_qubo_dwave(Q, num_reads):
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
	response_qpt = dwsampler.sample_qubo(Q_embeded, num_reads=10)
	client.close()

	# Transform the response from the embeded graph to our original architecture
	bqm = dimod.BinaryQuadraticModel.from_qubo(Q)
	unembedded = unembed_sampleset(response_qpt, embed, bqm, chain_break_method=majority_vote)

	# Order the solutions by lower energy and format the solutions to quboDict format
	solutions_list = sorted(unembedded.record, key=lambda x: +x[1])
	formated_sol_list = []
	for sol, energy, num_appereances in solutions_list:
		formated_sol_list = (rebuild_quboDict_from_vector(sol, num_reads), energy, num_appereances)
	return formated_sol_list


"""
Solve de novo assembly on D-Wave annealer
"""
def deNovo_on_DWave_QA(reads, print_solutions=False):
	tspAdjM = reads_to_tspAdjM(reads)
	#print('tspAdjM: {}'.format(tspAdjM))

	quboAdjM = tspAdjM_to_quboAdjM(tspAdjM, -1.6, 1.6, 1.6) # self-bias, multi-location, repetation
	#print('quboAdjM: {}'.format(quboAdjM))

	quboDict = quboAdjM_to_quboDict(quboAdjM)
	#print('quboDict: {}'.format(quboDict))

	solutions_list = solve_qubo_dwave(quboDict, len(reads))
	print('solutions_list: {}'.format(solutions_list))

	gen = rebuild_gen_from_response(solutions_list[0], reads)
	return gen

"""
	--------------------------------- EXPERIMENTS ---------------------------------
"""

def run_test_1():
	expected_sol = 'ATGGCGTGCAATGGCGTGC'
	reads = ['ATGGCGTGCA','GCGTGCAATG','TGCAATGGCG','AATGGCGTGC']
	reads = ['ATGGCGTGCA', 'CGTGCAATGG', 'CAATGGCGTG', 'GGCGTGC']

	sol = deNovo_on_DWave_QA(reads)
	print_test_results(1, expected_sol, sol, print_all=True)

run_test_1()

# quboDict = quboFile_to_quboDict("denovo_001.qubo")
# print(quboDict)
# solve_qubo_exact(quboDict, all=True)
# embed_qubo_chimera(quboAdjM)
