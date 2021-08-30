import numpy as np
import math
import matplotlib.pyplot as plt
import sys
import time
from collections import defaultdict

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
				result += str(n)

	# Shift the string so the 0 is on the first position
	double_result = result + result
	i = 0
	while True:
		if double_result[i] == '0':
			return double_result[i:i+num_reads]
		i += 1


# Only for 4-reads cycles, gets a dictionary cycle!
def get_cycle_type(cycle):
	result = transform_cycle_into_str(cycle)
	if result == '':
		return 'Invalid'
	elif result == '0123':
		return 'Type A'
	elif result == '0321':
		return 'Type B'
	elif result == '0213':
		return 'Type C'
	elif result == '0312':
		return 'Type D'
	elif result == '0132':
		return 'Type E'
	elif result == '0231':
		return 'Type F'
	else:
		sys.exit('----- Error in transform cycle! -----')

# ------------------------------------ PLOTTING --------------------------------

def read_response():
    s = open('results.txt', 'r').read()
    return eval(s)

def build_cycle(key):
	cycle = {}
	for elem in list(key):
		cycle[elem[0]] = elem[1]
	return cycle

# Returns a dict {(Type X, energy): occurrences, ...}
# If multiple invalids appear, multiple pairs for the obtained energy will appear
def get_types_to_occurences():
	response = read_response()
	types_to_occurrences = defaultdict(int)
	min_invalid_energy = 10000

	for key, value in response.items():
		cycle = build_cycle(key)
		cycle_type = get_cycle_type(cycle)
		types_to_occurrences[cycle_type] += value[0]
		if cycle_type == 'Invalid':
			min_invalid_energy = min(min_invalid_energy, value[1])

	return types_to_occurrences, min_invalid_energy

energies_dict = {
	'Type A': -7.9811,
	'Type B': -6.927,
	'Type C': -7.4541,
	'Type D': -7.4541,
	'Type E': -7.4541,
	'Type F': -7.4541,
}

def plot_occurrences_bars(show=False, filter_invalid=False):
	types_to_occurrences, min_invalid_energy = get_types_to_occurences()
	energies_dict['Invalid'] = '>{}'.format(min_invalid_energy)

	print('Cycle type\t\t| Ocurrences\t\t| Energy')
	for cycle_type, ocurrences in types_to_occurrences.items():
		print('{}\t\t\t| {}\t\t\t| {}'.format(cycle_type, ocurrences, energies_dict[cycle_type]))

	if filter_invalid:
		del types_to_occurrences['Invalid']

	fig, ax = plt.subplots()

	pos = np.arange(len(types_to_occurrences))

	colors = []
	for i in range(len(types_to_occurrences)):
		colors.append(plt.cm.tab10(i))

	ax.barh(pos, [x for x in types_to_occurrences.values()],
					tick_label=['{}\ne={}'.format(x, energies_dict[x]) 
									for x in types_to_occurrences.keys()],
					align='center', height=0.5, color=colors)

	#ax.set_title('Ocurrences of each type in a 10.000 samples experiment using Quantum Annealing')
	ax.set_xlabel('Ocurrences')


	plt.savefig('../../thesis/figures/experiments/experiment3.png', bbox_inches='tight')
	if show:
		plt.show()



plot_occurrences_bars(show=True, filter_invalid=True)

