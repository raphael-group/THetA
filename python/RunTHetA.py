 ###
 # 2013 Brown University, Providence, RI.
 #
 #                       All Rights Reserved
 #
 # Permission to use, copy, modify, and distribute this software and its
 # documentation for any purpose other than its incorporation into a
 # commercial product is hereby granted without fee, provided that the
 # above copyright notice appear in all copies and that both that
 # copyright notice and this permission notice appear in supporting
 # documentation, and that the name of Brown University not be used in
 # advertising or publicity pertaining to distribution of the software
 # without specific, written prior permission.
 #
 # BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
 # INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
 # PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
 # ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 # WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 # ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 # OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 # http://cs.brown.edu/people/braphael/software.html
 # 
 # @author Layla Oesper, Ahmad Mahmoody, Benjamin J. Raphael and Gryte Satas
 ###

from FileIO import *
from DataTools import *
from Misc import *
from Enumerator import Enumerator
from Optimizer import Optimizer

from multiprocessing import JoinableQueue, Queue, Process, Array

def task(queue, opt, returnQueue, sorted_index):
	ID = 1
	min_likelihood = float('inf') 
	best = []
	count = 0
	while True:
		count += 1
		C = queue.get()
		if C is 0:
			returnQueue.put(best)
			break

		soln = opt.solve(C)
		if soln is not None:
			(mu, likelihood,vals) = soln
					
			if isClose([likelihood],[min_likelihood]):
				C_new = reverse_sort_C(C,sorted_index)
				vals = reverse_sort_list(vals, sorted_index)
				best.append((C_new, mu, likelihood, vals))
			elif likelihood < min_likelihood:
				C_new = reverse_sort_C(C,sorted_index)
				vals = reverse_sort_list(vals, sorted_index)
				best = [(C_new, mu, likelihood, vals)]
				min_likelihood = likelihood
	

def start_threads(max_processes, queue, opt, returnQueue, sorted_index):
	processes = [Process(target=task, args=(queue, opt, returnQueue, sorted_index), name=i+1) for i in range(max_processes-1)]
	for p in processes:
		p.daemon = True
		p.start()
	return processes

def find_mins(best):
	min_likelihood = float('inf')
	trueBest = []
	for solns in best:
		likelihood = solns[0][2]
		if isClose([min_likelihood], [solns[0][2]]):
			trueBest += solns
		elif likelihood < min_likelihood:
			min_likelihood = likelihood
			trueBest = solns
	return trueBest

def main():
	###
	#  Read in arguments and data file
	##
	filename, n, k, tau, directory, prefix, max_normal, bound_heuristic, \
		normal_bound_heuristic,heuristic_lb, heuristic_ub = parse_arguments()
	print "Reading in query file..."
	tumorCounts, normCounts, m, upper_bounds, lower_bounds = read_interval_file(filename)


	###
	#  Process/sort read depth vectors and calculate bounds if necessary
	###
	print "Preprocessing data..."
	r,rN,sorted_index = sort_r(normCounts,tumorCounts)
	if bound_heuristic is not False or upper_bounds is None and lower_bounds is None:
		if bound_heuristic is False: bound_heuristic = 0.5

		upper_bounds,lower_bounds = calculate_bounds_heuristic(float(bound_heuristic),\
			 r, rN, m, tau, k)
	elif normal_bound_heuristic is not False:
		upper_bounds,lower_bounds = calculate_bounds_normal_heuristic( \
			normal_bound_heuristic, heuristic_lb, heuristic_ub, r, rN, m, k)
	else: 
		if upper_bounds is not None: upper_bounds = sort_by_sorted_index(upper_bounds,\
			sorted_index)
		if lower_bounds is not None: lower_bounds = sort_by_sorted_index(lower_bounds,\
			sorted_index)

	###
	#  Initialize optimizer and enumerator 
	###
	print "Performing optimization..."
	enum = Enumerator(n, m, k, tau, lower_bound=lower_bounds, upper_bound=upper_bounds)
	opt = Optimizer(r, rN, m, n,tau, upper_bound=max_normal)

	queue = Queue()
	returnQueue = Queue()

	#best = [[]]*max_threads
	
	processes = start_threads(max_threads, queue, opt, returnQueue, sorted_index)
	###
	#  Iterate through possible interval count matrices and perform optimization
	###
	
	C = enum.generate_next_C()
	while C is not False:
		queue.put(C)
		C = enum.generate_next_C()
	for i in range(max_threads-1):
		queue.put(0)

	for p in processes:
		  p.join()

	best = []
	while not returnQueue.empty():
		item = returnQueue.get()
		best.append(item)
	best = find_mins(best)

	###
	#  Write results out to file
	###
	upper_bounds_real = reverse_sort_list(upper_bounds, sorted_index)
	lower_bounds_real = reverse_sort_list(lower_bounds, sorted_index)
	write_out_bounds(directory, prefix, filename, upper_bounds_real, lower_bounds_real)
	write_out_result(directory, prefix, best)	

#import cProfile
if __name__ == '__main__':
	main()

	#global max_threads
	#j = 2
	#for i in range(10):
	#	print "-------------------------", j, "---------------------------"
	#	max_threads=j
	#	j += 2
	#	cProfile.run('main()')
