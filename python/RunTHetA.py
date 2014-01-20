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

def process_loop(queue, opt, returnQueue, sorted_index):
	"""
	Code that each child process executes. Repeatedly pops of new C
	values from queue until it reaches an exit signal. Then puts its results
	on the return queue and finishes

	Arguments:
	queue (multiprocessing.Queue): Task queue, containing C matrices
	opt (Optimizer): instance of an optimizer
	returnQueue (multiprocessing.Queue): Queue to put results in
	sorted_index (list): Array containing ordering information for sorting
	"""
	min_likelihood = float('inf') 
	best = []
	while True:
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
	

def start_processes(max_processes, queue, opt, returnQueue, sorted_index):
	"""
	Starts a max_processes number of processes, and starts them
	"""
	processes = [Process(target=process_loop, args=(queue, opt, returnQueue,\
			    sorted_index), name=i+1) for i in range(max_processes-1)]
	for p in processes:
		p.daemon = True
		p.start()
	return processes

def find_mins(best):
	"""
	Takes a the list of "best" C,mu pairs returned by each process and finds 
	the ones with the minimum likelihood
	"""
	min_likelihood = float('inf')
	true_best = []
	for solns in best:
		likelihood = solns[0][2]
		if isClose([min_likelihood], [solns[0][2]]):
			true_best += solns
		elif likelihood < min_likelihood:
			min_likelihood = likelihood
			true_best = solns
	return true_best

def do_optimization(n,m,k,tau,lower_bounds, upper_bounds, r, rN, \
		    max_normal, sorted_index, max_processes):
	"""
	Performs the optimization for the given parameters with max_proccesses
	number of processes
	Returns a list of the best C matrices and associated mu values 
	and likelihoods
	"""
	enum = Enumerator(n, m, k, tau, lower_bound=lower_bounds, \
			    upper_bound=upper_bounds)
	opt = Optimizer(r, rN, m, n,tau, upper_bound=max_normal)
	MAX_QUEUE_SIZE = int(10E6)
	queue = Queue(MAX_QUEUE_SIZE) #Task queue for the processes
	returnQueue = Queue(MAX_QUEUE_SIZE) #Shared queue for processes to return results

	processes = start_processes(max_processes, queue, opt, returnQueue, \
			    sorted_index)
	
	C = enum.generate_next_C()
	count = 0
	while C is not False:
		count += 1
		queue.put(C, True)
		C = enum.generate_next_C()
	if count == 0:
		print "Error: No valid Copy Number Profiles exist for these intervals within the bounds specified. Exiting..."
		sys.exit(1)


	# Send STOP signal to all processes
	for i in range(max_processes-1):
		queue.put(0)

	for p in processes:
		p.join()

	best = []
	while not returnQueue.empty():
		item = returnQueue.get()
		best.append(item)
	best = find_mins(best)
	return best

def do_optimization_single(n,m,k,tau,lower_bounds, upper_bounds, r, rN, \
		    max_normal, sorted_index):
	"""
	Performs the optimization for the given parameters with a single process
	Returns a list of the best C matrices and associated mu values 
	and likelihoods
	"""

	enum = Enumerator(n, m, k, tau, lower_bound=lower_bounds, \
			    upper_bound=upper_bounds)
	opt = Optimizer(r, rN, m, n,tau, upper_bound=max_normal)
	min_likelihood = float("inf")	
	best = []
	count = 0
	
	C = enum.generate_next_C()
	while C is not False:
		count += 1
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
		C = enum.generate_next_C()

	if count == 0: 
		print "Error: No valid Copy Number Profiles exist for these intervals within the bounds specified. Exiting..."
		sys.exit(1)
	return best

def time_estimate(n,m,k,tau,lower_bounds, upper_bounds, r, rN, \
		    max_normal, sorted_index, num_processes):
	"""
		Estimates the runtime with the specified bounds and parameters
	"""

	print "Estimating Runtime. This may take a few minutes for large numbers of intervals."
	if n is 3 and m > 30:
		print "With n=3 and", m, "intervals, the runtime would likely be greater than several days. Try reducing the number of intervals below 25"
		return

	TEST_NUM=20
	enum = Enumerator(n, m, k, tau, lower_bound=lower_bounds, \
			    upper_bound=upper_bounds)
	opt = Optimizer(r, rN, m, n,tau, upper_bound=max_normal)

	C = enum.generate_next_C()
	count = 0
	start = time.clock()
	dayCount = 0
	while C is not False:
		C = enum.generate_next_C()
		if count < TEST_NUM:
			soln = opt.solve(C)
		if count == TEST_NUM:
			end = time.clock()
			avgVal = float(end-start)/TEST_NUM
			dayCount = (86400 * num_processes / avgVal) * 1.25
		count += 1
		if dayCount > 0 and count > dayCount:
			print "With the current parameters, the estimated runtime is greater than a day. We suggest reducing the number of intervals, adding bounds or increasing the number of processes before re-running."
			return

	seconds = count * (float(end-start)/TEST_NUM) / num_processes
	print "Estimated Total Time:",
	if seconds < 60: print int(seconds + .5), "seconds"
	elif seconds < 3600: print int((seconds/60)+.5) , "minutes"
	else: print int((seconds/3600) + .5), "hours"

def main():
	###
	#  Read in arguments and data file
	##
	filename, n, k, tau, directory, prefix, max_normal, bound_heuristic, \
		normal_bound_heuristic,heuristic_lb, heuristic_ub, num_processes, \
		bounds_only, estimate_time,multi_event = parse_arguments()
	print "Reading in query file..."
	lengths, tumorCounts, normCounts, m, upper_bounds, lower_bounds = read_interval_file(filename)

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

	write_out_bounds(directory, prefix, filename, upper_bounds, lower_bounds)
	if bounds_only: sys.exit(0)

	if estimate_time: 
		time_estimate(n,m,k,tau,lower_bounds,upper_bounds,r,rN,max_normal,sorted_index, num_processes)
	###
	#  Initialize optimizer and enumerator 
	###
	print "Performing optimization..."

	if num_processes == 1:
		best = do_optimization_single(n,m,k,tau,lower_bounds,upper_bounds,
			r,rN,max_normal,sorted_index)
	else:
		best = do_optimization(n, m, k, tau, lower_bounds, upper_bounds, r, rN,\
			    max_normal, sorted_index, num_processes)

	###
	#  Write results out to file
	###
	upper_bounds = reverse_sort_list(upper_bounds, sorted_index)
	lower_bounds = reverse_sort_list(lower_bounds, sorted_index)
	write_out_result(directory, prefix, best)	

import time
if __name__ == '__main__':
	main()

