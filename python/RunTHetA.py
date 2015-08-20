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

from SelectIntervals import *
from FileIO import *
from DataTools import *
from Misc import *
from Enumerator import Enumerator
from Optimizer import Optimizer
from TimeEstimate import *
from CalcAllC import *
from plotResults import *
from RunBAFModel import run_BAF_model


from multiprocessing import JoinableQueue, Queue, Process, Array, current_process
import os

def process_loop(queue, opt, returnQueue, sorted_index, get_values):
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
	if get_values: solns = []
	while True:
		C = queue.get()
		if C is 0:
			returnQueue.put(best)
			break
		soln = opt.solve(C)
		if soln is not None:
			(mu, likelihood,vals) = soln
			if get_values: solns.append((C,mu,likelihood,vals))
			if isClose([likelihood],[min_likelihood]):
				C_new = reverse_sort_C(C,sorted_index)
				vals = reverse_sort_list(vals, sorted_index)
				best.append((C_new, mu, likelihood, vals))
			elif likelihood < min_likelihood:
				C_new = reverse_sort_C(C,sorted_index)
				vals = reverse_sort_list(vals, sorted_index)
				best = [(C_new, mu, likelihood, vals)]
				min_likelihood = likelihood

	if get_values:
		with open(pre+"."+"values" + str(current_process().name),'w') as f:
			for C,mu,likelihood,vals in solns:
				m,n = C.shape
				stringC = "".join((str(int(C[i][1])) for i in range(m)))
				valsStr = " ".join((str(v) for v in vals))
				f.write(stringC+"\t"+str(mu[0])+"\t"+str(likelihood)+"\t"+valsStr+"\n")
	

def start_processes(max_processes, queue, opt, returnQueue, sorted_index, get_values):
	"""
	Starts a max_processes number of processes, and starts them
	"""
	processes = [Process(target=process_loop, args=(queue, opt, returnQueue,\
			    sorted_index, get_values), name=i+1) for i in range(max_processes-1)]
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
		if len(solns) == 0: continue
		likelihood = solns[0][2]
		if isClose([min_likelihood], [solns[0][2]]):
			true_best += solns
		elif likelihood < min_likelihood:
			min_likelihood = likelihood
			true_best = solns
	return true_best

def do_optimization(n,m,k,tau,lower_bounds, upper_bounds, r, rN, \
		    max_normal, sorted_index, max_processes, multi_event, get_values):
	"""
	Performs the optimization for the given parameters with max_proccesses
	number of processes
	Returns a list of the best C matrices and associated mu values 
	and likelihoods
	"""
	enum = Enumerator(n, m, k, tau, lower_bounds, upper_bounds, multi_event)
	opt = Optimizer(r, rN, m, n,tau, upper_bound=max_normal)
	MAX_QUEUE_SIZE = int(10E6)
	try:
		queue = Queue(MAX_QUEUE_SIZE) #Task queue for the processes
	except OSError:
		MAX_QUEUE_SIZE = 2**15-1
		queue = Queue(MAX_QUEUE_SIZE)

	returnQueue = Queue(MAX_QUEUE_SIZE) #Shared queue for processes to return results

	processes = start_processes(max_processes, queue, opt, returnQueue, \
			    sorted_index, get_values)
	
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

	best = []
	for i in range(len(processes)):
		item = returnQueue.get()
		best.append(item)

	for p in processes:
		p.join()

	best = find_mins(best)
	return best

def do_optimization_single(n,m,k,tau,lower_bounds, upper_bounds, r, rN, \
		    max_normal, sorted_index, multi_event, get_values):
	"""
	Performs the optimization for the given parameters with a single process
	Returns a list of the best C matrices and associated mu values 
	and likelihoods
	"""

	enum = Enumerator(n, m, k, tau, lower_bounds, upper_bounds, multi_event)
	opt = Optimizer(r, rN, m, n,tau, upper_bound=max_normal)
	min_likelihood = float("inf")	
	best = []
	count = 0
	C = enum.generate_next_C()
	if get_values: solns = []
	while C is not False:
		count += 1
		soln = opt.solve(C)
		if soln is not None:
			(mu, likelihood,vals) = soln

			if get_values: solns.append((C,mu,likelihood))
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

	if get_values:
		with open(pre+"."+"likelihoods",'w') as f:
			for C,mu,likelihood in solns:
				m,n = C.shape
				stringC = "".join((str(int(C[i][1])) for i in range(m)))
				f.write(stringC+"\t"+str(mu[0])+"\t"+str(likelihood)+"\n")
	
	if count == 0: 
		print "Error: No valid Copy Number Profiles exist for these intervals within the bounds specified. Exiting..."
		sys.exit(1)
	return best

def best_near_max_contamination(best, max_normal):
	for C, mu, likelihood, vals in best:
		if abs(max_normal-mu[0]) < .01: return True
	return False


def main():
	###
	#  Read in arguments and data file
	##
	args = parse_arguments()
	intervals = read_interval_file(args[0])
	
	ub, lb = run_fixed_N(2, args, intervals)
	intervals[4] = ub
	intervals[5] = lb
	ub, lb = run_fixed_N(3, args, intervals)

def run_fixed_N(n, args, intervals):
	(filename, results, N, k, tau, directory, prefix, max_normal, bound_heuristic, \
		normal_bound_heuristic,heuristic_lb, heuristic_ub, num_processes, \
		bounds_only, multi_event, force, get_values, choose_intervals, num_intervals, \
		read_depth_file, graph_format, runBAF, tumorSNP, normalSNP, ratio_dev, min_frac,\
		cluster_bounds, density_bounds) = args

	global pre
	pre = prefix

	print "Reading in query file..."
	
	lengths, tumorCounts, normCounts, m, upper_bounds, lower_bounds = intervals

	###
	#	Determine is sample has enough copy number aberrations to run
	###
	frac = determine_frac_copy_num(normCounts, tumorCounts, lengths, ratio_dev)
	print "Frac with potential copy numbers:", frac
	if frac < min_frac:
		print "ERROR: This sample does not have enough large copy number aberrations to be a good candidate for tumor composition estimation using THetA.  See --RATIO_DEVIATION and --MIN_FRAC flags to modify how the potential presence of large copy number aberrations is determined.  Exiting..."
		exit(1)
	###
	#	Automatically Select Intervals
	#	note: This is the default behavior
	###

	if cluster_bounds is not None:
		#Add code here
		print "Setting bounds using clustering needs to be implemented."
	elif density_bounds is not None:
		# Add code here
		print "Setting bounds using density needs to be implemented."



	if choose_intervals:
		print "Selecting intervals..."
		allM, allLengths, allTumor, allNormal, allUpperBounds, allLowerBounds = (m, lengths, tumorCounts, normCounts, upper_bounds, lower_bounds)

		if n == 2:
			order, lengths, tumorCounts, normCounts = select_intervals_n2(lengths, tumorCounts, normCounts, m, k, force, num_intervals)
			
			# 07-09-15 - Comment out b/c we want to be able to set specific bounds
			#upper_bounds = None  
			#lower_bounds = None
		elif n == 3:
			if results is None: 
				print "ERROR: No results file supplied. Unable to automatically select intervals for n=3 without results of n=2 analysis. See --RESULTS flag, or --NO_INTERVAL_SELECTION to disable interval selection. Exiting..."
				exit(1)
			else: 
			# Need to read in original file, bounds file and results file. Original file needed because copy numbers are based on 
				copy = read_results_file(results)
				order, lengths, tumorCounts, normCounts, upper_bounds, lower_bounds, copy = select_intervals_n3(lengths, tumorCounts, normCounts, m, upper_bounds, lower_bounds, copy, tau, force, num_intervals)

		m = len(order)

	sum_r = sum(tumorCounts)
	sum_rN = sum(normCounts)
	set_total_read_counts(sum_r, sum_rN)
	###
	#  Process/sort read depth vectors and calculate bounds if necessary
	###
	print "Preprocessing data..."

	r,rN,sorted_index = sort_r(normCounts,tumorCounts)


	if normal_bound_heuristic is not False:
		upper_bounds,lower_bounds = calculate_bounds_normal_heuristic( \
			normal_bound_heuristic, heuristic_lb, heuristic_ub, r, rN, m, k)
	elif bound_heuristic is not False or upper_bounds is None and lower_bounds is None:
		if bound_heuristic is False: bound_heuristic = 0.5
		upper_bounds,lower_bounds = calculate_bounds_heuristic(float(bound_heuristic),\
			 r, rN, m, tau, k)
	else: 
		if upper_bounds is not None: upper_bounds = sort_by_sorted_index(upper_bounds,\
			sorted_index)
		if lower_bounds is not None: lower_bounds = sort_by_sorted_index(lower_bounds,\
			sorted_index)
		


	###Bounds files in their original orders
	ub_out = reverse_sort_list(upper_bounds, sorted_index)
	lb_out = reverse_sort_list(lower_bounds, sorted_index)
	if choose_intervals:
		write_out_bounds(directory, prefix, filename, ub_out, lb_out, n, order)
	else: 
		write_out_bounds(directory, prefix, filename, ub_out, lb_out, n)

	if bounds_only: sys.exit(0)



	enum = time_estimate(n,m,k,tau,lower_bounds,upper_bounds,r,rN,max_normal,sorted_index, num_processes, multi_event, force)
	###
	#  Initialize optimizer and enumerator 
	###
	print "Performing optimization..."

	if num_processes == 1:
		best = do_optimization_single(n,m,k,tau,lower_bounds,upper_bounds,
			r,rN,max_normal,sorted_index, multi_event, get_values)
	else:
		best = do_optimization(n, m, k, tau, lower_bounds, upper_bounds, r, rN,\
			    max_normal, sorted_index, num_processes, multi_event, get_values)
	if best == []:
		print "ERROR: Maximum Likelihood Solution not found within given bounds."
		exit(1)

	if n == 2 and best_near_max_contamination(best, max_normal):
		print "WARNING: At least one of the top solutions is near the upper bound on normal contamination. Further analysis may required as the sample likely falls into one of the following categories:\n\t1. This sample has high normal contamination. Consider re-running with an increased normal contamination upper bound. See --MAX_NORMAL option\n\t2. This sample may not satisfy the assumption that most of the tumor genome retains the normal expected copynumber (e.g. a genome duplication event has occurred). See THetA optional parameters in changing the expected copy number.\n\t3. This sample may not be a good candidate for THetA analysis (i.e. does not contain large copy number aberrations that distinguish populations)."
	r = reverse_sort_list(r, sorted_index)
	rN = reverse_sort_list(rN, sorted_index)
	if choose_intervals:
		if n == 2:
			best = calc_all_c_2(best, r, rN, allTumor, allNormal, order)
		elif n == 3 and not multi_event: 
			best = calc_all_c_3(best, r, rN, allTumor, allNormal, order)
		else:
			best = calc_all_c_3_multi_event(best, r, rN, allTumor, allNormal, order)

		#lko 7-6-15 Fix multiple solutions to only return the solution with the overal min NLL
		best = find_mins(best)	

	#run BAF model on results to determine most likely solution
	if runBAF and tumorSNP is not None and normalSNP is not None:
		if len(best) != 1:
			BAFprefix = prefix + ".preliminary"
			write_out_result(directory, BAFprefix, best, n)

			resultsFile = BAFprefix + ".n"+str(n)+".results"
			resultsPath = os.path.join(directory, resultsFile)
			try:
				run_BAF_model(tumorSNP, normalSNP, filename, resultsPath, prefix=prefix + ".n" + str(n), directory=directory)
			except IOError:
				print "ERROR: Invalid locations for tumor and normal SNP files. The BAF model will not be run. You can try running the BAF model again directly from the runBAFModel.py script."
		else:
			write_out_result(directory, prefix, best, n)
	elif runBAF and (tumorSNP is None or normalSNP is None):
		print "ERROR: Need file location for tumor and normal SNP files to run the BAF model. The BAF model will not be run. You can try running the BAF model again directly from the runBAFModel.py script."
		write_out_result(directory, prefix, best, n)
	else:
		write_out_result(directory, prefix, best, n)

	###
	# Make Results Plots
	###
	print "Plotting results as a " + graph_format + "..."
	#plot_results(directory, filename, prefix, n)
	plot_results(directory, filename, prefix, read_depth_file, n, graph_format)

	if n == 2: write_out_N3_script(directory, prefix, filename)

	return ub_out, lb_out

import time
if __name__ == '__main__':
	main()

