from FileIO import *
from DataTools import *
from Misc import *
from Enumerator import Enumerator
from Optimizer import Optimizer

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
	C = enum.generate_next_C()
	min_likelihood = float("inf")	
	best = []
	count = 0
	
	###
	#  Iterate through possible interval count matrices and perform optimization
	###
	while C is not False:
		count += 1
		if count%5000 == 0:
			print "\tIteration #",count
		
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

	###
	#  Write results out to file
	###
	upper_bounds_real = reverse_sort_list(upper_bounds, sorted_index)
	lower_bounds_real = reverse_sort_list(lower_bounds, sorted_index)
	write_out_bounds(directory, prefix, filename, upper_bounds_real, lower_bounds_real)
	write_out_result(directory, prefix, best)	

if __name__ == '__main__':
	  main()
