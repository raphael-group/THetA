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

import numpy

def calculate_bounds_heuristic(x,r,rN,m,tau,k):
	print "Calculating bounds using bound heuristic"
	# Normalize tumor and normal counts
	sum_r = sum(r)
	r_norm = [float(i)/sum_r for i in r]
	sum_rN = sum(rN)
	rN_norm = [float(i)/sum_rN for i in rN]

	ratios = [(t/n) for (t,n) in zip(r_norm,rN_norm)] 
	# Mean and standard deviation of the ratio of tumor to normal counts
	mean = (1.0/m) * sum(ratios)
	std_dev =((1.0/(m-1)) * sum([(mean - ratio)**2 for ratio in ratios]))**(.5)
	c = mean + (x*std_dev)

	lower_bounds = [0]*m
	upper_bounds = [tau]*m
	for i, ratio in enumerate(ratios):
		if ratio > c:
			y = round(tau*ratio)
			lower_bounds[i] = max(tau, y-1)
			upper_bounds[i] = max(k, y+1)
	return upper_bounds, lower_bounds

def calculate_bounds_normal_heuristic(normal_bound_heuristic, heuristic_lb, heuristic_ub, r, rN, m, k):
	print "Calculating bounds using normal bound heuristic"
	# Normalize tumor and normal counts
	sum_r = sum(r)
	r_norm = [float(i)/sum_r for i in r]
	sum_rN = sum(rN)
	rN_norm = [float(i)/sum_rN for i in rN]

	ratios = [(t/n) for (t,n) in zip(r_norm,rN_norm)] 

	upper_bounds = [normal_bound_heuristic]*m
	lower_bounds = [normal_bound_heuristic]*m

	for j, ratio in enumerate(ratios):
		if ratio < heuristic_lb: 
			lower_bounds[j] = 0
			upper_bounds[j] = normal_bound_heuristic
		elif ratio > heuristic_ub:
			if ratio > 2:
				y = round(normal_bound_heuristc * ratio)
				lower_bounds[j] = y-1
				upper_bounds[j] = max(k, y+1)
			else:
				lower_bounds[j] = normal_bound_heuristic
				upper_bounds[j] = k
	return upper_bounds, lower_bounds

def sort_r(rN, r):
	"""
	Sort r in ascending order by the ratio r[i]:rN[i]

	Args:
		rN (list of ints): normal read depth vector
		r (list of ints): tumor read depth vector
	"""
	sum_n = 1.0*sum(rN)
	sum_t = 1.0*sum(r)
	#ratio is the ratio between rN and r after normalization
	ratio = [(t*1.0/n)*(sum_n/sum_t) for (n,t) in zip(rN, r)]
	#ratioNum just appends the index, so we can keep track of it after sorting
	#and be able to reverse the sort at the end
	ratioNum = [(ratio[i],i) for i in range(len(ratio))]
	ratioNum = sorted(ratioNum, key = lambda k: k[0])
	sorted_index = [v[1] for v in ratioNum]
	#sort tumor counts based on the ordering we calculated above
	r = [r[sorted_index[i]] for i in range(len(sorted_index))]
	rN = [rN[sorted_index[i]] for i in range(len(sorted_index))]
	return (r, rN, sorted_index)

def sort_by_sorted_index(vec, sorted_index):
	"""
	Sort vector by the ordering specified in sorted_index.
	That is, an object in location i in the new list, will be the object that
		was in the location sorted_index[i] in vec

	Args:
		vec (list): list with len(list) = len(sorted_index)
		sorted_index (list of ints): desired order for vec
	"""
	return [vec[sorted_index[i]] for i in range(len(sorted_index))]

def reverse_sort_C(C, sorted_index):
	"""
	Perform the reverse of sort described in sort_by_sorted_index, on rows in a numpy array.
	
	Args:
		C (numpy.array): array with C.shape[0] = len(sorted_index)
		sorted_index (list of ints): desired order for rows of C
	"""
	m,n = C.shape
	C_new = numpy.zeros(C.shape)
	for i in range(len(sorted_index)):
		row = sorted_index[i]
		for j in range(n):
			C_new[row][j] = C[i][j]
	return C_new
	
def reverse_sort_list(vec, sorted_index):
	"""
	Perform the reverse of sort described in sort_by_sorted_index
	
	Args:
		vec (list): list with len(list) = len(sorted_index)
		sorted_index (list of ints): desired order for vec
	"""
	new_vec = [0]*len(sorted_index)
	for i in range(len(sorted_index)):
		new_vec[sorted_index[i]] = vec[i]
	return new_vec
