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
from CalcAllC import weighted_C, L2, L3
import sys
import math


def set_total_read_counts(r, rN):
	global sum_r, sum_rN
	sum_r = float(r)
	sum_rN = float(rN)

def calculate_bounds_heuristic(x,r,rN,m,tau,k):
	print "Calculating bounds using bound heuristic..."
	global sum_r, sum_rN
	# Normalize tumor and normal counts
	r_norm = [float(i)/sum_r for i in r]
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
	print "Calculating bounds using normal bound heuristic..."
	global sum_r, sum_rN
	# Normalize tumor and normal counts
	r_norm = [float(i)/sum_r for i in r]
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
				y = round(normal_bound_heuristic * ratio)
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

	global sum_r, sum_rN
	#ratio is the ratio between rN and r after normalization
	ratio = [(t*1.0/n)*(sum_rN/sum_r) for (n,t) in zip(rN, r)]
	## TODO: Why is it necessary to normalize here? The order should be the
	## same regardless

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


def determine_frac_copy_num(rN, r, lengths, dev):
	"""
	Determines the fraction of the genome that potentially
	contains copy number aberrations.  Used to determine if a sample
	is a reasonable sample to run THetA on.

	Args:
		rN (list of ints): normal read depth vector
		r (list of ints): tumor read depth vector
		lengths (list of ints): length of the intervals
		dev (float): the deviation away from 1.0 for a ratio to
					 potentially be a copy number aberrations.

	Returns:
		frac (float): the fraction of the genome that has a ratio
					  outside of dev away from 1.0.
	"""
	sum_r = sum(r)
	sum_rN = sum(rN)
	m=len(r)
	indexes = range(m)

	low = 1.0 - dev
	up = 1.0 + dev

	tot_len = sum(lengths)
	dev_lens = []

	for i in indexes:
		t = r[i]
		n = rN[i]
		if rN[i] == 0: continue
		ratio = (t*1.0/n)*(1.0*sum_rN/sum_r)
		if ratio > up or ratio < low: dev_lens.append(lengths[i])


	frac = float(sum(dev_lens))/float(tot_len)
	return frac

def un_meta_cluster_bounds(bounds, order, intervalMap):
	"""
	This function will expand the given bounds based on the rows in order
	to the full set of rows in the intervalMap.
	"""
	new_bounds = []
	new_order = []
	for i,v in enumerate(order):
		cur_bound = bounds[i]
		rows = intervalMap[v]

		for r in rows:
			new_order.append(r)
			new_bounds.append(cur_bound)

	return new_bounds, new_order


def un_meta_cluster_results_N2(best, meta_order, intervalMap, allTumor, allNormal):

	newBest = []
	rev_meta_cluster = []
	new_order = []
	r = []
	rN = []
	for i, v in enumerate(meta_order):
		rows = intervalMap[v]
		num_orig_clustered = len(rows)

		rev_meta_cluster = rev_meta_cluster + num_orig_clustered*[i]
		new_order = new_order + rows

	#Reset r and rN
	new_m = len(rev_meta_cluster)
	for x in range(new_m):
		r.append(allTumor[new_order[x]])
		rN.append(allNormal[new_order[x]])

	for c, mu, NLL, p in best:
		m,n = c.shape
		#new_m = len(rev_meta_cluster)

		c_new = numpy.zeros((new_m,n))
		for x in range(new_m):
			for y in range(n):
				c_new[x][y] = c[rev_meta_cluster[x]][y]
			#r.append(allTumor[new_order[x]])
			#rN.append(allNormal[new_order[x]])

		c_weight = weighted_C(c_new, rN)


		likelihood, vals = L2(mu[0], c_weight, len(r), r)

		newBest.append((c_new,mu,likelihood,vals))
		
	return newBest, r, rN


def un_meta_cluster_results_N3(best, meta_order, intervalMap, allTumor, allNormal, n):


	newBest = []
	rev_meta_cluster = []
	new_order = []
	r = []
	rN = []
	for i, v in enumerate(meta_order):
		rows = intervalMap[v]
		num_orig_clustered = len(rows)

		rev_meta_cluster = rev_meta_cluster + num_orig_clustered*[i]
		new_order = new_order + rows

	#Reset r and rN
	new_m = len(rev_meta_cluster)
	for x in range(new_m):
		r.append(allTumor[new_order[x]])
		rN.append(allNormal[new_order[x]])

	for c, mu, NLL, p in best:
		m,n = c.shape
		#new_m = len(rev_meta_cluster)

		c_new = numpy.zeros((new_m,n))
		for x in range(new_m):
			for y in range(n):
				c_new[x][y] = c[rev_meta_cluster[x]][y]
			#r.append(allTumor[new_order[x]])
			#rN.append(allNormal[new_order[x]])

		c_weight = weighted_C(c_new, rN)


		likelihood, vals = L3(mu, c_weight, len(r), r, n)

		newBest.append((c_new,mu,likelihood,vals))
		
	return newBest, r, rN

def score_clusters(intervalMap, lengths, rd, baf, m):
	"""
	This function will score the given cluster assignments in the intervalMap based on the average
	distance to the cluster center.
	"""

	cluster_scores=[float('inf') for x in range(m)]
	for key in intervalMap.keys():

		if key == -1:
			#cluster_scores[key] = float('inf')
			continue

		rows = intervalMap[key]

		cluster_lens = lengths[rows]
		cluster_rd = rd[rows]
		cluster_baf = baf[rows]
		tot_len = sum(cluster_lens)

		#Give small clusters high distance scores
		if tot_len < 1000000: 
			cluster_scores[key] = float('inf')
			continue

		weighted_rd = [p*q for p,q in zip(cluster_lens, cluster_rd)]
		weighted_baf = [p*q for p,q in zip(cluster_lens, cluster_baf)]

		rd_mean = sum(weighted_rd)/float(tot_len)
		baf_mean = sum(weighted_baf)/float(tot_len)

		dists=[math.sqrt((rd_mean -x)**2 + (baf_mean-y)**2) for x,y in zip(cluster_rd, cluster_baf)]
		score = sum(p*q for p,q in zip(cluster_lens,dists))/float(tot_len)

		cluster_scores[key] = score

	return cluster_scores




