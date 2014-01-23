import FileIO
import numpy
import sys

# Take an input file
# 1 method: Takes input, sorts all inputs by lengths, returns top 100 values expected 

def get_top_intervals_by_length(lengths, tumor_counts, norm_counts, m, upper_bounds, lower_bounds, num):


	if upper_bounds is None or lower_bounds is None:
		lines = [[i, lengths[i], tumor_counts[i], norm_counts[i]] for i in range(m)]
	else:
		lines = [(i, lengths[i], tumor_counts[i], norm_counts[i], upper_bounds[i], lower_bounds[i]) for i in range(m)]
	lines.sort(key=lambda x: x[1])


	lim = min(num, m)
	topLines = lines[-lim:]

	if upper_bounds is None or lower_bounds is None:
		return [column(topLines, i) for i in range(len(topLines[0]))] + [None, None]
	else:
		return [column(topLines, i) for i in range(len(topLines[0]))]

def column(array, column):
	return [m[column] for m in array]

##############################################

#filename = sys.argv[1]
#lengths, tumor_counts, norm_counts, m, upper_bounds, lower_bounds = FileIO.read_interval_file(filename)
#upper_bounds = [3]*m
#lower_bounds = [0]*m
#get_top_100_intervals_by_length(lengths, tumor_counts, norm_counts, m, upper_bounds, lower_bounds)
#exit(0)
#
def L2New(mu, C, m, r, row, rRow):
	total_sum = 0
	mu1 = 1-mu
	denomC = sum([C[j][0]*mu + C[j][1]*mu1 for j in range(m)])
	denom = denomC + row[0]*mu + row[1]*mu1
	for i in range(m):
		numer = C[i][0]*mu + C[i][1]*mu1
		total_sum += r[i] * numpy.log(numer/denom)
	total_sum += rRow * numpy.log((row[0]*mu + row[1]*mu1)/denom)
	return -total_sum

#
#
#
#
#### Given a matrix C for some subset of intervals and a fixed mu, as well as rn
#### for the subset, and all intervals, calculate the max likelihood copy number 
#### for the rest of the intervals.
## Suppose you have C for the top n intervals by length:
#array = [[1, 0],[1, 1],[1, 1], [1, 2],[1, 1]]
#C = numpy.array(array)
#
#length = 5
#mu = .5
#filename = sys.argv[1]
#lengths, tumor_counts, norm_counts, m, upper_bounds, lower_bounds = FileIO.read_interval_file(filename)
#
#print lengths
#print tumor_counts
#print norm_counts
#print upper_bounds
#print lower_bounds
#
## Sort all by length
#tumor_counts = [x for (y,x) in sorted(zip(lengths, tumor_counts), key=lambda pair: pair[0])]
#norm_counts = [x for (y,x) in sorted(zip(lengths, norm_counts), key=lambda pair: pair[0])]
#if upper_bounds: upper_bounds = [x for (y,x) in sorted(zip(lengths, upper_bounds), key=lambda pair: pair[0])]
#if lower_bounds: lower_bounds = [x for (y,x) in sorted(zip(lengths, lower_bounds), key=lambda pair: pair[0])]
#lengths = sorted(lengths)
#ratios = [float(t)/n for (t,n) in zip(tumor_counts, norm_counts)]
#
#print lengths
#print tumor_counts
#print norm_counts
#print upper_bounds
#print lower_bounds
#print ratios
#

def calc_all_c(best, r, topNum, m, order, all_tumor, all_normal, lower_bounds, upper_bounds, tau):
	C, mu, likelihood, vals = best

	ratios = [float(tumor)/normal for (tumor, normal) in zip(all_tumor, all_normal)]

	copyNums = [0 for _ in range(m)]
	for i, val in enumerate(order):
		copyNums[val] = C[i][1]

	for i in range(m):
		if copyNums[i] is not 0: continue
		minVal = float("inf")
		minNum = None
		for j in range(lower_bounds[i], upper_bounds[i]):
			row = (tau, j)
			val = L2New(mu[0], C, topNum, r, row, ratios[i])
			if val < minVal: 
				minVal = val
				minNum = j
		copyNums[i] = minNum
