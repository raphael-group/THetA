import FileIO
import numpy
import sys

# Take an input file
# 1 method: Takes input, sorts all inputs by lengths, returns top 100 values expected 

def get_top_intervals_by_length(lengths, tumor_counts, norm_counts, m, upper_bounds, lower_bounds, num, k):
	indexes = filter_intervals_n2(lengths, tumor_counts, norm_counts, m, k)

	total_length = float(sum(lengths))
	if upper_bounds is None or lower_bounds is None:
		lines = [[i, lengths[i], tumor_counts[i], norm_counts[i]] for i in indexes]
	else:
		lines = [(i, lengths[i], tumor_counts[i], norm_counts[i], upper_bounds[i], lower_bounds[i]) for i in indexes]
	lines.sort(key=lambda x: x[1])

	lim = min(num, len(indexes))
	topLines = lines[-lim:]
	print "LINES:", lim

	new_total_length = float(sum([topLines[i][1] for i in range(len(topLines))]))
	if new_total_length < 0.1*total_length:
		print "WARNING: This sample isn't a good candidate for THetA analysis. The longest ", len(indexes), "intervals chosen for analysis represent <10% of the combined length of all provided intervals." 

	if upper_bounds is None or lower_bounds is None:
		return [column(topLines, i) for i in range(len(topLines[0]))] + [None, None]
	else:
		return [column(topLines, i) for i in range(len(topLines[0]))]

def filter_intervals_n2(lengths, tumor_counts, norm_counts, m,  k):
	MIN_LENGTH = 1000000
	total_length = float(sum(lengths))
	total_tumor = float(sum(tumor_counts))
	total_normal = float(sum(norm_counts))
	
	indexes = range(m)
	indexes = [i for i in indexes if tumor_counts[i] > 0 and norm_counts[i] > 0 and lengths[i] >= MIN_LENGTH]
	indexes = [i for i in indexes if ((tumor_counts[i]/total_tumor) / (norm_counts[i]/total_normal)) < float(k+1)/2]
	
	return indexes

def column(array, column):
	return [m[column] for m in array]

def L2New(mu, C, m, r, row, rRow):
	total_sum = 0
	mu1 = 1-mu
	denomC = sum([C[j][0]*mu + C[j][1]*mu1 for j in range(m)])
	denom = denomC + row[0]*mu + row[1]*mu1
	for i in range(m):
		numer = C[i][0]*mu + C[i][1]*mu1
		total_sum += r[i] * numpy.log(numer/denom)
	numer = row[0]*mu + row[1]*mu1
	total_sum += rRow * numpy.log(numer/denom)
	return -total_sum

def normalize_C(C, row, m, n):
	#Sum up columns
	sum_C = [sum([C[i][j] for i in range(m)]) for j in range(n)]
	sum_C[0] += row[0]
	sum_C[1] += row[1]

	C_new = numpy.zeros((m,n))
	for i in range(m):
		for j in range(n):
			C_new[i][j] = C[i][j]/sum_C[j]
	row_new = (row[0]/sum_C[0], row[1]/sum_C[1])
	return C_new, row_new

def weighted_C(C, rN):
	m,n = C.shape
	C_new = numpy.zeros((m,n))
	for row in range(m):
		for col in range(n):
			C_new[row][col] = rN[row] * C[row][col]
	return C_new
		

def calc_all_c(best, r,rN, topNum, m, order, all_tumor, all_normal, lower_bounds, upper_bounds, tau):

	print r
	ratios = [float(normal)/tumor for (tumor, normal) in zip(all_tumor, all_normal)]
	newBest = []

	for C,mu,likelihood,vals in best:
		C_w = weighted_C(C, rN)
		Cnew = numpy.zeros((m,2))

		for i in range(m): Cnew[i][0] = tau
		for i, val in enumerate(order):
			Cnew[val][1] = C[i][1]
		for i in range(m):
			if Cnew[i][1] != 0: continue
			minVal = float("inf")
			minNum = None
			for j in range(lower_bounds[i], upper_bounds[i]+1):
				row = (tau, j)
				row_w = (all_normal[i] * tau, all_normal[i] * j)
				#CNorm, rownorm = normalize_C(C_w, row_w, topNum, 2)
				#print row, ratios[i]
				val = L2New(mu[0], C_w, topNum, r, row_w, all_tumor[i])
				if val < minVal: 
					minVal = val
					minNum = j
				print j, val
			Cnew[i][1] = minNum
			print minNum
		newBest.append((Cnew,mu,likelihood,vals))
	return newBest
		
