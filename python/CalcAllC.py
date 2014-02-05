import FileIO
import numpy
import sys
import math

###
#	Params for Interval Selection
###
MAX_INTERVALS_N2 = 100
MIN_LENGTH_N2 = 1000000 # 1Mb
MAX_INTERVALS_N3 = 20 
MIN_LENGTH_N3 = 5000000 # 5Mb

def select_intervals_n3(lengths, tumor_counts, norm_counts, m, upper_bounds, lower_bounds, copy, tau, force, num_intervals):
	"""
	Selects num intervals from the input file for n=3 analysis
	Returns lines in original input order
	"""
	
	if tau != 2: 
		print "ERROR: For automatic interval selection with 3 subpopulations, the default copy number (--TAU) must be 2. To run with other values, bounds must be provided in the input file."
		exit(1)


	b = int(math.ceil(num_intervals*.75))
	c = int(num_intervals - b)

	indexes = [i for i in range(len(lengths)) if lengths[i] >= MIN_LENGTH_N3]

	lines = [[i, lengths[i], tumor_counts[i], norm_counts[i], upper_bounds[i], lower_bounds[i], copy[i]] for i in indexes]
	lines.sort(key=lambda x: -x[1])

	### Select b longest intervals that have C != 2
	### Select up to c longest intervals that have C == 2 && ub == 2 
	intervals = []
	for i, line  in enumerate(lines):
		if c > 0 and line[6] == 2 and line[4] == 2:
			intervals.append(i)
			c -= 1
		elif b > 0 and line[6] in [0,1,3]:
			intervals.append(i)
			b -= 1

	for i, line  in enumerate(lines):
		if c > 0 and line[6] == 2 and line[4] > 2:
			intervals.append(i)
			c -= 1


	if c > 0 or b > 0:
		if not force:
			print "WARNING: This sample isn't a good candidate for THetA analysis with 3 subpopulations: There aren't a sufficient number of intervals that fit the criteria for interval selection. Run with --FORCE flag to ignore this warning. Exiting..."
			exit(1)
		else:
			print "WARNING: This sample isn't a good candidate for THetA analysis with 3 subpopulations: There aren't a sufficient number of intervals that fit the criteria for interval selection."
	topLines = [lines[i] for i in intervals]
	### Calculate new bounds
	for line in topLines:
		c = line[6]
		# If c = 0 then lb and ub don't change
		if c == 0: 
			pass 
		# If c = 1 then lb = 1 and ub doesn't change
		elif c == 1:
			line[5] = 1
		# If c = 2 then lb = 1 and ub = min(3, ub)
		elif c == 2:
			line[5] = 1
			line[4] = min(3, line[4])
		# If c = 3 then lb = lb and ub = 3
		else:
			line[4] = 3



	topLines.sort(key=lambda x: x[0])

	print "\tSelected", len(intervals), "intervals for analysis."
	return [column(topLines, i) for i in range(len(topLines[0]))]

def select_intervals_n2(lengths, tumor_counts, norm_counts, m, upper_bounds, lower_bounds, k, force, num_intervals):
	"""
	Selects num intervals for n=2 analysis
	Returns lines in original input order
	"""
	indexes = filter_intervals_n2(lengths, tumor_counts, norm_counts, m, k)

	total_length = float(sum(lengths))
	if upper_bounds is None or lower_bounds is None:
		lines = [[i, lengths[i], tumor_counts[i], norm_counts[i]] for i in indexes]
	else:
		lines = [(i, lengths[i], tumor_counts[i], norm_counts[i], upper_bounds[i], lower_bounds[i]) for i in indexes]
	lines.sort(key=lambda x: x[1])

	lim = min(num_intervals, len(indexes))
	topLines = lines[-lim:]

	new_total_length = float(sum([topLines[i][1] for i in range(len(topLines))]))
	if new_total_length < 0.1*total_length:
		if not force:
			print "WARNING: This sample isn't a good candidate for THetA analysis. The longest ", len(indexes), "intervals chosen for analysis represent <10% of the combined length of all provided intervals. Run with --FORCE flag to ignore this warning. Exiting..." 
			exit(1)
		else:
			print "WARNING: This sample isn't a good candidate for THetA analysis. The longest ", len(indexes), "intervals chosen for analysis represent <10% of the combined length of all provided intervals." 
	topLines.sort(key=lambda x: x[0])

	
	print "\tSelected", len(topLines), "intervals for analysis."
	if upper_bounds is None or lower_bounds is None:
		return [column(topLines, i) for i in range(len(topLines[0]))] + [None, None]
	else:
		return [column(topLines, i) for i in range(len(topLines[0]))]

def filter_intervals_n2(lengths, tumor_counts, norm_counts, m,  k):
	"""
	Filters out intervals that aren't long enough or overly amplified
	"""
	total_length = float(sum(lengths))
	total_tumor = float(sum(tumor_counts))
	total_normal = float(sum(norm_counts))
	indexes = range(m)
	indexes = [i for i in indexes if tumor_counts[i] > 0 and norm_counts[i] > 0 and lengths[i] >= MIN_LENGTH_N2]
	indexes = [i for i in indexes if ((tumor_counts[i]/total_tumor) / (norm_counts[i]/total_normal)) < float(k+1)/2]
	return indexes

def column_filter(lines, column, intervals):
	return [lines[i][column] for i in intervals]

def column(lines, column):
	return [line[column] for line in lines]

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
		
