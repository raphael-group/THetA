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


import FileIO
import numpy
import sys
import math

###
#	Params for Interval Selection
###
MIN_LENGTH_N2 = 1000000 # 1Mb
MIN_LENGTH_N3 = 5000000 # 5Mb

def select_intervals_n3(lengths, tumor_counts, norm_counts, m, upper_bounds, lower_bounds, copy, tau, force, num_intervals):
	"""
	Selects num intervals from the input file for n=3 analysis
	Returns lines in original input order
	"""

	if tau != 2: 
		print "ERROR: For automatic interval selection with 3 subpopulations, the default copy number (--TAU) must be 2. To run with other values, bounds must be provided in the input file."
		exit(1)


	# Filter out intervals not used in n=2 analysis (those that don't have bounds set)


	interval_used = [x != "X" for x in upper_bounds]
	real_indexes = [i for i in range(m) if interval_used[i]]
	lengths = [v for i,v in enumerate(lengths) if interval_used[i]]
	tumor_counts = [v for i,v in enumerate(tumor_counts) if interval_used[i]]
	norm_counts = [v for i,v in enumerate(norm_counts) if interval_used[i]]
	upper_bounds = [int(v) for i,v in enumerate(upper_bounds) if interval_used[i]]
	lower_bounds = [int(v) for i,v in enumerate(lower_bounds) if interval_used[i]]
	copy = [int(v) for i,v in enumerate(copy) if interval_used[i]]

	b = int(math.ceil(num_intervals*.75))
	c = int(num_intervals - b)

	indexes = [i for i in range(len(lengths)) if lengths[i] >= MIN_LENGTH_N3]
	
	lines = [[real_indexes[i], lengths[i], tumor_counts[i], norm_counts[i], upper_bounds[i], lower_bounds[i], copy[i]] for i in range(len(real_indexes)) if lengths[i] >= MIN_LENGTH_N3]
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

def select_intervals_n2(lengths, tumor_counts, norm_counts, m, k, force, num_intervals):
	"""
	Selects num intervals for n=2 analysis
	Returns lines in original input order
	"""
	indexes = filter_intervals_n2(lengths, tumor_counts, norm_counts, m, k)

	total_length = float(sum(lengths))
	lines = [[i, lengths[i], tumor_counts[i], norm_counts[i]] for i in indexes]
	lines.sort(key=lambda x: x[1])

	lim = min(num_intervals, len(indexes))
	topLines = lines[-lim:]

	new_total_length = float(sum([topLines[i][1] for i in range(len(topLines))]))
	if new_total_length < 0.1*total_length:
		if not force:
			print "WARNING: This sample isn't a good candidate for THetA analysis. The longest ", lim, "intervals chosen for analysis represent <10% of the combined length of all provided intervals. Run with --FORCE flag to ignore this warning. Exiting..." 
			exit(1)
		else:
			print "WARNING: This sample isn't a good candidate for THetA analysis. The longest ", lim, "intervals chosen for analysis represent <10% of the combined length of all provided intervals." 
	topLines.sort(key=lambda x: x[0])

	
	print "\tSelected", len(topLines), "intervals for analysis."
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

def column(lines, column):
	return [line[column] for line in lines]

		
