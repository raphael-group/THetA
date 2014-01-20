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

from Optimizer import Optimizer
from Enumerator import Enumerator
import time

def time_estimate(n,m,k,tau,lower_bounds, upper_bounds, r, rN, \
		    max_normal, sorted_index, num_processes, multi_event):
	"""
		Estimates the runtime with the specified bounds and parameters
	"""

	print "Estimating Runtime. This may take a few minutes for large numbers of intervals."
	if n is 3 and m > 30:
		print "With n=3 and", m, "intervals, the runtime would likely be greater than several days. Try reducing the number of intervals below 25"
		return

	enum = Enumerator(n, m, k, tau, lower_bounds, upper_bounds, multi_event)
	opt = Optimizer(r, rN, m, n,tau, upper_bound=max_normal)

	if n == 2:
		TEST_NUM=100
		C = enum.generate_next_C()
		if C is False:
			print "Error: No valid Copy Number Profiles exist for these intervals within the bounds specified. Exiting..."
			sys.exit(1)
		start = time.clock()
		for i in range(TEST_NUM):
			try:
				soln = opt.solve(C)
				C = enum.generate_next_C()
			except:	break
		end = time.clock()
		count = count_number_matrices_2(m,upper_bounds, lower_bounds)
	else:
		TEST_NUM=20
		C = enum.generate_next_C()
		count = 1
		start = time.clock()
		dayCount = 0
		while C is not False:
			try:
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
			except: break
		if count == 0:
			print "Error: No valid Copy Number Profiles exist for these intervals within the bounds specified. Exiting..."
			sys.exit(1)
	avgVal = float(end-start)/TEST_NUM
	seconds = count * (float(end-start)/TEST_NUM) / num_processes
	print "Estimated Total Time:",
	if seconds < 60: print int(seconds + .5), "seconds"
	elif seconds < 3600: print int((seconds/60)+.5) , "minutes"
	else: print int((seconds/3600) + .5), "hours"


def count_number_matrices_2(m, upper_bounds, lower_bounds):
	"""
	Calculates the number of possible matrices for n=2
	"""
	upper_bounds = [int(v) for v in upper_bounds]
	lower_bounds = [int(v) for v in lower_bounds]
	possValues = [0]*(max(upper_bounds)+1)
	for i in range(lower_bounds[0], upper_bounds[0]+1): possValues[i] += 1

	for i in range(m-1):
		row1=i
		row2=i+1
		possValuesNew = [0]*(max(upper_bounds)+1)
		for j, v in enumerate(possValues):
			if v > 0:
				minVal = max(j, lower_bounds[row2])
				maxVal = upper_bounds[row2]
				for k in range(minVal, maxVal+1):
					possValuesNew[k] += v
		possValues = possValuesNew
	return sum(possValues)

