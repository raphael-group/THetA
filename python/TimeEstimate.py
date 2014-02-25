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
		    max_normal, sorted_index, num_processes, multi_event, force):
	"""
		Estimates the runtime with the specified bounds and parameters
	"""


	print "Estimating time..."
	if n is 3 and m > 30 and not force:
		print "\tWARNING: With n=3 and", m, "intervals, the runtime would likely be excessive. Try reducing the number of intervals below 25. Run with --FORCE to continue."
		exit(1)

	enum = Enumerator(n, m, k, tau, lower_bounds, upper_bounds, multi_event)
	opt = Optimizer(r, rN, m, n,tau, upper_bound=max_normal)

	if n == 2:
		TEST_NUM=100
		count = count_number_matrices_2(m,upper_bounds, lower_bounds)
	else:
		TEST_NUM=20
		count = count_number_matrices_3(m,upper_bounds, lower_bounds, enum)
	C = enum.generate_next_C()
	if C is False:
		print "ERROR: No valid Copy Number Profiles exist for these intervals within the bounds specified. Exiting..."
		sys.exit(1)
	start = time.clock()
	for i in range(TEST_NUM):
		try:
			soln = opt.solve(C)
			C = enum.generate_next_C()
		except:	break

	end = time.clock()
	avgVal = float(end-start)/TEST_NUM
	seconds = count * (float(end-start)/TEST_NUM) / num_processes
	print "\tEstimated Total Time:",
	if seconds < 60: print int(seconds + .5), "second(s)"
	elif seconds < 3600: print int((seconds/60)+.5) , "minute(s)"
	else: 
		hours = int((seconds/3600) + .5)

		print int((seconds/3600) + .5), "hour(s)"
		if hours > 200: 
			if not force:
				print "WARNING: With the current settings, the runtime is likely excessive. To reduce runtime, try:\n\t1) Increase the number of processes used with the --NUM_PROCESSES flag.\n\t2) Reduce the number of intervals chosen using the --NUM_INTERVALS flag.\n\t3) Disable automatic interval selection using --NO_INTERVAL_SELECTION, and hand-select a smaller number of intervals, or set tighter bounds on the current intervals.\n\t Run with --FORCE to continue with current settings."
				exit(1)
	return enum

		


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

def count_number_matrices_3(m, upper_bounds, lower_bounds, enum):
	"""
	Calculates the number of possible matrices for n=3
	Unlike for n=2, this isn't an exact value, and is probably 
	an overestimate
	It accounts for matrices with columns switched by just 
	dividing the final value in half (when the real value is 
	slightly over half), and it doesn't take into account at all
	matrices which would be filtered out because there aren't valid mu values.
	"""

	rows, edges = enum.get_graph()
	rowBounds = [(min(row), max(row)) for row in rows]

	upper_bounds = [int(v) for v in upper_bounds]
	lower_bounds = [int(v) for v in lower_bounds]
	possValues = [0]*(len(rows))
	for i in range(len(rows)):
		if rowBounds[i][0] >= lower_bounds[0] and rowBounds[i][1] <= upper_bounds[0]:
			possValues[i] += 1
	for i in range(m-1):
		possValuesNew = [0]*(len(rows))
		for j, v in enumerate(possValues):
			if v > 0:
				for k in edges[j]:
					if all([a >= lower_bounds[i+1] and a <= upper_bounds[i+1] for a in rows[k]]):
						possValuesNew[k] += v
		possValues = possValuesNew

	return sum(possValues)/2

