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
from Optimizer import *

def L2(mu, C, m, r):
	vals = []
	total_sum = 0
	mu1 = 1-mu
	valid_rows = [i for i in range(m) if C[i][0] != 0]
	if mu == 0.0:
		valid_rows = [i for i in range(m) if C[i][1] != 0]
	denom = sum([C[j][0]*mu + C[j][1]*mu1 for j in valid_rows])
	for i in range(m):
		if i in valid_rows:
			numer = C[i][0]*mu + C[i][1]*mu1
			total_sum += r[i] * numpy.log(numer/denom)
			vals.append(numer/denom)
		else:
			vals.append("X")
	return (-total_sum, vals)



def L3(mu, C, m, r, n):
	total_sum = 0
	vals = []
	valid_rows = [i for i in range(m) if C[i][0] != 0]

	denom = sum([C[h][j]*mu[j] for j in range(n) for h in valid_rows])
	for i in range(m):
		if i in valid_rows:
			numer = sum([C[i][j]*mu[j] for j in range(n)])
			total_sum += r[i] * numpy.log(numer/denom)
			vals.append(numer/denom)
		else:
			vals.append("X")
	return (-total_sum, vals)

def calculateX(tumorI, normalI, sumR, sumAll, mu, n, row, h):
	# Calculate the (real) value for C[m][h] that minimizes the likelihood
	# by calculating the value that makes d/dx(L)=0

	# nR is the normalized read count for interval i.
	row = [r*normalI for r in row]
	nR = float(tumorI)/(sumR + tumorI)

	# val0 = sum(x_g,k * mu_k)
	sumRow = sum([row[i] * mu[i] for i in range(n) if i != h]) 

	return float(nR*(sumAll+sumRow) - sumRow)/((1-nR)*mu[h])


def calc_all_c_2(best, r, rN, all_tumor, all_normal, intervals_used):
	bestNew = []
	num_intervals = len(all_tumor)
	for c,mu,likelihood, vals in best:
		m,n = c.shape
		c_new = numpy.zeros((m+1,n))
	
		for x in range(m):
			for y in range(n):
				c_new[x][y] = c[x][y]

		c_new = weighted_C(c_new, rN + [0,])
		c_all = numpy.zeros((len(all_tumor), n))
		for i, val in enumerate(intervals_used):
			c_all[val] = c[i]

		sum_all = sum([c_new[j][k]*mu[k] for j in range(m) for k in range(n)]) 
		sum_r = sum(r)

		for i in range(num_intervals):
			if i not in intervals_used:
				if all_normal[i] == 0:
					c_all[i][0] = 2
					c_all[i][1] = -1
					continue

				c_all[i][0] = 2

				wX = calculateX(all_tumor[i], all_normal[i], sum_r, sum_all, mu, n, [2,0], 1)
				x = wX/all_normal[i]
				if x < 0: 
					c_all[i][1] = 0 
					continue

				bot = math.floor(x)
				top = math.ceil(x)
			
				c_new[m][0] = 2*all_normal[i]
				c_new[m][1] = bot*all_normal[i]
				lBot = L2(mu[0], c_new, m+1, r+[all_tumor[i],])
				c_new[m][1] = top*all_normal[i]
				lTop = L2(mu[0], c_new, m+1, r+[all_tumor[i],])
			
				if lBot[0] < lTop[0]:
					c_all[i][1] = int(bot)
				else:
					c_all[i][1] = int(top)

		c_all_w = weighted_C(c_all, all_normal)
		likelihood, vals = L2(mu[0], c_all_w, len(all_tumor), all_tumor)
		bestNew.append((c_all, mu, likelihood, vals))
	return bestNew

def calc_all_c_3(best, r, rN, all_tumor, all_normal, intervals_used):
	#Calculate optima where you set either one to two
	#Then iterate up for both the ssame
	bestNew = []
	num_intervals = len(all_tumor)

	for c,mu,likelihood, vals in best:
		m,n = c.shape

		# C_new is a copy of the orignal C, with an extra row at the end.
		c_new = numpy.zeros((m+1,n))
		for x in range(m):
			for y in range(n):
				c_new[x][y] = c[x][y]
		c_new = weighted_C(c_new, rN + [0,])

		# C_all is the copy numbers profile including all of the intervals, 
		# which we're calculating below.
		c_all = numpy.zeros((len(all_tumor), n))
		for i, val in enumerate(intervals_used):
			c_all[val] = c[i]

		# |C_new*mu|, excluding the unknown row.
		sum_all = sum([c_new[j][k]*mu[k] for j in range(m) for k in range(n)]) 
		sum_r = sum(r)

		for i in range(num_intervals):
			if i not in intervals_used:
				c_all[i][0] = 2
				if all_normal[i] == 0:
					# If the number of reads aligning to the interval in the
					# normal sample is 0, there's isn't an optimal value.
					c_all[i][0] = 2
					c_all[i][1] = -1
					c_all[i][2] = -1
					continue

				candidates = []
				c_new[m][0] = 2*all_normal[i]
				c_new[m][2] = 2*all_normal[i]

				#FIX X, Calculate min/max

				wX = calculateX(all_tumor[i], all_normal[i], sum_r, sum_all, mu, n, [2,0,2], 1)
				x = wX/all_normal[i]
				top = int(max(0,math.ceil(x)))
				bot = int(max(0,math.floor(x)))

				c_new[m][1] = bot*all_normal[i]
				xBot = L3(mu, c_new, m+1, r+[all_tumor[i],],n)[0]


				c_new[m][1] = top*all_normal[i]
				xTop = L3(mu, c_new, m+1, r+[all_tumor[i],],n)[0]

				candidates.append((xBot,[bot,2]))
				candidates.append((xTop,[top,2]))

				#FIX Y, Calculate min/max

				wY = calculateX(all_tumor[i], all_normal[i], sum_r, sum_all, mu, n, [2,2,0], 2)
				y = wY/all_normal[i]
				top = int(max(0,math.ceil(y)))
				bot = int(max(0,math.floor(y)))

				c_new[m][1] = 2*all_normal[i]
				
				c_new[m][2] = bot*all_normal[i]
				yBot = L3(mu, c_new, m+1, r+[all_tumor[i],],n)[0]

				c_new[m][2] = top*all_normal[i]
				yTop = L3(mu, c_new, m+1, r+[all_tumor[i],],n)[0]

				candidates.append((yBot,[2,bot]))
				candidates.append((yTop,[2,top]))


				#FIX X==Y, Calculate all until we find min
				prev = float('inf')
				j = 0
				while True:
					c_new[m][1] = j*all_normal[i]
					c_new[m][2] = j*all_normal[i]
					l = L3(mu, c_new, m+1, r+[all_tumor[i],],n)[0]
					candidates.append((l, [j,j]))
					j += 1
					if l > prev: break
					prev = l

				candidates.sort()
				rowMin = candidates[0][1]
				c_all[i][1] = rowMin[0]
				c_all[i][2] = rowMin[1]

		c_all_w = weighted_C(c_all, all_normal)

		likelihood, vals = L3(mu, c_all_w, len(all_tumor), all_tumor, n)
		bestNew.append((c_all, mu, likelihood, vals))
	return bestNew

def calc_all_c_3_multi_event(best, r, rN, all_tumor, all_normal, intervals_used):
	bestNew = []
	num_intervals = len(all_tumor)

	for c,mu,likelihood, vals in best:
		m,n = c.shape

		# C_new is a copy of the orignal C, with an extra row at the end.
		c_new = numpy.zeros((m+1,n))
		for x in range(m):
			for y in range(n):
				c_new[x][y] = c[x][y]
		c_new = weighted_C(c_new, rN + [0,])

		# C_all is the copy numbers profile including all of the intervals, 
		# which we're calculating below.
		c_all = numpy.zeros((len(all_tumor), n))
		for i, val in enumerate(intervals_used):
			c_all[val] = c[i]

		# |C_new*mu|, excluding the unknown row.
		sum_all = sum([c_new[j][k]*mu[k] for j in range(m) for k in range(n)]) 
		sum_r = sum(r)

		for i in range(num_intervals):
			if i not in intervals_used:
				c_all[i][0] = 2
				if all_normal[i] == 0:
					# If the number of reads aligning to the interval in the
					# normal sample is 0, there's isn't an optimal value.
					c_all[i][1] = -1
					c_all[i][2] = -1
					continue

				# Weighted optimal x value (float)
				wX = calculateX(all_tumor[i], all_normal[i], sum_r, sum_all, mu, n, [2,0,0], 1)
				maxX = math.ceil(wX/all_normal[i])

				c_new[m][0] = 2*all_normal[i]

				lMin = float('inf')
				rowMin = None

				if maxX < 0: maxX = 0

				for x in range(int(maxX)+1):
					wX = x * all_normal[i]
					c_new[m][1] = wX

					wY = calculateX(all_tumor[i], all_normal[i], sum_r, sum_all, mu, n, [2,x,0], 2)
					y = wY/all_normal[i]

					bot = int(max(0,math.floor(y)))
					top = int(max(0,math.ceil(y)))

					if x < 2:
						bot = min(bot,2)
						top = min(top,2)
					elif x > 2:
						bot = max(2,bot)
						top = max(2,top)

					c_new[m][2] = bot*all_normal[i]
					lBot = L3(mu, c_new, m+1, r+[all_tumor[i],],n)
					if lBot[0] < lMin:
						lMin = lBot[0]
						rowMin = [2,x,bot]

					c_new[m][2] = top*all_normal[i]
					lTop = L3(mu, c_new, m+1, r+[all_tumor[i],], n)

					if lTop[0] < lMin: 
						lMin = lTop[0]
						rowMin = [2,x,top]


				c_all[i][1] = rowMin[1]
				c_all[i][2] = rowMin[2]

		c_all_w = weighted_C(c_all, all_normal)

		likelihood, vals = L3(mu, c_all_w, len(all_tumor), all_tumor, n)
		bestNew.append((c_all, mu, likelihood, vals))
	return bestNew

