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
 # @author Layla Oesper, Ahmad Mahmoody, Benjamin J. Raphael, and Gryte Satas
 ###

import numpy

class Enumerator:
	def __init__(self, n, m, k, tau, lower_bound = None, upper_bound = None, multi_event = False):
		"""
		Initialize Enumerator

		Args:
			n (int): number of subpopulations
			m (int): number of intervals
			k (int): maximum copy number
			tau (int): number of copies in the normal genome
			upper_bound(list of ints): max copy number for each interval
			lower_bound(list of ints): min copy number for each interval
		"""
		self.m = m
		self.n = n-1 #N-1 since the first column is fixed. We are only generating n-1 columns
		self.k = k
		self.tau = tau
		self.iter = [0]*m
		self.allow_multi_event = True
		
		self.lower_bound, self.upper_bound = self._check_bound_order(lower_bound,upper_bound)

		if n == 2:
	
			if self.lower_bound is None: self.lower_bound = [0]*m
			if self.upper_bound is None: self.upper_bound = [k]*m
			self.iter = self.lower_bound[:]	
		elif n == 3:
			if self.upper_bound is None: self.upper_bound = [k]*m
			if self.lower_bound is None: self.lower_bound = [0]*m

			#TODO: Implement self.low_val
			#	else: self.low_val = min(self.lower_bound)
			self.generator = self._generate_next_C_3()
			self.rows, self.edges = self._create_graph()
			
	def generate_next_C(self):
		"""
		Generate next valid interval count matrix

		Returns:
			C (numpy.array): Interval count matrix
		"""
		if self.n == 1:
			return self._generate_next_C_2()
		else:
			try:
				return next(self.generator)
			except StopIteration:
				return False


	def _check_bound_order(self,lower_bound, upper_bound):
		"""
		Check that the lower and upper bounds are consistent with the ordering,
		and adjust if needed

		Args:
			lower_bound (list of ints): minimum copy number at each interval
			upper_bound (list of ints): maximum copy number at each interval
		Returns:
			lower_bound (list of ints): adjusted lower_bound
			upper_bound (list of ints): adjusted upper_bound
		"""
		# For lower bounds, each item must be at least as large as theone before it
		if lower_bound is not None:
			for i in range(1,len(lower_bound)):
				if lower_bound[i] < lower_bound[i-1]:
					lower_bound[i] = lower_bound[i-1]
		# For upper bounds, each item must be no larger than the one after it
		if upper_bound is not None:
			for i in reversed(range(len(upper_bound)-1)):
				if upper_bound[i] > upper_bound[i+1]:
					upper_bound[i] = upper_bound[i+1]

		return lower_bound, upper_bound
	
	###
	#	n = 2
	###

	def _generate_next_C_2(self):
		"""
		In the n=2 case, generate next valid interval count matrix

		Returns:
			C (numpy.array): Interval count matrix
		"""
		# Using list self.iter, which is a copy of the n=1 column
		# finds an index in the array where the value is less than
		# the one immediately following, and adds one
		for i in range(self.m-1):
			if self.iter[i] < self.iter[i+1] and self.iter[i] < self.upper_bound[i]:
				self.iter[i] += 1
				for j in range(i):
					self.iter[j]= self.lower_bound[j]
				return self._C_to_array()
	
		# If no such value exists (as described above), increment
		# the value at the last index, and reset all others to 0
		if self.iter[self.m-1] < self.upper_bound[self.m-1]:
			self.iter[self.m-1] += 1
			for j in range(self.m-1):
				self.iter[j] = self.lower_bound[j]
			return self._C_to_array()
	
		# If we can't increment the last value further then there are
		# no more valid arrays
		else:
			return False
	
	def _C_to_array(self):
		""" Converts self.iter into a numpy matrix """
		C = numpy.zeros((self.m,self.n+1))
		for i in range(self.m):
			C[i][0] = self.tau
			C[i][1] = self.iter[i]
		return C

	###
	#   n = 3
	###

	def get_graph(self):
		"""
		Used for time estimate for n=3
		"""
		return self.rows, self.edges
	
	def _generate_next_C_3(self):
		C = [0]*self.m
		count = 0
		for i in range(len(self.rows)):
			if not self._in_bounds(self.rows[i], 0): continue
			# The following check prevents duplicate matrices, which just have 
			# columns in a different order. 
			row = self.rows[i]
			switch = False
			if row[0] > row[1]: continue
			elif row[0] == row[1]: switch = True
			else: switch = False
	
			C[0] = i
			for val in self._generate_next_C_3_recurse(C, 0, float("-inf"), float("inf"), switch):
				yield self._to_matrix(C)
	
	def _generate_next_C_3_recurse(self, C, depth, minVal, maxVal, switch):
		if depth == self.m-1:
			yield C
			return
		for child in self.edges[C[depth]]:
			if not self._in_bounds(self.rows[child], depth+1): continue

			row = self.rows[child]
			# The following check prevents duplicate matrices, which just have 
			# columns in a different order. 
			switchNew = False
			if switch:
				if row[0] > row[1]: continue
				elif row[0] == row[1]: switchNew = True
	
			C[depth+1] = child
			newMin, newMax = self._get_mu_bounds(self.rows[C[depth]], row)
			if newMin: 
				newMin = max(minVal, newMin)
				newMax = maxVal
			elif newMax: 
				newMax = min(maxVal, newMax)
				newMin = minVal
			if newMin <= newMax:
				for val in self._generate_next_C_3_recurse(C, depth+1, newMin, newMax, switchNew):
					yield val
	
	def _to_matrix(self, C):
		CNew = numpy.zeros((self.m, self.n+1))
		for i in range(self.m):
			row = self.rows[C[i]]
			CNew[i][0] = self.tau
			for j in range(1, self.n+1):
				CNew[i][j] = row[j-1]
		return CNew
	
	def _get_mu_bounds(self, row1, row2):
		"""
		Returns the minimum/maximum ratio between mu[1] and mu[2]
		"""
		#TODO: Memoize
		xChange = float(row2[0] - row1[0])
		yChange = float(row2[1] - row1[1])
		if xChange == 0 or yChange == 0:
			return float("-inf"), float("inf")
		elif xChange > 0:
			minVal = yChange/(-xChange)
			return minVal, False
		else:
			maxVal = yChange/(-xChange)
			return False, maxVal

	def _in_bounds(self, row, i):
		return all([a >= self.lower_bound[i] and a <= self.upper_bound[i] for a in row])
	

	###
	#	Initialization methods for n=3
	###
	def _add_one(self, row, maxVal):
		"""
		Used for enumerating possible rows
		Recursively adds one
		"""
		if len(row) == 0: return []
		if row[0] < maxVal:
			row[0] += 1
			return row
		else:
			return [0] + self._add_one(row[1:], maxVal)
	
	def _is_valid_edge(self, row1, row2):
		if row1 == row2: return True
		else: return any([i < j for i,j in zip(row1,row2)])

	def _is_valid_row(self, row):
		# Checks for multi-event cases
		return (self.tau - row[0])*(self.tau - row[1]) >= 0

	def _no_multi_event(self, row):
		# Checks for multi-event cases
		return row[0] == row[1] or row[0] == self.tau or row[1] == self.tau
	
	def _create_graph(self):
		"""
		For enumerating matrices for n=3
		Nodes in the graph are rows, and an edge from v to w in the graph, means
		that row w can follow row v in a matrix.
		"""
		start = [0,0]
		row = [0,0]
		rows = [row[:],]
		badCount = 0
		while True:
			row = self._add_one(row, self.k)
			if start == row: break

			if self.allow_multi_event and self._is_valid_row(row):
				# Where is valid row checks that tehre isnt amp/del in same interval
				rows.append(list(row))
			elif self._no_multi_event(row):
				rows.append(list(row))


	
		edges = []
		for i, row in enumerate(rows):
			edges.append([j for j,row2 in enumerate(rows) if self._is_valid_edge(row, row2)])
		
		return rows, edges

