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
	def __init__(self, n, m, k, tau, upper_bound = None, lower_bound = None):
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
		
		self.lower_bound, self.upper_bound = self._check_bound_order(lower_bound,upper_bound)

		if n == 2:
	
			if self.lower_bound is None: self.lower_bound = [0]*m
			if self.upper_bound is None: self.upper_bound = [k]*m
			self.iter = self.lower_bound[:]	
		elif n == 3:
			
			if self.upper_bound is None: self.upper_bound = [k]*m
			if self.lower_bound is None: self.lower_bound = [0]*m
			else: self.low_val = min(self.lower_bound)
			self._gen_nvectors()
			
			# Upper bounds for n=3 are the index of the last row vector that 
			# should be used for that row, so the start index of the largest value
			self.upper_bound = [self.start_indexes[self.upper_bound[i]+1] - 1 \
				if self.upper_bound[i] < k else (len(self.nVectors) - 1) \
				for i in range(m)]
			
	def generate_next_C(self):
		"""
		Generate next valid interval count matrix

		Returns:
			C (numpy.array): Interval count matrix
		"""
		if self.n == 1:
			return self._generate_next_C_2()
		else:
			return self._generate_next_C_3()


	def _check_bound_order(self,lower_bounds, upper_bounds):
		"""
		Check that the lower and upper bounds are consistent with the ordering,
		and adjust if needed

		Args:
			lower_bounds (list of ints): minimum copy number at each interval
			upper_bounds (list of ints): maximum copy number at each interval
		Returns:
			lower_bounds (list of ints): adjusted lower_bounds
			upper_bounds (list of ints): adjusted upper_bounds
		"""
		# For lower bounds, each item must be at least as large as theone before it
		if lower_bounds is not None:
			for i in range(1,len(lower_bounds)):
				if lower_bounds[i] < lower_bounds[i-1]:
					lower_bounds[i] = lower_bounds[i-1]
		# For upper bounds, each item must be no larger than the one after it
		if upper_bounds is not None:
			for i in reversed(range(len(upper_bounds)-1)):
				if upper_bounds[i] > upper_bounds[i+1]:
					upper_bounds[i] = upper_bounds[i+1]

		return lower_bounds, upper_bounds
	
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
	def _generate_next_C_3(self):
		"""
		In the n=3 case, generate next valid interval count matrix

		Returns:
			C (numpy.array): Interval count matrix
		"""
	
		# Get the next valid combination of vectors
		valid = False
		while not valid and self.iter is not False:
			self._add_one_with_reset(self.m-1)
			# Check if every value in the row is at least the lower bound for that row
			if self.iter is not False and True not in \
				[min(self.nVectors[self.iter[i]]) < self.lower_bound[i] for i in range(self.m)]:
				valid = True

		# Generate a matrix from self.iter	
		if self.iter is not False:
			C = numpy.zeros((self.m,self.n+1))
			for i in range(self.m):
				C[i][0] = self.tau
				vec = self.nVectors[self.iter[i]]
				for j in range(1,self.n+1):
					C[i][j] = vec[j-1]
			return C
		else:
			return False
   
	def _add_one_with_reset(self, index):
		"""
		Recursively add one to self.iter and carry down if necessary
		"""

		###
		# The simple version of the code below which may be easier to understand
		# looks like:
		#
		# if self.iter[index] < self.upper_bound[index]:
		#	self.iter[index] += 1
		# elif index > 0:
		# 	self._add_one_with_reset(index-1)
		#	if self.iter is not False: self._reset_high_indexes(index)
		# else: 
		#	self.iter = False
		#
		# Everything in addition to that is to skip past "mirrored" matrices,
		# that is, pairs of matrices in which columns 1 and 2 are switched,
		# since they're mathematically equivalent. (This only works for n=3)
		#
		# The two boolean arrays used below:
		# self.mirror: True if the row vector at that index is palindromic
		#		else False
		# self.unique: False if there is a row vector that comes earlier 
		#		in the list nVectors that is that vectors reverse.
		#		else True
		###

		hitMax = False
		
		if self.iter[index] < self.upper_bound[index] and index != 0: 
			# This index hasn't hit the max yet

			mirrors = [self.mirror[i] for i in self.iter[:index]]
			if all(mirrors):
				# If all of the rows that come before this one are
				# mirrored, then this one must be unique
				self.iter[index] += 1
				while (not self.unique[self.iter[index]]):
					if (self.iter[index] < self.upper_bound[index]):
						self.iter[index] += 1
					else:
						hitMax = True; break

				if hitMax: 
					self._add_one_with_reset(index-1)
					if self.iter is not False: self._reset_high_indexes(index)

			else:
				self.iter[index] += 1
		elif self.iter[index] < self.upper_bound[index] and index == 0: 
			# We're on the first row and haven't hit the max yet
			self.iter[0] += 1

			# The first row MUST be unique
			while (not self.unique[self.iter[0]]):
				if self.iter[0] < self.upper_bound[0]:
					self.iter[0] += 1
				else: 
					hitMax = True; break
			if hitMax:
				self.iter = False
		elif index > 0:
			self._add_one_with_reset(index-1)
			if self.iter is not False: self._reset_high_indexes(index)
		else:
			# None of the rows can be incremented further
			self.iter = False

	def _reset_high_indexes(self, index):
		"""
		Reset all the values with indexes >= index, while maintaining ordering 
		constraints
		"""
		valueMin = 0
		for i in range(index):
			tempMin = min(self.nVectors[self.iter[i]])
			valueMin = tempMin if tempMin >= valueMin else max(self.nVectors[self.iter[i]])

		for i in range(index, len(self.iter)):
			self.iter[i] = self.start_indexes[valueMin]
    
	###
	#	Initialization methods for n=3
	###

	def _gen_nvectors(self):
		"""
		Generate the list of all n-vectors (vectors in list form), sorted by the maximum
		value they contain
		"""
   	 
		# Generate all possible n-vectors
		vector = [self.low_val]*self.n
		set_vectors = [vector]
		vector = self._add_with_carry(0,vector[:])
		while vector is not False:
			set_vectors.append(vector)
			vector = self._add_with_carry(0,vector[:])
		
		# Sorts the vectors by the maximum value they contain
		# Basic radix sort
		counts = [0]*(self.k + 1)
		self.nVectors = [0]*len(set_vectors)
		for item in set_vectors:
			max_item = max(item)
			counts[max_item] += 1

		self.start_indexes = [0]*(self.k + 1)
		for i in range(1,self.k+1):
			self.start_indexes[i] = counts[i-1] + self.start_indexes[i-1]
		counts = self.start_indexes[:]
				
		for item in set_vectors:
			max_item = max(item)
			self.nVectors[counts[max_item]] = item
			counts[max_item] += 1
		self.unique = mark_unique(self.nVectors) 
		self.mirror = mark_mirror(self.nVectors)
	def _add_with_carry(self,i,vector):
		"""
		Recursively add one to the vector and carry up if needed
		"""
		if vector[i] < self.k:
			vector[i] += 1
		elif i+1 < self.n:
			vector[i] = self.low_val
		
			vector = self._add_with_carry(i+1,vector)
		else:
			return False
		return vector

def mark_unique(vectors):
	"""
	Creates a boolean array corresponding to the list vectors
	which is false for an vector v at index i, if there there is
	another vector at a position earlier than i that is the reverse
	of v
	"""
	mask = [True]*len(vectors)
	for i,v1 in enumerate(vectors):
		for v2 in vectors[:i]:
			if v1[0] == v2[1] and v1[1] == v2[0]:
				mask[i] = False
	return mask

def mark_mirror(vectors):
	"""
	Creates a boolean array corresponding to the list vectors which is 
	true if the row vector is palindromic and false otherwise
	"""
	mask = [False]*len(vectors)
	for i,v1 in enumerate(vectors):
		if v1[0] == v1[1]: mask[i] = True
	return mask
