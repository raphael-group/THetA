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

from Misc import *

import numpy
from scipy import optimize
import warnings
class Optimizer:
	
	def __init__(self, r, rN, m, n, tau, lower_bound=0, upper_bound=1):
		"""
		Initialize Optimizer

		Args:
			r (list of ints): tumor read depth vector
			rN (list of ints): normal read depth vector
			m (int): number of intervals
			n (int): number of subpopulations
			lower_bound (float): min fraction for normal
			upper_bound (float): max fraction for normal
		"""
		self.r = r
		self.m = m
		self.n = n
		self.rN = rN
		self.lB = lower_bound
		self.uB = upper_bound

		if self.n > 2:
			# Small optimization for dLambda_dMu
			global dLambda_dMu_numers
			dLambda_dMu_numers = [[]]*(self.n)


	def solve(self, C):
		"""
		Run optimization for matrix C

		Args:
			C (numpy.array): Possible interval count matrix
		Returns:
			mu (n-tuple of floats): Optimum value for mu
			likelihood (float):	The likelihood at the optimum mu
			vals (list of floats):
		"""
		
		# Some of the optimization functions will raise runtime warnings on things
		# that do not affect the correctness of the code. This filters them out so
		# they don't clutter up output
		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			if self.n == 2:	
				return self._solve_n2(C)
			elif self.n > 2:
				return self._solve_n3plus(C)

	def _solve_n2(self, C):
		"""
		For the n=2 case, find the optimum value for mu, given an interval 
			count matrix C

		Args:
			C (numpy.array): Possible interval count matrix
		Returns:
			mu (n-tuple of floats): Optimum value for mu
			likelihood (float):	The likelihood at the optimum mu
			vals (list of floats):

		"""

		# Find a root of dL_dMu between 0 and 1

		C_w = weighted_C(C,self.rN)
		C_hat = normalize_C(C_w, self.m, self.n)
		
		lB = self.lB
		uB = self.uB
		if lB != 0: lB = M2_Rev(C_w,self.lB, self.m, self.n)
		if uB != 1: uB = M2_Rev(C_w,self.uB, self.m, self.n)

		global dLdMu_numers
		dLdMu_numers = []

		try:
			val = optimize.brenth(dL_dMu, lB, uB, args = (C_hat,self.m,self.r))
		except Exception:
			#Case that there isn't a root in the interval
			return None

		mu = M2(C_w,val,self.m, self.n)
		likelihood,vals = L2(mu,C_w,self.m, self.r)
		return ((mu, 1-mu), likelihood, vals)

	def _solve_n3plus(self, C):
		"""
		For the n=3 case, find the optimum value for mu, given an interval
			count matrix C

		Args:
			C (numpy.array): Possible interval count matrix
		Returns:
			mu (n-tuple of floats): Optimum value for mu
			likelihood (float):	The likelihood at the optimum mu
			vals (list of floats):
		"""
	
		global dLambda_dMu_numers
		dLambda_dMu_numers = [dLambda_dMu_numers[0]] + [[]]*(self.n)

		# Find a root for derivative functions
		C_w = weighted_C(C,self.rN)
		C_hat = normalize_C(C_w,self.m,self.n)
		start = [1.0/self.n]*(self.n) + [1]
		val = optimize.fsolve(equations, start, args = (self.r,self.m,C_hat, self.n),\
				    fprime = jacobian)
		mu = val[:self.n]
		if not inRange(mu):
			#In the case that we find the wrong root (one of the values is negative),
			# directly minimize the function
			start = [1.0/self.n] * (self.n-1)
			mu = optimize.fmin_bfgs(L3_hat, start, fprime = dL3_hat, args = \
					    (C_hat, self.r, self.m, self.n), disp=0)
			mu = mu.tolist()
			mu.append(1-sum(mu))
			if not inRange(mu): 
				#Case that a minimum doesn't exist
				return None
		answer = M3(C_w,mu,self.m,self.n)
		likelihood, vals = L3(answer,C_w,self.r, self.m, self.n)
		
		return (answer, likelihood, vals)

def normalize_C(C, m, n):
	#Sum up columns
	sum_C = [sum([C[i][j] for i in range(m)]) for j in range(n)]
	C_new = numpy.zeros((m,n))
	for i in range(m):
		for j in range(n):
			C_new[i][j] = C[i][j]/sum_C[j]
	return C_new

def weighted_C(C, rN):
	m,n = C.shape
	C_new = numpy.zeros((m,n))
	for row in range(m):
		for col in range(n):
			C_new[row][col] = rN[row] * C[row][col]
	return C_new

###
# Equations for the n=2 case
###
def L2(mu, C, m, r):
	vals = []
	total_sum = 0
	mu1 = 1-mu
	denom = sum([C[j][0]*mu + C[j][1]*mu1 for j in range(m)])
	for i in range(m):
		numer = C[i][0]*mu + C[i][1]*mu1
		total_sum += r[i] * numpy.log(numer/denom)
		vals.append(numer/denom)
	return (-total_sum, vals)

def L_hat(mu, C_hat, m, r):
	total_sum = 0
	mu1 = 1-mu
	for i in range(m):
		term1 = (C_hat[i][0]*mu)
		term2 = (C_hat[i][1]*mu1)
		total_sum += r[i] * numpy.log(term1+term2)
	total_sum += (1 - sum(mu))
	return -total_sum	

dLdMu_numers = []
def dL_dMu(mu, C_hat, m, r):
	
	# The values in the numerators are going to be the same for every call
	# for a given C, so we can calculate these once and then reuse them
	global dLdMu_numers
	if len(dLdMu_numers) == 0:
		dLdMu_numers = [r[i] * (C_hat[i][0] - C_hat[i][1]) for i in range(m)]

	total_sum = 0
	mu1 = 1-mu
	for i in range(m):
		total_sum += dLdMu_numers[i]/((C_hat[i][0] * mu) + (C_hat[i][1] * mu1))
	return -total_sum

def M2(C, mu, m, n):
	numer = -mu * sum([C[i][1] for i in range(m)])
	denom = (mu - 1)* sum([C[i][0] for i in range(m)]) + numer
	return numer/denom
	
def M2_Rev(C, mu, m, n):
	numer = -mu * sum([C[i][0] for i in range(m)])
	denom = (mu - 1)* sum([C[i][1] for i in range(m)]) + numer
	return numer/denom

###
# Equations for the n>=3 case
###
def L3(mu, C, r, m, n):
	total_sum = 0
	vals = []
	for i in range(m):
		numer = sum([C[i][j]*mu[j] for j in range(n)])
		denom = sum([C[h][j]*mu[j] for j in range(n) for h in range(m)])
		total_sum += r[i] * numpy.log(numer/denom)
		vals.append(numer/denom)
	return (-total_sum, vals)

def L3_hat(mu, C_hat, r, m, n):
	total_sum = 0
	munew = mu.tolist()
	munew.append(1-sum(mu))
	for i in range(m):
		total_sum += r[i] * numpy.log(sum([C_hat[i][j] * munew[j] for j in range(n)]))
	return -total_sum


def dL3_hat(mu, C_hat, r, m, n):
	# Gradient of L3_hat (for n=3 only) 
	vals = numpy.zeros((2))
	for i in range(m):
		numer0 = (C_hat[i][0] - C_hat[i][2])
		numer1 = (C_hat[i][1] - C_hat[i][2])
		denom = (C_hat[i][0] - C_hat[i][n-1])*mu[0] + \
			(C_hat[i][1] - C_hat[i][n-1])*mu[1] + C_hat[i][2] 
		vals[0] += r[i] * (numer0/denom)
		vals[1] += r[i] * (numer1/denom)
	return vals

def dLambda_dL(x, args):
	
	r, m, n, C_hat = args
	mu = x[:n]
	return 1 - sum(mu)

def dLambda_dMu(x, args, k):
	r, m, n, C_hat = args
	mu = x[:n]
	L = x[n]

	global dLambda_dMu_numers	
	if len(dLambda_dMu_numers[k]) == 0:
		dLambda_dMu_numers[k] = [r[i] * C_hat[i][k] for i in range(m)]
	
	total_sum = 0
	for i in range(m):
		total_sum += dLambda_dMu_numers[k][i]/ \
			sum([C_hat[i][j] * mu[j] for j in range(n)])
	return (-total_sum) - L

def jacobian(x, r, m, C_hat, n):
	# jacobian for system of equations defined in equations
	mu = x[:n]
	L = x[n]
	jac = numpy.zeros((n+1,n+1))
	for i in range(n+1):
		jac[n][i] = -1
		jac[i][n] = -1
	jac[n][n] = 0

	for i in range(n):
		for j in range(n):
			jac[i][j] = second_deriv(x,r, m, C_hat,n,i,j)
	return jac

def second_deriv(x, r, m, C_hat, n, k, h):
	mu = x[:n]
	L = x[n]
	total_sum = 0
	for i in range(m):
		numer = r[i] * C_hat[i][k] * C_hat[i][h]
		denom = sum([C_hat[i][j] * mu[j] for j in range(n)])**2
		total_sum += numer/denom
	return total_sum

def equations(x, r, m, C_hat, n):
	args = (r,m,n,C_hat)
	eqs = [dLambda_dMu(x,args,i) for i in range(n)] + [dLambda_dL(x,args)] 
	return eqs

def M_eq(mu_new, C, mu, m, n):
	csums = [sum([C[i][h] for i in range(m)]) for h in range(n)]
	eqs = [0]*(n+1)
	for j in range(n):
		temp = sum([mu_new[h] * csums[h] for h in range(n)])
		eqs[j] = (mu[j] * temp) - (mu_new[j] * csums[j]) - mu_new[n]
	eqs[n] = sum(mu_new[:n]) - 1
	return eqs

def M3(C, mu, m, n):
	start = [.33]*n + [0]
	val = optimize.fsolve(M_eq, start, args = (C,mu,m,n))
	return val[:n]
