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

import string
import os
import argparse

N_VALS = [2,3]
MAX_K = 7

def parse_arguments():
	"""
	Parse command line arguments 
	
	Returns:
		query file: full path to the location of the input file
		n: number of subpopulations
		k: maximum value of k to be considered
		tau: expected copy number for normal genome
		directory: target directory for output files
		prefix: prefix for output files
		max_normal: maximum fraction to consider for normal (only enforced for n=2)
		bound_heuristic: float parameter for bound heuristic if supplied. False by default
		normal_bound_heuristic: int parameter for normal bound heuristic if supplied,
			False by default
		heuristic_lb: lower bound for normal bound heuristic
		heuristic_ub: upper bound for normal bound heuristic
	"""

	parser = argparse.ArgumentParser()
	parser.add_argument("QUERY_FILE", help="Interval file", metavar="QUERY_FILE")
	parser.add_argument("-n","--N", help="Number of subpopulations", metavar="N", \
			type=int, default=2, required=False)
	parser.add_argument("-k","--MAX_K", help="The maximum value allowed for \
			entries in C", metavar="K", default=3, type=int, required=False)
	parser.add_argument("-t","--TAU", help="Expected number of copies in normal \
			genome", default=2, metavar="TAU", type=int, required=False)
	parser.add_argument("-d","--DIR", help="Directory where result file is written\
			to", default="./", metavar="DIR", required=False)
	parser.add_argument("-p","--OUTPUT_PREFIX", help="Prefix for output files created. By\
			default, it will be the beginning of the input file name (i.e.if input\
			filename were example.input, output filed would be example.output and\
			example.withbounds)", default=None, metavar="PRE", required=False)
	parser.add_argument("-m","--MAX_NORMAL", help="The maximum fraction to consider\
			for normal. Only enforced for n=2", default=.5, type=float,\
			metavar="MAX_NORMAL", required=False)
	parser.add_argument("--BOUND_HEURISTIC", metavar="BH", default=False, required=False)
	parser.add_argument("--NORMAL_BOUND_HEURISTIC", metavar="NBH", type=int,\
			default=False,required=False)
	parser.add_argument("--HEURISTIC_LB", metavar="LB", type=float, default=0.9, \
			required=False)
	parser.add_argument("--HEURISTIC_UB", metavar="UB", type=float,default=1.1, \
			required=False)

	args = parser.parse_args()

	filename = args.QUERY_FILE
	
	n = args.N

	if n not in N_VALS:
		err_msg = "Invalid value entered for n: "+str(n)+". Currently supported values for n: "+str(N_VALS)
		raise ValueError(err_msg)
	
	k = args.MAX_K
	
	if k not in range(MAX_K):
		err_msg = "Invalid value entered for k: "+str(k)+". Supported values for k: 0-"+str(MAX_K)
		raise ValueError(err_msg)

	tau = args.TAU
	
	if tau < 0:
		err_msg = "Invalid value for tau: "+str(tau)+". Tau must be non-negative"
		raise ValueError(err_msg)

	directory = args.DIR
	prefix = args.OUTPUT_PREFIX
	if prefix == None: prefix = os.path.basename(filename).split(".")[0]
	max_normal = args.MAX_NORMAL
	if max_normal < 0 or max_normal > 1:
		err_msg = "Invalid value for max_normal: "+str(max_normal)+". Max_normal must be between 0 and 1"
		raise ValueError(err_msg)
	
	bound_heuristic = args.BOUND_HEURISTIC
	normal_bound_heuristic = args.NORMAL_BOUND_HEURISTIC
	heuristic_lb = args.HEURISTIC_LB
	heuristic_ub = args.HEURISTIC_UB

	
	print "================================================="
	print "Arguments are:"
	print "\tQuery File:", filename
	print "\tn:", n
	print "\tk:", k
	print "\ttau:", tau
	print "\tOutput Directory:", directory
	print "\tOutput Prefix:", prefix
	print "\tMax Normal:", max_normal
	if bound_heuristic is not False:
		print "\tBound Heuristic:", bound_heuristic
	if normal_bound_heuristic is not False:
		print "\tNormal Bound Heuristic:", normal_bound_heuristic
		print "\tHeuristic Lower Bound:", heuristic_lb
		print "\tHeuristic Upper Bound:", heuristic_ub
	print "================================================="
	
	return filename,n,k,tau,directory,prefix,max_normal,bound_heuristic, \
			normal_bound_heuristic, heuristic_lb, heuristic_ub


def read_interval_file(filename):
	"""
	Read in input file

	Args:
		filename (string): full path to the input file 
	Returns:
		tumor_counts: tumor read depth vector
		norm_counts: normal read depth vector
		m: number of intervals
		upper_bounds: list of upper_bounds if supplied, None otherwise
		lower_bounds: list of lower_bounds if supplied, None otherwise
	"""
	  
	f = open(filename)
	lines = f.readlines()
	f.close()

	if "#" in lines[0]: lines = lines[1:]
	m = 0
	tumor_counts = []
	norm_counts = []
	upper_bounds = []
	lower_bounds = []
	for line in lines:
		line = line.split()	
		
		if len(line) < 6 or len(line) > 8:
			raise Exception("Invalid input file format")

		tumor_counts.append(int(line[4]))
		norm_counts.append(int(line[5]))
		if len(line) > 6:
			upper_boundsSupplied = True
			upper_bounds.append(int(line[6]))

		if len(line) > 7:
			lower_boundsSupplied = True
			lower_bounds.append(int(line[7]))
		m += 1

	if len(upper_bounds) == 0: upper_bounds = None	
	if len(lower_bounds) == 0: lower_bounds = None	
	return (tumor_counts, norm_counts, m, upper_bounds, lower_bounds)

def write_out_result(directory, prefix, results):
	"""
	Writes out the file containing the optimum C,mu pairs

	Args:
		directory (string): Target directory for output files
		prefix (string): Prefix for output files. Output file will be named
			prefix.results
		results (list of tuples): List containing optimum C,mu pairs as well
			as the likelihood associated with those pairs
	"""

	filename = prefix + ".results"
	path = os.path.join(directory,filename)
	
	print "Writing results file to", path
	
	f = open(path, 'w')

	# Header
	f.write("#NLL\tmu\tC\tp*\n")

	for C,mu,L,vals in results:
		l_str = str(L) + "\t"
		mu_str = string.join([str(m) for m in mu],",") + "\t"
		m,n = C.shape
		C_str = ""
		for i in range(m):
			for j in range(1,n):
				C_str = C_str + str(int(C[i][j])) +","
			C_str = C_str[:-1]
			C_str += ":"
		C_str = C_str[:-1] + "\t"
		
		val_str = string.join([str(val) for val in vals],",")
		f.write(l_str)
		f.write(mu_str)
		f.write(C_str)
		f.write(val_str)
		f.write("\n\n")
	f.close()

def write_out_bounds(directory, prefix, inputFile, upper_bounds, lower_bounds):
	"""
	Writes out a copy of the input file with the bounds included

	Args:
		directory (string): Target directory for output files
		prefix (string): Prefix for output files. Output file will be named
			prefix.withBounds
		inputFile (string): full path to original input file
		upper_bounds (list of ints): array of upper bounds
		lower_bounds (list of ints): array of lower bounds

	"""
	f = open(inputFile)
	lines = f.readlines()
	f.close()

	outputFile = os.path.join(directory,prefix+".withBounds")
	f = open(outputFile,"w")


	print "Writing bounds file to", outputFile
	length = len(lines[1].split())
	
	if "#" in lines[0]: lines = lines[1:]
	
	# Header
	f.write("#ID\tchrm\tstart\tend\ttumorCount\tnormalCount\tUpperBound\tLowerBound\n")
	
	for i,line in enumerate(lines):
		f.write(line.strip())	
		if length == 6:				# Input file does not contain upper bounds
			  f.write("\t" + str(int(upper_bounds[i])))
		if length in [6,7]:			# Input file does not contain lower bounds
			  f.write("\t" + str(int(lower_bounds[i])))
		f.write("\n")

	f.close()
