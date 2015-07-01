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
import sys
import argparse

N_VALS = [2,3]
MAX_K = 7

def parse_arguments(silent=False):
	"""
	Parse command line arguments 
	
	Returns:
		query file: full path to the location of the input file
		results: for n=3 automatic interval selection, must provide results 
			of n=2 analysis
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
		num_processes: number of processes for THetA to use
		bounds_only: flag specifying to write out the bounds then exit
		time_estimate: flag to include the time estimate
		multi_event: flag to include rows with multi-events
		force: ignores certain warnings and forces THetA to run
		get_values: collects and prints out values for C, mu and likelihood for
			all Cs considered, for development purposes
		runBAF: flag indicating if the BAF post-processing model should be run.
		tumorSNP: file location for tumor SNP file used in the BAF post-processing model.
		normalSNP: file location for the normal SNP file used in the BAF post-processing model.
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
	parser.add_argument("--NUM_PROCESSES", help="The number of processes to be used",
			default=1, type=int, metavar="NUM_PROCESSES", required=False)
	parser.add_argument("--NUM_INTERVALS", help="The maximum number of intervals used by automatic interval selection.", default=100, type=int, metavar="NUM_INTERVALS", required=False)
	parser.add_argument("--BOUND_HEURISTIC", metavar="BH", default=False, required=False)
	parser.add_argument("--NORMAL_BOUND_HEURISTIC", metavar="NBH", type=int,\
			default=False,required=False)
	parser.add_argument("--HEURISTIC_LB", metavar="LB", type=float, default=0.9, \
			required=False)
	parser.add_argument("--HEURISTIC_UB", metavar="UB", type=float,default=1.1, \
			required=False)
	parser.add_argument("--BOUNDS_ONLY", action='store_true', default=False, required=False)
	parser.add_argument("--NO_MULTI_EVENT", action='store_true', default=False, required=False)
	parser.add_argument("--RESULTS", metavar = "filename", default=None, required=False)
	parser.add_argument("--FORCE", action = "store_true", default=False, required=False)
	parser.add_argument("--GET_VALUES", action = "store_true", default=False, required=False)
	parser.add_argument("--NO_INTERVAL_SELECTION", action = "store_true", default=False, required=False)
	parser.add_argument("--READ_DEPTH_FILE", metavar="FILENAME",  default=None, required=False)
	parser.add_argument("--GRAPH_FORMAT", help = "Options are .pdf, .jpg, .png, .eps" , default = ".pdf", required=False)
	parser.add_argument("--BAF", help="Option to run the BAF model.", , action='store_true', default=False, required=False, metavar="BAF")
	parser.add_argument("--TUMOR_SNP", help="File location for tumor SNP file used in the BAF model.", default=None, metavar="TUMOR_SNP", required=False)
	parser.add_argument("--NORMAL_SNP", help="File location for the normal SNP file used in the BAF model.", default=None, metavar="NORMAL_SNP", required=False)
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

	num_processes = args.NUM_PROCESSES
	bound_heuristic = args.BOUND_HEURISTIC
	normal_bound_heuristic = args.NORMAL_BOUND_HEURISTIC
	heuristic_lb = args.HEURISTIC_LB
	heuristic_ub = args.HEURISTIC_UB
	bounds_only = args.BOUNDS_ONLY
	multi_event = not(args.NO_MULTI_EVENT)
	results = args.RESULTS
	force = args.FORCE	
	get_values = args.GET_VALUES
	interval_selection = not(args.NO_INTERVAL_SELECTION)
	num_intervals = args.NUM_INTERVALS
	read_depth_file = args.READ_DEPTH_FILE
	graph_format = args.GRAPH_FORMAT
	runBAF = args.BAF
	tumorSNP = args.TUMOR_SNP
	normalSNP = args.NORMAL_SNP
	if n == 3 and num_intervals == 100: num_intervals = 20
	
	if not silent:
		print "================================================="
		print "Arguments are:"
		print "\tQuery File:", filename
		if n == 3 and results is not None: print "\tResults File:", results
		print "\tn:", n
		print "\tk:", k
		print "\ttau:", tau
		print "\tOutput Directory:", directory
		print "\tOutput Prefix:", prefix
		if n == 2: print "\tMax Normal:", max_normal
		if not(interval_selection): print "\tInterval Selection:", interval_selection
		if bound_heuristic is not False:
			print "\tBound Heuristic:", bound_heuristic
		if normal_bound_heuristic is not False:
			print "\tNormal Bound Heuristic:", normal_bound_heuristic
			print "\tHeuristic Lower Bound:", heuristic_lb
			print "\tHeuristic Upper Bound:", heuristic_ub
		print "\tNum Processes:", num_processes
		print "\tGraph extension:", graph_format
		if bounds_only: print "\tBounds Only:", bounds_only
		if force: print "\tForce:", force
		if get_values: print "\tGet Values:", get_values
		if read_depth_file is not None: print "Read depth file:", read_depth_file
		if runBAF:
			print "\t Tumor SNP File Location: ", tumorSNP
			print "\t Normal SNP File Location: ", normalSNP
		print "================================================="


	# NOTE: If you add an argument here, you must update both RunTHetA.py and GetPrefix.py to accept
	# the new argument
	return filename,results,n,k,tau,directory,prefix,max_normal,bound_heuristic, \
			normal_bound_heuristic, heuristic_lb, heuristic_ub, num_processes, \
			bounds_only, multi_event, force, get_values, interval_selection, \
			num_intervals, read_depth_file, graph_format, runBAF, tumorSNP, normalSNP
def read_binned_file(filename):
	"""
	Parses the data in a .binned.txt file.

	Arguments:
		filename (str): the location of the .binned.txt file

	Returns:
		data (2D list): A 2D list, where each row is [chromosome number, start interval,
						        	end interval, tumor counts,
								normal counts, count ratio, mean BAF,
								number SNPs]
	"""

	data = []

	with open(filename) as f:
		for line in f:
			if line.startswith("#"): continue

			chrm, start, end, tumorCounts, normalCounts, corrRatio, meanBAF, numSNPs = line.split("\t")

			chrm = int(chrm)
			start = int(start)
			end = int(end)
			tumorCounts = int(tumorCounts)
			normalCounts = int(normalCounts)
			corrRatio = float(corrRatio)
			meanBAF = float(meanBAF)
			numSNPs = int(numSNPs)

			if (corrRatio == -1) or (meanBAF == -1): continue

			data.append([chrm, start, end, tumorCounts, normalCounts, corrRatio, meanBAF, numSNPs])

	return data

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

	tumor_counts = []
	norm_counts = []
	upper_bounds = []
	lower_bounds = []
	lengths = []
	numLine = 0
	with open(filename) as f:
		for line in f:
			if line.startswith("#"): continue


			line = line.strip().replace(" ","\t").split()	
			numLine += 1
			
			if len(line) < 6 or len(line) > 8:
				sys.stderr.write("Invalid input file format in interval file line #"+str(numLine)+":\n" + str(line)+"\nToo few/many columns. Exiting...\n")
				sys.exit(1)
			# Read Lengths
			start = int(line[2])
			end = int(line[3])
			lengths.append(end-start)
		
			# Read Tumor Counts
			tumor_counts.append(int(line[4]))
			norm_counts.append(int(line[5]))

			# Read Bounds
			if len(line) > 6:
				upper_bounds.append(line[6])
			else: 
				upper_bounds.append("X")

			if len(line) > 7:
				lower_bounds.append(line[7])
			else: 
				lower_bounds.append("X")

	if numLine == 1:
		sys.stderr.write("Number of intervals must be greater than 1. Exiting...\n") 
		sys.exit(1)

	if all([x == "X" for x in upper_bounds]): upper_bounds = None	
	if all([x == "X" for x in lower_bounds]): lower_bounds = None	
	m = len(lengths)

	return (lengths, tumor_counts, norm_counts, m, upper_bounds, lower_bounds)

def read_interval_file_BAF(filename):
	"""
	Parses an interval file for information relevant to the BAF model.

	Args:
		filname: location of results file

	Returns:
		dataArray: a 2D array, with rows defined as (chromosome number, start position, end position).
				   Rows are stored as tuples.
	"""

	print "Reading interval file at " + filename
	chrmArray = []
	startPosArray = []
	endPosArray = []
	chrmsToUse = set()
	with open(filename) as f:
		for line in f:
			if line.startswith("#"): continue

			iden, chrm, startPos, endPos, tCount, nCount = line.strip().split("\t")
			chrmArray.append(int(chrm))
			chrmsToUse.add(int(chrm))
			startPosArray.append(int(startPos))
			endPosArray.append(int(endPos))

	dataArray = zip(chrmArray, startPosArray, endPosArray)
	chrmsToUse = list(chrmsToUse)
	return chrmsToUse, dataArray

def read_results_file(filename):
	"""
	For n=3 with automatic interval selection, reads in the results file to get C
	Args:
		filename: location of results file
	Returns:
		C: list form of the second column of the result copy number profile
	"""

	with open(filename) as f:
		lines = f.readlines()
	if lines[0].startswith("#"):
		lines = lines[1:]
	if len(lines) == 0:
		print "ERROR: The result file provided appears to be empty. Exiting..."
	elif len(lines) > 1:
		print "WARNING: The results file contains more than one solution. THetA will use the first provided solution."
	
	soln = lines[0].strip().split("\t")
	copy = [i for i in soln[2].split(":")]
	return copy

def read_results_file_full(filename):
	"""
	Parses a results file.
	Args:
		filename: location of results file
	Returns:
		results: a dictionary of the results with the following keys:
			-NLL: an array of negative log likelihoods.
			-mu: an array of mu vectors
			-C: an array of C matrices
			-p: an array of p* vectors
			-k: number of solutions contained in results file
	"""

	print "Reading results file at " + filename
	negLLArray = []
	muArray = []
	cMatArray = []
	pArray = []
	k = 0
	with open(filename) as f:
		for line in f:
			if line.startswith("#"): continue

			negLL, mu, c, p = line.strip().split("\t")
			negLLArray.append(float(negLL))

			mu = map(float, mu.split(","))
			muHead = [mu[0]]
			muTail = mu[1:]
			n = len(muTail)
			muTail = zip(range(n), muTail)
			muTail = sorted(muTail, key=lambda x: x[1], reverse=True)
			ind, muTail = zip(*muTail)
			mu = muHead + list(muTail)
			muArray.append(mu)

			c = c.split(":")
			c = map((lambda x: x.split(",")), c)
			for i in range(len(c)):
				temp = range(n + 1)
				if c[i][0] == "X":
					temp = [-1] * (n + 1)
				else:
					temp[0] = 2
					for j in range(n):
						temp[j + 1] = int(c[i][ind[j]])
				c[i] = temp
			cMatArray.append(c)

			p = map(lambda x: -1 if x == "X" else float(x), p.split(","))
			pArray.append(p)

			k += 1

	results = {'NLL': negLLArray, 'mu': muArray, 'C': cMatArray, 'p': pArray, 'k': k}
	return results



def read_snp_file(filename, option="normal"):
	"""
	Converts an SNP file to a 2D array.

	Args:
		filename: location of SNP file
	Returns:
		data: data from SNP file in a 2D array, where each row is [chromosome, position, ref count, mut count]
	"""

	chrmInd = 0
	posInd = 1
	if option == "het":
		refInd = 2
		mutInd = 3
	else:
		refInd = 7
		mutInd = 8

	print "Reading SNP file at " + filename

	data = []
	with open(filename) as f:
		for line in f:
			if line.startswith("#"): continue

			vals = line.split("\t")

			if vals[0] == "X":
				chrm = 23
			elif vals[0] == "Y":
				chrm = 24
			else:
				chrm = int(vals[chrmInd])

			position = int(vals[posInd])
			refCount = float(vals[refInd])
			mutCount = float(vals[mutInd])
			data.append([chrm, position, refCount, mutCount])

	return data

def write_out_result(directory, prefix, results, n):
	"""
	Writes out the file containing the optimum C,mu pairs

	Args:
		directory (string): Target directory for output files
		prefix (string): Prefix for output files. Output file will be named
			prefix.results
		results (list of tuples): List containing optimum C,mu pairs as well
			as the likelihood associated with those pairs
	"""

	filename = prefix + ".n"+str(n)+".results"
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
				if int(C[i][j]) == -1: 
					C_str = C_str + "X" + ","
				else:
					C_str = C_str + str(int(C[i][j])) +","
			C_str = C_str[:-1]
			C_str += ":"
		C_str = C_str[:-1] + "\t"
		
		val_str = string.join([str(val) for val in vals],",")
		f.write(l_str)
		f.write(mu_str)
		f.write(C_str)
		f.write(val_str)
		f.write("\n")
	f.close()

def write_out_NLL_result(directory, prefix, results):
	"""
	Writes out the file containing the results from the BAF model

	Args:
		directory (string): Target directory for output files
		prefix (string): Prefix for output files. Output file will be named
			prefix.BAF.NLL.results
		results (dictionary): dictionary of results produced by runBAFModel.py
	"""
	NLL = results['NLL']
	mu = results['mu']
	C = results['C']
	p = results['p']
	BAF_NLL = results['BAF_NLL']

	filename = prefix + ".BAF.NLL.results"
	path = os.path.join(directory, filename)

	print "Writing results file to", path

	f = open(path, 'w')

	f.write("#NLL\tmu\tC\tp*\tBAF_NLL\n")
	for i in range(results['k']):
		to_csv = lambda x: ",".join(map(lambda y: str(y) if y != -1 else "X", x))

		currNLL = NLL[i]
		NLLstr = str(currNLL)
		f.write(NLLstr + "\t")

		currMu = mu[i]
		muStr = to_csv(currMu)
		f.write(muStr + "\t")

		currC = C[i]
		commaSV = map(lambda x: to_csv(x[1:]), currC)
		colonSV = ":".join(commaSV)
		f.write(colonSV + "\t")

		currP = p[i]
		pStr = to_csv(currP)
		f.write(pStr + "\t")

		currBAF_NLL = BAF_NLL[i]
		BAF_NLLStr = str(currBAF_NLL)
		f.write(BAF_NLLStr + "\n")

	f.close()

def write_out_bounds(directory, prefix, inputFile, upper_bounds, lower_bounds, n, order=None):
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

	outputFile = os.path.join(directory,prefix+".n"+str(n)+".withBounds")
	f = open(outputFile,"w")

	print "Writing bounds file to", outputFile
	length = len(lines[1].split())
	
	if "#" in lines[0]: lines = lines[1:]
	
	# Header
	f.write("#ID\tchrm\tstart\tend\ttumorCount\tnormalCount\tUpperBound\tLowerBound\n")

	if order is not None:
		orderMap = {}
		for i, v in enumerate(order):
			orderMap[v] = i
		for i,line in enumerate(lines):
			line = "\t".join(line.strip().split("\t")[:6])
			f.write(line.strip())	
			if i in orderMap:
				f.write("\t" + str(int(upper_bounds[orderMap[i]])))
				f.write("\t" + str(int(lower_bounds[orderMap[i]])))
			else: 
				f.write("\tX")
				f.write("\tX")
			f.write("\n")
	else:	
		for i,line in enumerate(lines):
			line = "\t".join(line.strip().split("\t")[:6])
			f.write(line.strip())	
			f.write("\t" + str(int(upper_bounds[i])))
			f.write("\t" + str(int(lower_bounds[i])))
			f.write("\n")
	f.close()

def write_out_N3_script(directory, prefix, inputFile):
	# All the arguments are the same except -n 3 instead of 2 and --bounds prefix.n2.withbounds	

	filename = os.path.join(directory, prefix+".RunN3.bash")
	print "Writing script to run N=3 to ", filename
	with open(filename, 'w') as f:
		argString = " ".join(sys.argv)
		boundsFile = os.path.join(directory, prefix + ".n2.withBounds")
		resultsFile = os.path.join(directory, prefix + ".n2.results")

		string = "python "+ argString.replace("-n 2", "").replace(inputFile, boundsFile) +" -n 3" + " --RESULTS " + resultsFile
		f.write("#!/bin/bash\n")
		f.write(string)
