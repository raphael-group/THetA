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


import sys
import os.path
import math
from FileIO import read_interval_file, load_results


def ModelSelection(inputFile, n2Result, n3Result):

	###
	#	Read Result File to Obtain #Tumor Reads, #Normal Reads, #Intervals
	###
	numTumor = 0
	numNormal = 0
	numIntervals = 0
	with open(inputFile) as f:
		for line in f:
			if line.startswith("#"): continue
			tumor, normal = line.strip().split("\t")[4:6]

			if int(normal) > 0:
				numTumor += int(tumor)
				numNormal += int(normal)
				numIntervals += 1

	###
	#	Obtain NLL from n = 2
	###
	n2NLL= float("inf")
	with open(n2Result) as f:
		for line in f:
			if line.startswith("#"): continue
			likelihood = float(line.strip().split("\t")[0])
			if likelihood < n2NLL: n2NLL = likelihood

	###
	#	Obtain NLL from n = 3
	###
	n3NLL= float("inf")
	with open(n3Result) as f:
		for line in f:
			if line.startswith("#"): continue
			likelihood = float(line.strip().split("\t")[0])
			if likelihood < n3NLL: n3NLL = likelihood


	###
	#	Penalized NLLs
	#	P-NLL = 2*L + (m+1)(n-1)log(T+N)
	###
	P_NLL_N2 = 2*n2NLL + (numIntervals+1) * math.log(numTumor + numNormal)
	P_NLL_N3 = 2*n3NLL + (numIntervals+1) * 2 * math.log(numTumor + numNormal)

	from subprocess import call

	selected_num=2
	selected_res=n2Result
	if P_NLL_N3 <= P_NLL_N2:
		selected_num, selected_res=additional_criteria(n2Result, n3Result, inputFile)

	postfix=".n"+str(selected_num)+".results"
	filename = selected_res.replace(postfix, ".BEST.results")
	print "Selected n="+str(selected_num)+" solution.  Writing to",filename
	with open(filename, 'w') as out, open(selected_res) as nIn:
		for line in nIn:
			out.write(line)
	pdfFileN = selected_res + ".pdf"
	pdfFileBest = filename + ".pdf"
	if os.path.isfile(pdfFileN):
		call(["cp", pdfFileN, pdfFileBest])
		print ",", pdfFileBest
	else: print ""

def additional_criteria(n2Result, n3Result, inputFile, min_pop=0.05, min_clonal=0.0, max_ratio=5, min_ratio=0.05):
	###
	# Will determine if the n=3 result passes a set of additional criteria including:
	# (1) All tumor pops >= min_pop
	# (2) Clonal population > min_clonal
	# (3) Subclonal/Clonal ratio <= max_ratio
	# (4) Subclonal/Clonal ratio >= min_ratio
	# (5) If multiple solutions exist, at least on must pass criteria (1)-(4).
	###

	selected_num=2
	selected_res=n2Result

	intervals = read_interval_file(inputFile)
	lengths, _, _, _, _, _ = intervals

	results = load_results(n3Result)
	isValid = False

	for nll, C, mu in results:
		noCNA, clonal, subclonal = get_frac_breakdown(C,lengths)

		pop_is_big = all(i > min_pop for i in mu[1:])
		clonal_is_big = clonal > min_clonal

		ratio_is_small_enough = False
		ratio_is_big_enough = True
		if clonal > 0:
			ratio = float(subclonal)/float(clonal)
			ratio_is_small_enough = ratio < max_ratio

			if ratio < min_ratio: ratio_is_big_enough = False

		if pop_is_big and clonal_is_big and ratio_is_small_enough and ratio_is_big_enough:
			isValid = True

	if isValid:
		selected_num=3
		selected_res=n3Result

	return selected_num, selected_res

def get_frac_breakdown(C,lengths):
	"""
	Returns a vector with the breakdown of the genome indicator
	three categories:
	1. No copy number aberrations.
	2. Clonal aberrations.
	3. Subclonal aberrations.

	Note, this may not sum to 1 as there may be portions of the genome where we cannot call
	copy numbers at all.
	"""
	tot_len = 0
	tot_norm = 0
	tot_clonal = 0
	tot_subclonal = 0

	for i,row in enumerate(C):
		cur_len = lengths[i]
		tot_len += cur_len

		row_str = [str(x) for x in row[1:]]

		if checkEqual(row_str,'X'):
			continue
		elif checkEqual(row_str,'2'):
			tot_norm += cur_len
			continue
		elif checkEqual(row_str,row_str[0]):
			tot_clonal += cur_len
		else:
			tot_subclonal += cur_len

	norm_frac = float(tot_norm)/float(tot_len)
	clonal_frac = float(tot_clonal)/float(tot_len)
	subclonal_frac = float(tot_subclonal)/float(tot_len)

	return norm_frac, clonal_frac, subclonal_frac

def checkEqual(vec, val):
	"""
	Checks if all entries are equal to the supplied value.
	"""
	allEqual = True

	for v in vec:
		if v != val: return False

	return allEqual

if __name__ == '__main__':
	inputFile = sys.argv[1]
	n2Result = sys.argv[2]
	n3Result = sys.argv[3]

	ModelSelection(inputFile, n2Result, n3Result)
