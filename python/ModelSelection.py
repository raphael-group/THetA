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
import math 
import os.path

inputFile = sys.argv[1]
n2Result = sys.argv[2]
n3Result = sys.argv[3]

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

if P_NLL_N2 < P_NLL_N3:
	filename = n2Result.replace(".n2.results", ".BEST.results")
	print "Selected n=2 solution. Writing result to", filename,
	with open(filename, 'w') as out, open(n2Result) as n2In:
		for line in n2In:
			out.write(line)
	pdfFileN2 = n2Result + ".pdf"
	pdfFileBest = filename + ".pdf"
	if os.path.isfile(pdfFileN2):
		call(["cp", pdfFileN2, pdfFileBest])
		print ",", pdfFileBest
	else: print ""
else:
	filename = n3Result.replace(".n3.results", ".BEST.results")
	print "Selected n=3 solution. Writing result to", filename,
	with open(filename, 'w') as out, open(n3Result) as n3In:
		for line in n3In:
			out.write(line)
	pdfFileN3 = n3Result + ".pdf"
	pdfFileBest = filename + ".pdf"
	if os.path.isfile(pdfFileN3):
		call(["cp", pdfFileN3, pdfFileBest])
		print ",", pdfFileBest
	else: print ""
