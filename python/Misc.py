 ###
 # Copyright 2012, 2013, 2014, 2015, 2017 Brown University, Providence, RI.
 #
 # All Rights Reserved
 # 
 # Permission to use this software, and any documentation, for non-commercial academic research
 # purposes only is hereby granted with the following terms and conditions:
 # 
 # (1) the above copyright notice and this permission notice shall be preserved in all instances
 # of the software and in any supporting documentation;
 # 
 # (2) the name of Brown University shall not be used in advertising or publicity pertaining 
 # to the use of the software without specific, written prior permission;
 # 
 # (3) the rights granted herein are individual and personal to the recipient and may not be
 # sublicensed or distributed to any third party without specific, written prior permission; and
 # 
 # (4) the permitted user acknowledges that all commercial rights are licensed to Medley
 # Genomics, Inc., and any inquiries related to commercial use shall be directed to Medley
 # Genomics, Inc.
 # 
 # BROWN UNIVERSITY PROVIDES THIS SOFTWARE AND ANY DOCUMENTATION
 # “AS IS” AND DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE
 # AND ANY DOCUMENTATION, INCLUDING ALL IMPLIED WARRANTIES OF
 # MERCHANTABILITY AND FITNESS FOR ANY PARTICULAR PURPOSE. IN NO
 # EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR ANY SPECIAL, INDIRECT OR
 # CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 # LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
 # NEGLIGENCE OR OTHER ACTION BASED ON ANY OTHER LEGAL THEORY,
 # ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
 # SOFTWARE.
 # 
 # @author Layla Oesper, Ahmad Mahmoody, Benjamin J. Raphael and Gryte Satas
 ###

def isClose(v1, v2, margin = 10e-4):
	"""
	Check if every pair of floats (v1[i],v2[i] are within the margin value of eachother)
	
	Args:
		v1, v2 (list of floats): lists with len(v1) = len(v2)
		margin (float): margin distance to judge "closeness"
	"""
	for (a,b) in zip(v1,v2):
		if abs(a-b) > margin: 
			return False
	return True

def inRange(v1, minVal = 0, maxVal = 1):
	"""
	Check if every value in v1 is in the range (minVal, maxVal)
	"""
	for v in v1: 
		if v < minVal or v > maxVal: 
			return False

	return True
