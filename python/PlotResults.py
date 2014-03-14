###
# 2014 Brown University, Providence, RI.
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
# http://cs.brown.edu/software/
# 
# @author Layla Oesper, Ahmad Mahmoody, Benjamin J. Raphael and Gryte Satas
###

import string
import os
import sys
import argparse
import subprocess

from FileIO import *


def plot_results(out_dir, interval_file, n, prefix):
	"""
	Create results plots if matlab exists on system.
	
	Args:
		out_dir (string) - the output directory.
		interval_file (string) - the interval file used as input.
		n (int) - the number of populations to infer.
		prefix (string) - the output prefix used.
	"""
	
	#Check for matlab
	try:
		matlab = subprocess.check_output("which matlab",shell=True)
	except subprocess.CalledProcessError:
		print "Warning!  Matlab is required to make results plots."
		exit(1)
	

	#.m file
	plot_filename = "plotResults.m"
	plot_path = os.path.abspath(os.path.join(out_dir, plot_filename))

	#interval file
	interval_path = os.path.abspath(interval_file)

	#results file
	results_filename = prefix + ".n" + str(n) + ".results"
	results_path = os.path.abspath(os.path.join(out_dir, results_filename))

	#matlab directory
	matlab_dir=os.path.join(os.path.dirname(__file__),'..','matlab','mainMethod')

	#Check for results and interval files
	if os.path.isfile(results_path) and os.path.isfile(interval_path):
	
		fid = open(plot_path, 'w')

		fid.write("path(path, '"+matlab_dir+"')\n")
		fid.write("plotResultsFromFile('"+results_path+"','"+interval_path+"','"+prefix+"',12,3)\n")
		fid.write("exit")
		fid.close()

		#Call matlab code from here
		var = subprocess.check_output("stty -g", shell=True) #save old configuration
		command = 'matlab -nodisplay -nosplash -r "run '+plot_path+'"'
		os.system(command)
		os.system('stty '+var) #reset the old configuration
		print('\n') 

	os.remove(plot_path)
	
	

if __name__=='__main__':
	main()

