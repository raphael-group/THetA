#!/usr/bin/python

import argparse
import gzip
import sys
import os

def parse_arguments():
	"""
	Parse command line arguments.

	Returns:
		input_filename: path to the location of the input file.
		output_filename: path to the location of the output file.
	"""

	parser = argparse.ArgumentParser()
	parser.add_argument("BAF_FILE", help="File containing BAF data.", metavar="BAF_FILE")
	parser.add_argument("-p", help="Prefix for output file. By default, it will \
		be the beginning of the input file name.", default=None, metavar="p", required=False)
	parser.add_argument("-o", help="Directory where result file is written.", 
		metavar="o", default="./", required=False)
	args = parser.parse_args()

	input_filename = args.BAF_FILE
	
	output_prefix = args.p
	if output_prefix is None:
		output_prefix = os.path.basename(input_filename).split(".")[0]
	
	output_filename = args.o + output_prefix + ".withCounts"

	return input_filename, output_filename

def convert_file(input_filename, output_filename):
	"""
	Converts a CSV BAF file to a TSV BAF file.

	Args:
		input_filename: path to the location of the input file.
		output_filename: path to the location of the output file.
	"""

	suffix = os.path.basename(input_filename).split(".")[-1]
	try:
		if suffix == "gz":
			f = gzip.open(input_filename)
		else:
			f = open(input_filename)
	except IOError:
		print "An error occured while opening the input file. Exiting program..."
		sys.exit(1)

	try:
		o = open(output_filename, "w")
	except IOError:
		print "An error occured while opening the output file. Exiting program..."
		print output_filename
		sys.exit(1)

	o.write("#Chrm\tpos\tA\tC\tG\tT\ttotal\trefCount\tmutCount\n")
	linenum = 0
	for line in f:
		linenum += 1
		if not line == "\n":
			vals = line.split(",")
			if len(vals) != 5:
				print "Invalid input file; insufficient number of values at line %i. Exiting program..." % linenum 
				os.remove(output_filename)
				sys.exit(1)
			chrm, pos, refCount, mutCount, val = vals
			o.write("%s\t%s\t0\t0\t0\t0\t0\t%s\t%s\n" % (chrm, pos, refCount, mutCount))
	o.close()

if __name__ == "__main__":
	input_filename, output_filename = parse_arguments()
	convert_file(input_filename, output_filename)
