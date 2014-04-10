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

# This code is to convert output from Exome Segmented data into THetA input.

import string
import os
import sys
import argparse
import subprocess


def parse_arguments():
	"""
	Parse all command line arguments.

	Returns:
		segments - full file path to segment file
		tumor - full file path to tumor bam
		normal - full file path to normal bam
		fasta - full file path to correct FASTA file
		exons - full file path to bed file of all exons
		directory - target directory for all output files
		prefix - prefix for output files
	"""

	parser = argparse.ArgumentParser()
	parser.add_argument("-s", "--SEGMENT_FILE", help="Segment file", metavar="s")
	parser.add_argument("-t", "--TUMOR_BAM", help="Tumor BAM", metavar="t")
	parser.add_argument("-n", "--NORMAL_BAM", help="Normal BAM", metavar="n")
	parser.add_argument("--OUTPUT_PREFIX", help="Output Prefix",metavar="OUTPUT_PREFIX",required=False)
	parser.add_argument("--DIR",help="Directory where results file is written.", \
						 default="./", metavar="DIR",required=False)
	parser.add_argument("--FA",help="Fasta File", metavar="FA")
	parser.add_argument("--EXON_FILE",help="Exon bed File", metavar="EXON_FILE")
	parser.add_argument("--QUALITY",help="Min Alignment Quality", metavar="QUALITY", default=30, \
						type=int, required=False)
	parser.add_argument("--EXCAVATOR",help="If running with EXCAVATOR segment file, \
						a file of chromosome ends must be provided.", metavar="EXCAVATOR", required=False)
	args = parser.parse_args()

	segment = args.SEGMENT_FILE
	tumor = args.TUMOR_BAM
	normal = args.NORMAL_BAM
	fasta = args.FA
	exons = args.EXON_FILE
	directory = args.DIR
	quality = args.QUALITY

	prefix = args.OUTPUT_PREFIX
	if prefix == None: prefix = os.path.basename(segment).split(".")[0]

	excavator = args.EXCAVATOR



	print "================================================="
	print "Arguments are:"
	print "\tSegment File:", segment
	print "\tTumor BAM:", tumor
	print "\tNormal BAM:", normal
	print "\tFASTA File:", fasta
	print "\tExon Bed File:", exons
	print "\tOutput Directory:", directory
	print "\tOutput Prefix:", prefix
	print "\tMin Alignment Quality:", quality
	if excavator is not None: print "\tChrm Ends File (for EXCAVATOR):", excavator
	print "================================================="

	return segment, tumor, normal, fasta, exons, directory, prefix, quality, excavator

def read_seg_file(segmented):

	f = open(segmented)
	lines = f.readlines()
	f.close()

	seg_data =[]

	for l in lines:

		#check for header
		if l.startswith("#"):
			continue

		line = l.strip().replace(" ","\t").split("\t")

		#Check for chrm format
		chrm = get_formatted_chrm(line[0])

		if chrm != -1:
			seg_data.append((chrm,int(line[1]),int(line[2])))

	return seg_data

def read_excavator_seg_file(segment, excavator):

	#Load Chrm ends
	chrm_ends = getChrmEnds(excavator)


	f = open(segment)
	lines = f.readlines()
	f.close()

	seg_data=[]
	chrm_start = 1
	prevChrm = 0
	prevPos = 1

	for l in lines:

		#skip header
		if l.startswith("#"):
			continue

		line = l.strip().replace(" ","\t").split("\t")
        
		chrm = get_formatted_chrm(line[0])
		
		#skip non-integer values		
		if chrm == -1:
			continue

		start = int(line[1])
		end = int(line[2])

		#Check new chrm
		if chrm != prevChrm and prevChrm != 0:

			end_of_chrm = chrm_ends[prevChrm]
			seg_data.append((prevChrm, prevPos, end_of_chrm))
			prevPos = chrm_start
			#prevChrm = chrm
			prevChrm+=1

		#Add any missing chromosomes
		while chrm != prevChrm and prevChrm != 0:
			end_of_chrm = chrm_ends[prevChrm]
			seg_data.append((prevChrm, chrm_start, end_of_chrm))
			prevChrm+=1
                        
		#Add previous interval if necessary
		if (start > prevPos):
			seg_data.append((chrm,prevPos, start-1))

		#Add current intervals
		seg_data.append((chrm, start, end))
		prevPos=end+1
		prevChrm = chrm

	#Add final segment
	end_of_chrm = chrm_ends[prevChrm]
	seg_data.append((prevChrm, prevPos, end_of_chrm))

	#Add missing chrms up to 22
	prevChrm+=1
	while prevChrm < 23:
		end_of_chrm=chrm_ends[prevChrm]
		seg_data.append((prevChrm, chrm_start, end_of_chrm))
		prevChrm+=1

	return seg_data



# Creates a map from chromosome number to the length of the chromosome
# for the provided build
def getChrmEnds(chrm_end_file):
	"""
	Creates a map from chomosome number to the length of the chromosome
	from the provided file.
	"""

	f = open(chrm_end_file)
	lines = f.readlines()
	f.close()

	chrm_ends = dict()
        
	for l in lines:

		#skip header
		if l.startswith("#"):
			continue

		line = l.strip().replace(" ","\t").split("\t")
		chrm = int(line[1])
		end = int(line[3])
		chrm_ends.update({chrm:end})

	return chrm_ends


def get_formatted_chrm(chr_string):
	"""
	Returns chromosome as an integer or -1 if invalid.
	"""

	#Check for start with chr
	if chr_string.lower().startswith("chr"):
		chr_string=chr_string[3:]

	#Check for X
	if chr_string.lower() == "x":
		chr_string = 23

	if chr_string.lower() == "y":
		chr_string = 24

	if chr_string.isdigit():
		return int(chr_string)
	else:
		return -1
		


def count_reads(seg_data, length, pileup, col):
	###
	# Counts the "number of reads" for the provided segmented intervals
	# for the provided pileup field in the specified column.
	###

	num_segs = len(seg_data)
	counts = [0]*num_segs
	cur_idx = 0
	cur_seg=seg_data[cur_idx]
	cur_chrm=cur_seg[0]
	cur_start=cur_seg[1]
	cur_end=cur_seg[2]

	#open pileup file
	with open(pileup,'r') as f:
		for line in f:

			vals=line.strip().replace(" ","\t").split("\t")
			
			#Process chromosome info, continue is not numeric
			chrm_str=vals[0]
			chrm=get_formatted_chrm(chrm_str)
			if chrm == -1:
				continue

			position=int(vals[1])
			count=int(vals[col])


			valid = 0
			done = 0

			while valid == 0:
				if cur_chrm > chrm:
					valid = 1
				elif chrm > cur_chrm: #moved beyond

					#check end
					if cur_idx == num_segs-1:
						done=1
						break

					#update
					cur_idx+=1
					cur_seg=seg_data[cur_idx]
					cur_chrm = cur_seg[0]
					cur_start = cur_seg[1]
					cur_end = cur_seg[2]
				elif cur_end >= position:
					valid = 1
				else:
					
					if cur_idx == num_segs-1:
						done=1
						break
					cur_idx+=1
					cur_seg=seg_data[cur_idx]
					cur_chrm = cur_seg[0]
					cur_start = cur_seg[1]
					cur_end = cur_seg[2]

			if done == 1:
				break
									

			#Check if contained in current interval
			if chrm == cur_chrm and cur_start <= position and position<= cur_end:
				counts[cur_idx] = counts[cur_idx] + count

	#normalize all counts by length
	for i in range(num_segs):
		curCount = counts[i]
		normCount = round(curCount/length)
		counts[i] = normCount

	return counts


def write_out_results(directory, prefix, seg_data, tumor, norm):

	outFile = os.path.join(directory,prefix+".input")
	f = open(outFile,"w")

	#Header
	f.write("#ID\tchrm\tstart\tend\ttumorCount\tnormalCount\n")

	#Enumerate and print lines
	for i,(chrm,start,end) in enumerate(seg_data):
		id="start_"+str(chrm)+"_"+str(start)+":end_"+str(chrm)+"_"+str(end)
		f.write(id+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(end)+"\t"+str(int(tumor[i]))+"\t"+str(int(norm[i]))+"\n")

	f.close()

def create_pileup(bam, exons, fasta, quality):
	"""
	Creates a pileup file for the specified bam file
	
	Args:
		bam - the bam file
		exons - the bed file with the list of exons
		fasta - the correct fasta build file
		quality - the minimum alignment quality to use

	Returns:
		pileup - the associated pileup file
		col - the column in the pileup file that contains coverage (diff for pileup vs mpileup)
	"""

	#determine samtools version (and quit if not installed)
	try:
		samtools = subprocess.check_output("which samtools",shell=True)
	except subprocess.CalledProcessError:
		print "Warning!  samtools is required."
		exit(1)
	

	pileup=os.path.abspath(bam).split(".bam")[0]+".pileup"
	col = 3

	
	e = os.system("samtools mpileup -f "+fasta+" -l "+exons+" -q "+str(quality)+" "+bam+" > "+ pileup)

	#if invalid try just pileup command
	if e != 0:

		#try pileup
		bam_quality=os.path.abspath(bam).split(".bam")[0]+".q"+str(quality)+".bam"
		e = os.system("samtools view -bh -q "+str(quality)+ " "+ bam + " > " + bam_quality)

		if e != 0:
			print "Wanring!  samtools unable to filter for quality."
			exit(1)


		e = os.system("samtools pileup -c -f "+fasta+" -r 0.0000007 -l "+exons+" "+bam+" > "+ pileup)

		if e != 0:
			print "Warning! samtools unable to make pileup file."
			exit(1)
		col=7

	return pileup, col

def get_read_length(bam):
	"""
	Returns the read length indicated by the first read in the bam file
	"""

	#check for samtools
	try:
		samtools = subprocess.check_output("which samtools",shell=True)
	except subprocess.CalledProcessError:
		print "Warning!  samtools is required."
		exit(1)

	#use subprocess to get the 10th column (index 9) of an entry.  Length is read length
	out1=subprocess.Popen(['samtools', 'view', bam], stdout=subprocess.PIPE)
	out2=subprocess.check_output(['head', '-n 1'] ,stdin=out1.stdout)
	length = len(out2.strip().replace(" ","\t").split("\t")[9])

	return length
		

def main():
	###
	# Read in arguments
	###
	segment, tumor, normal, fasta, exons, directory, prefix, quality, excavator  = parse_arguments()

	###
	# Create Pileup Files
	###
	pileup_tumor, t_col = create_pileup(tumor, exons, fasta, quality)
	pileup_normal, n_col = create_pileup(normal, exons, fasta, quality)

	###
	# Read in segmentation data
	###
	if excavator is None:
		seg_data = read_seg_file(segment)
	else:
		seg_data = read_excavator_seg_file(segment, excavator)

	###
	# Get Read Length
	###
	t_len=get_read_length(tumor)
	print("Tumor Read Length: " + str(t_len))
	n_len=get_read_length(normal)
	print("Normal Read Length: " + str(n_len))
	
	norm_reads = count_reads(seg_data, n_len, pileup_normal,n_col)
	tumor_reads = count_reads(seg_data, t_len, pileup_tumor, t_col)

	#directory="./"
	write_out_results(directory, prefix, seg_data, tumor_reads, norm_reads)


if __name__ == '__main__':
	main()
