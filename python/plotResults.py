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
 # http://cs.brown.edu/people/braphael/software.html
 # 
 # @author Layla Oesper, Ahmad Mahmoody, Benjamin J. Raphael, Gryte Satas, David Liu
 ###

import matplotlib
matplotlib.use('Agg')
import sys, traceback
import matplotlib.pyplot as plt
plt.ioff()
import csv
from collections import defaultdict
import os

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Given a results, interval, and concordant read file, this script outputs the graph as a pdf.

Usage:

python <results_file> <interval_file> <concordant_file> <prefix> <n_subpops>

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Read input results, interval, concordant files
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def plot_results(out_dir, filename, prefix, concordant_file, n_subpops, extension):
	#results file
	results_file = os.path.join(out_dir,prefix + ".n"+str(n_subpops)+".results")
	results_path = os.path.abspath(results_file)

	#interval file
	interval_file = os.path.join(out_dir,prefix+".n"+str(n_subpops)+".withBounds")
	interval_path = os.path.abspath(interval_file)

	concordant_path = None
	if concordant_file:
		#concordant file
		concordant_path = os.path.abspath(concordant_file)

	#Output file
	output_filename = prefix + ".n" + str(n_subpops) + ".graph" + extension
	output_path = os.path.abspath(os.path.join(out_dir, output_filename))

	sample = prefix

	#Dictionary of (k,v): (Chromsome#: (intervalstart, intervalend))
	intervals = defaultdict(list)
	#Dictionary of (k,v): (Chromosome#: length)
	chromosome_lengths = dict()

	#Interval File
	#Read in interval file only once
	with open(interval_path, "r+") as interval_file:
		#Store the bounds data. Note that these are both 1-indexed. Would work for x/y chromosome, as well. Also note that indices are strings.
		columns = defaultdict(list)
		reader = csv.DictReader(interval_file, delimiter = "\t") # read rows into a dictionary format
		for row in reader: # read a row as {column1: value1, column2: value2,...}
		    for (k,v) in row.items(): # go over each column name and value 
		        columns[k].append(v) # append the value into the appropriate list based on column name k

		#Read in the intervals/chromosomes
		for i in range(len(columns["#ID"])):
			intervals[columns["chrm"][i]].append((int(columns["start"][i]), int(columns["end"][i])))

		for i in range(1, len(intervals) + 1):
			i = str(i)
			chromosome_lengths[i] = intervals[i][-1][1]


	#For plotting purposes.
	names = chromosome_lengths.keys()
	chromosome_names = []
	for name in names:
		try:
			chromosome_names.append(int(name))
		except:
			chromosome_names.append(name)
	chromosome_names.sort()

	lengths = [chromosome_lengths[str(i)] for i in chromosome_names]
	chromosome_cummulative = [sum(lengths[:i]) for i in range(1, len(chromosome_lengths) + 1)]

	minor_locations = []
	for i,c in enumerate(chromosome_names):
		minor_locations.append(chromosome_cummulative[i] - chromosome_lengths[str(c)]/2)

	
	#Concordant files
	#Dictionary of (k,v): (Chromosome#: [( (start + end)/2, tumor_count, normal_count)])
	bins = defaultdict(list)
	#Totals, reset on each chromosome
	tumor_total = 0
	normal_total = 0
	#Dictionary of (k,v): (Chromosome#: (tumor total, normal total))
	totals = dict()
	num_plots = 0
	current_chrm = '1'

	#If the read depth file exists
	if concordant_path:
		with open(concordant_path, "r+") as concordant_read_file:
			reader = csv.reader(concordant_read_file, delimiter = "\t")
			num_plots = len(concordant_read_file.readline()) - 6
			for row in reader:
				if row[1] != current_chrm:
					totals[current_chrm] = (int(tumor_total), (normal_total))
					current_chrm = row[1]
					tumor_total = 0
					normal_total = 0
				bins[current_chrm].append( ((int(row[2]) + int(row[3]))/2 , int(row[4]), int(row[5])) )
				tumor_total += int(row[4])
				normal_total += int(row[5])
			#At the very end
			totals[current_chrm] = (tumor_total, normal_total)

	def make_subplot(a_plt, line, number):

		parts = line.split("\t")
		#Store the copy number data
		#The first part doesn't seem to mean anything.
		del parts[0]
		mu = parts[0].split(",")
		num_subpop = len(mu) - 1

		#Split the comma-seperated values.
		def split_seperate(alist):
			#returnArray = [[]] * num_subpop
			returnArray = []
			for i in range(num_subpop):
				returnArray.append([])
			for entry in alist:
				values = entry.split(",")
				for n, value in enumerate(values):
					returnArray[n].append(int(value))
			return returnArray

		#C is a list of lists. Index 0 maps to subpopulation 1.
		C = []

		if num_subpop ==1:
			C = [[int(i) for i in parts[1].split(":")]]
		elif num_subpop > 1:
			C = split_seperate(parts[1].split(":"))

		#p = [float(p) for p in parts[2].split(",")]

		max_copy_number = 0
		max_copy_number = max(max(pop) for pop in C)


		"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		Make axes
		"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

		#Actual Plot
		ax = fig.add_subplot(len(lines), 1, number + 1)
		#Axes Labels
		ax.set_xlabel("Chromosome")
		ax.set_ylabel("Copy Number")

		#Adjust ticks
		xaxis = ax.get_xaxis()
		xaxis.set_ticklabels(chromosome_names, minor = True)
		xaxis.set_ticklabels([])
		xaxis.set_ticks(chromosome_cummulative)
		xaxis.set_ticks(minor_locations, minor = True)
		xaxis.set_tick_params(which = 'minor', labelsize = 8)
		xaxis.grid(True, which = 'major', linestyle = "-")
		yaxis = ax.get_yaxis()
		yaxis.set_tick_params(size = 0)
		ax.set_ylim(0, 6)
		ax.set_xlim(0, chromosome_cummulative[-1])


		#Title
		subtitle = "Normal:" + str(round(float(mu[0]) * 100, 1)) + "%, "
		for i in range(num_subpop):
			subtitle += "Tumor" + str(i + 1) + ":" + str(round(float(mu[i + 1]) * 100, 1)) + "%"
			if i != num_subpop - 1:
				subtitle += ", "
		ax.set_title(subtitle)
		# plt.text(0.5, 1.08, subtitle,
  #        horizontalalignment='center',
  #        fontsize=20,
  #        transform = ax.transAxes)


		#Colors, assuming there are fewer than 9 subpopulations
		#black = k
		colors = "brgcmyw"

		labels = ["Normal"]
		#Make labels
		for i in range(num_subpop):
			labels.append("Tumor " + str(i + 1))


		"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		Plot grey dots: tumor to normal read ratio
		"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		if concordant_path:
			names = [ ]
			namez = totals.keys()
			for name in namez:
				try:
					names.append(int(name))
				except:
					names.append(name)
			names.sort()
			for n, name in enumerate(names):
				this_chrom_bins = bins[str(name)]
				this_tumor_total = totals[str(name)][0]
				this_normal_total = totals[str(name)][1]
				x = [] #Interval midpoints
				y = [] #Ratios
				for bin in this_chrom_bins:
					if n == 0:
						x.append(bin[0])
					else:
						x.append(bin[0] + chromosome_cummulative[n - 1])
					try:
						y.append(2 * (bin[1]/float(this_tumor_total))/(bin[2]/float(this_normal_total)) )
					except:
						del x[-1]
						continue
				ax.scatter(x,y, marker = '.', facecolor='0.75', lw = 0, s = 5)


		"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		Plot Copy Numbers
		"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

		#Plot reference data

		interval_num = 0
		x1 = 0
		x2 = 0

		for n, name in enumerate(chromosome_names):
			name = str(name)
			chrom_intervals = intervals[name]
			for i in chrom_intervals:
				if i[1] - i[0] < 10000:
					interval_num += 1
					continue
				else: 
					if n == 0:
						x1 = i[0]
						x2 = i[1]
					else:
						x1 = i[0] + chromosome_cummulative[n - 1]
						x2 = i[1] + chromosome_cummulative[n - 1]
					if interval_num == 0:
						ax.plot((x1,x2), (2,2), color = 'k', linewidth = 3, label = labels[0], solid_capstyle = "butt")
					ax.plot((x1, x2), (2, 2), color = 'k', linewidth = 3, solid_capstyle = "butt")
					interval_num += 1



		#Plot tumor bars

		interval_num = 0
		x1 = 0
		x2 = 0

		for n, name in enumerate(chromosome_names):
			name = str(name)
			chrom_intervals = intervals[name]
			for i in chrom_intervals:
				#Don't plot intervals with less than 100000bp
				if i[1] - i[0] < 10000:
					interval_num += 1
					continue
				else:
					if n == 0:
						x1 = i[0]
						x2 = i[1]
					else:
						x1 = i[0] + chromosome_cummulative[n - 1]
						x2 = i[1] + chromosome_cummulative[n - 1]
					for j in range(num_subpop):
						j = j + 1
						# + 0.1 for no overlap.
						#Add the label only once.
						if interval_num == 0:
							ax.plot((x1, x2), (C[j - 1][interval_num] + 0.10 *j, C[j - 1][interval_num] + 0.10 *j), color = colors[j - 1], linewidth = 3, label = labels[j], solid_capstyle = "butt")
						else:
							ax.plot((x1, x2), (C[j - 1][interval_num] + 0.10 *j, C[j - 1][interval_num] + 0.10 *j), color = colors[j - 1], linewidth = 3, solid_capstyle = "butt")
					interval_num += 1



		#Legend
		# Shink current axis by 5%
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.95, box.height])
		ax.legend(loc='center left', bbox_to_anchor=(1.03, 0.5), prop = {'size' : 8}, borderpad=1.5, labelspacing=1.5)





	"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	Main method + Show
	"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


	with open(results_path, "r+") as results_file:
		global lines
		lines = results_file.readlines()
		#Skip header
		del lines[0]

		#Only one plot, format this way.
		if len(lines) == 1:

			fig = plt.figure(facecolor='w', dpi = 150, edgecolor='k', figsize=(12, len(lines) * 3))
			# fig, axarr = plt.subplots(len(lines))
			title_text = sample
			fig.suptitle(title_text, fontsize = 16, x = 0.45)		
			for i, result in enumerate(lines):
				make_subplot(fig, result, i)

			plt.subplots_adjust(hspace = 0.4, left = 0.05, right = 0.85, top = 0.82, bottom = 0.15)

		#More than one plot.
		elif len(lines) > 1:

			fig = plt.figure(facecolor='w', dpi = 150, edgecolor='k', figsize=(12, len(lines) * 3))
			# fig, axarr = plt.subplots(len(lines))
			title_text = sample
			fig.suptitle(title_text, fontsize = 16, x = 0.45)	
			for i, result in enumerate(lines):
				make_subplot(fig, result, i)

			plt.tight_layout()

			plt.subplots_adjust(hspace = 0.45, left = 0.05, right = 0.85, top = 0.86, bottom = 0.15)

		
	plt.savefig(output_path)


#plot_results(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], ".pdf")
