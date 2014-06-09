import sys, traceback
import matplotlib.pyplot as plt
import csv
from collections import defaultdict
import os

__author__ = "David Liu"

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Given a results file and a bound file, this script renders the copy number graph as a pdf.

Usage:

python <results_file> <results_bound_file>

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Read input interval file
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def plot_results(out_dir, interval_file, prefix, n_subpops):

	#results file
	results_filename = prefix + ".n" + str(n_subpops) + ".results"
	results_path = os.path.abspath(os.path.join(out_dir, results_filename))

	#interval file
	interval_path = os.path.abspath(interval_file)

	#Output file
	output_filename = prefix + ".n" + str(n_subpops) + ".graph.pdf"
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



	def make_subplot(a_plt, line, number):

		parts = line.split("\t")

		#Store the copy number data
		#The first part doesn't seem to mean anything.
		del parts[0]
		mu = parts[0].split(",")
		num_subpop = len(mu) - 1

		#Split the comma-seperated values.
		def split_seperate(alist):
			returnArray = [[]] * num_subpop
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

		p = [float(p) for p in parts[2].split(",")]

		max_copy_number = 0
		max_copy_number = max(max(pop) for pop in C)


		"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		Make axes
		"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


		#Actual Plot
		ax = fig.add_subplot(2, 1, number + 1)
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


		#Title
		subtitle = "Normal:" + str(round(float(mu[0]) * 100, 1)) + "%, "
		for i in range(num_subpop):
			subtitle += "Tumor" + str(i + 1) + ":" + str(round(float(mu[i + 1]) * 100, 1)) + "%"
			if i != num_subpop - 1:
				subtitle += ", "
		ax.set_title(subtitle)


		#Colors, assuming there are fewer than 9 subpopulations
		#black = k
		colors = "bgrcmyw"

		labels = ["Normal"]
		#Make labels
		for i in range(num_subpop):
			labels.append("Tumor " + str(i + 1))


		"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		Plot Stuff
		"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

		#Plot reference data

		interval_num = 0
		x1 = 0
		x2 = 0

		for n, name in enumerate(chromosome_names):
			name = str(name)
			chrom_intervals = intervals[name]
			for i in chrom_intervals:
				if i[1] - i[0] < 100000:
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
						ax.plot((x1,x2), (2,2), color = 'k', linewidth = 2, label = labels[0], solid_capstyle = "butt")
					ax.plot((x1, x2), (2, 2), color = 'k', linewidth = 2, solid_capstyle = "butt")
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
				if i[1] - i[0] < 100000:
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
							ax.plot((x1, x2), (C[j - 1][interval_num] + 0.09 *j, C[j - 1][interval_num] + 0.09 *j), color = colors[j - 1], linewidth = 2, label = labels[j], solid_capstyle = "butt")
						else:
							ax.plot((x1, x2), (C[j - 1][interval_num] + 0.09 *j, C[j - 1][interval_num] + 0.09 *j), color = colors[j - 1], linewidth = 2, solid_capstyle = "butt")
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
		lines = results_file.readlines()
		#Skip header
		del lines[0]
		fig = plt.figure(facecolor='w', edgecolor='k', figsize=(12, len(lines) * 3))
		title_text = sample
		fig.suptitle(title_text, fontsize = 16, x = 0.45)		
		for i, result in enumerate(lines):
			make_subplot(fig, result, i)

	plt.subplots_adjust(hspace = 0.4, left = 0.05, right = 0.85, top = 0.80)
	#Save to output path		
	plt.savefig("resultsFile.pdf")


plot_results("", sys.argv[1], sys.argv[2], sys.argv[3])