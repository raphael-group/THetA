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
 # "AS IS" AND DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE
 # AND ANY DOCUMENTATION, INCLUDING ALL IMPLIED WARRANTIES OF
 # MERCHANTABILITY AND FITNESS FOR ANY PARTICULAR PURPOSE. IN NO
 # EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR ANY SPECIAL, INDIRECT OR
 # CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 # LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
 # NEGLIGENCE OR OTHER ACTION BASED ON ANY OTHER LEGAL THEORY,
 # ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
 # SOFTWARE.

 #
 # @author Layla Oesper, Ahmad Mahmoody, Benjamin J. Raphael, Gryte Satas, and Alex Ashery
 ###

from FileIO import *
import matplotlib.pyplot as plt
import numpy as np
import os

def plot_chromosome_clustering(data, mus, sigmas, clusterAssignments, ax):
	data = np.array(data)

	xs = data[:,0]
	ys = data[:,1]

	for i in range(len(mus)):
		currMu = mus[i]
		currSigma = sigmas[i]
		xvals = [xs[j] for j in range(len(xs)) if clusterAssignments[j] == i]
		yvals = [ys[j] for j in range(len(xs)) if clusterAssignments[j] == i]

		plot_gaussian(ax, currMu, currSigma, 'black')
		ax.plot(xvals, yvals, 'o', color='blue', zorder=1)

	ax.set_xlim([0, 5])
	ax.set_ylim([0, 0.5])

def plot_gaussian(ax, mu, Sigma, color):
	"""
	Plots a gaussian on an axis object. This code comes from Michael Hughes' bnpy
	package (see plotGauss2DContour in /bnpy/viz/GaussViz.py)
	"""

	radiusLengths = [0.31863936396437514, 0.67448975019608171, 1.1503493803760079]
	sqrtSigma = np.sqrt(Sigma)

	# Prep for plotting elliptical contours
	# by creating grid of (x,y) points along perfect circle
	ts = np.arange(-np.pi, np.pi, 0.03)
	x = np.sin(ts)
	y = np.cos(ts)
	Zcirc = np.vstack([x, y])

	# Warp circle into ellipse defined by Sigma's eigenvectors
	Zellipse = np.dot(sqrtSigma, Zcirc)

	# plot contour lines across several radius lengths
	for r in radiusLengths:
		Z = r * Zellipse + mu[:, np.newaxis]
		ax.plot(
			Z[0], Z[1], '.', markerfacecolor=color, markeredgecolor=color, zorder=2)

def plot_classifications(mus, sigmas, intervals, clusterAssignments, numClusters, sampleName, singleCopyParamInds, zeroCopyParamInds, ampParamInds, diploidInd, outdir):
	print "Plotting classifications..."

	fig = plt.figure()
	ax = fig.add_subplot(111)

	def color_map(num):
		if num == diploidInd:
			return 'green'
		elif num in singleCopyParamInds:
			return 'red'
		elif num in zeroCopyParamInds:
			return 'orange'
		else:
			return 'blue'

	for i in range(numClusters):
		currMu = mus[i]
		currSigma = sigmas[i]
		currColor = color_map(i)

		currX = np.array([intervals[j][5:7] for j in range(len(intervals)) if clusterAssignments[j] == i])
		if list(currX) == []: continue
		xvals = currX[:,0]
		yvals = currX[:,1]

		plot_gaussian(ax, currMu, currSigma, currColor)
		ax.plot(xvals, yvals, 'o', color=currColor, zorder=1)

	ax.set_title(sampleName + " meta Clustering")
	ax.set_xlim([0, 5])
	ax.set_ylim([0, 0.5])
	fig.savefig(os.path.join(outdir, sampleName + "_classifications.png"))

def plot_clusters(intervals, clusterAssignments, numClusters, sampleName, amp_upper, stepSize, diploidRDR, clonalsingleCopyRDR, outdir):
	print "Plotting clusters..."
	cmap = plt.get_cmap('gist_rainbow')
	colors = [cmap(i) for i in np.linspace(0, 1, numClusters)]

	fig = plt.figure()
	ax = fig.add_subplot(111)
	xs = map(lambda row: row[5], intervals)
	ys = map(lambda row: row[6], intervals)
	colorAssignment = map(lambda assignment: colors[assignment], clusterAssignments)

	ax.scatter(xs, ys, c=colorAssignment)
	ax.plot([clonalsingleCopyRDR, clonalsingleCopyRDR], [0.0, 0.5], color='red')
	ax.plot([diploidRDR, diploidRDR], [0.0, 0.5], color='green')
	if amp_upper != []:
		maxStep = int(max(amp_upper) - 1)
	else:
		maxStep = 1
	for scale in range(1, maxStep):
		barX = (scale * stepSize) + diploidRDR
		ax.plot([barX, barX], [0.0, 0.5], color='blue')
	ax.set_ylim([0, 0.5])
	ax.set_xlim([0, ((maxStep * stepSize) + diploidRDR)])
	fig.savefig(os.path.join(outdir, sampleName + "_assignment.png"))

def parse_preprocessed_data(filename):
	sampleList = []
	numClustersList = []
	clusterLengths = []
	clusterRDRs = []
	clusterMeanBAFs = []
	clusterClassifications = []
	with open(filename) as f:
		f.readline()
		while True:
			sample = f.readline().strip()
			if sample == "": break

			sampleList.append(sample)
			numClusters = int(f.readline().strip())
			numClustersList.append(numClusters)

			lengths = []
			RDRs = []
			meanBAFs = []
			classifications = []
			for i in range(numClusters):
				length, RDR, BAF, classification = f.readline().strip().split("\t")
				lengths.append(int(length))
				RDRs.append(float(RDR))
				meanBAFs.append(float(BAF))
				classifications.append(classification)
			clusterLengths.append(lengths)
			clusterRDRs.append(RDRs)
			clusterMeanBAFs.append(meanBAFs)
			clusterClassifications.append(classifications)

			f.readline()

	return sampleList, numClustersList, clusterLengths, clusterRDRs, clusterMeanBAFs, clusterClassifications

def show_histogram(filename):
	sampleList, numClustersList, clusterLengths, clusterRDRs, clusterMeanBAFs, clusterClassifications = parse_preprocessed_data(filename)

	clusterLengths = [length for sublist in clusterLengths for length in sublist]
	clusterMeanBAFs = [BAF for sublist in clusterMeanBAFs for BAF in sublist]
	clusterRDRs = [RDR for sublist in clusterRDRs for RDR in sublist]

	scaledMeanBAFs = []
	for BAF, length in zip(clusterMeanBAFs, clusterLengths):
		BAFlist = map(lambda x: BAF, range(length / 100000))
		scaledMeanBAFs.append(BAFlist)
	scaledMeanBAFs = [BAF for sublist in scaledMeanBAFs for BAF in sublist]

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.hist(scaledMeanBAFs, bins=35)
	fig.savefig('all_sample_clusters_scaled.png')

def plot_two_largest_from_preprocessed(filename):
	sampleList, numClustersList, clusterLengths, clusterRDRs, clusterMeanBAFs, clusterClassifications = parse_preprocessed_data(filename)

	twoLargestRDRs = []
	twoLargestBAFs = []
	for i in range(len(sampleList)):
		lengths = clusterLengths[i]
		RDRs = clusterRDRs[i]
		meanBAFs = clusterMeanBAFs[i]

		sortedlist = sorted(zip(lengths, RDRs, meanBAFs), key=lambda (a, b, c): a)
		lengths, RDRs, meanBAFs = zip(*sortedlist)
		lengths = list(lengths)[-2:]
		RDRs = list(RDRs)[-2:]
		meanBAFs = list(meanBAFs)[-2:]

		twoLargestRDRs.append(RDRs)
		twoLargestBAFs.append(meanBAFs)

	scatterFig = plt.figure()
	scatterAx = scatterFig.add_subplot(111)
	for RDR, BAF in zip(twoLargestRDRs, twoLargestBAFs):
		if BAF[-1] < 0.1: continue
		scatterAx.plot(RDR, BAF, 'black', zorder=1)
		scatterAx.scatter(RDR, BAF, c=['red', 'blue'], zorder=2)
	scatterFig.savefig('two_largest_mass_clusters.png')

def plot_BAF_by_chrm(intervalfile, resultsfile, clusterAssignments, outdir):
	sampleName = os.path.basename(intervalfile).split(".")[0]

	missingData, intervals = read_interval_RD_BAF_file(intervalfile)
	results = read_results_file_full(resultsfile)
	Carray = results['C']
	muArray = results['mu']
	numResults = len(Carray)

	BAFbyChrm = [[] for i in range(24)]
	for i in range(len(intervals)):
		row = intervals[i]
		row.append(clusterAssignments[i])
		chrm = row[0]
		BAFbyChrm[chrm].append(row)

	fig, Ax = plt.subplots(nrows=numResults, ncols=1)

	for i in range(numResults):
		C = Carray[i]
		mu = muArray[i]
		if numResults > 1:
			ax = Ax[i]
		else:
			ax = Ax

		delta = generate_delta(C, mu)

		for row, deltaj in zip(intervals, delta):
				row[6] = abs(row[6] - deltaj)



		cmap = plt.get_cmap('gist_rainbow')
		colors = [cmap(i) for i in np.linspace(0, 1, max(clusterAssignments) + 1)]

		offset = 0
		xlabelPoints = []
		xlabels = []
		for i in range(24):
			chrmData = BAFbyChrm[i]

			if chrmData == []: continue

			minPos = min(row[1] for row in chrmData)

			for chrm, start, end, tumorCounts, diploidCounts, corrRatio, BAF, numSNPs, assignment in chrmData:
				color = colors[assignment]
				ax.plot([start + offset - minPos, end + offset - minPos], [BAF, BAF], color=color, linewidth=2, solid_capstyle="butt")
			chrmEnd = max([row[2] for row in chrmData])
			labelPoint = (offset + offset + chrmEnd) / 2
			xlabelPoints.append(labelPoint)
			xlabels.append(i)
			ax.plot([offset, offset], [0, 0.5], color='black')
			offset += chrmEnd - minPos


		ax.set_title('BAF for ' + sampleName)
		ax.set_xticks(xlabelPoints)
		ax.set_xticklabels(xlabels)
		ax.set_xlabel('Chromosome')
		ax.set_ylabel('BAF')
		ax.set_xlim([0, offset])
		ax.tick_params(axis='x', labelsize=8)

	fig.tight_layout()
	N = len(muArray[0])
	fig.savefig(outdir + sampleName + "_BAF_by_chrm_N" + str(N) + ".png")

def generate_delta(C, mu):
	def phi(a):
		if a == 0:
			return 0.0
		elif a == 3:
			return 2.0
		else:
			return 1.0

	delta = []
	for row in C:
		if all(item <= 3 and item >= 0 for item in row):
			numerator = sum(map(lambda (a, b): phi(a) * b, zip(row, mu)))
			denominator = sum(map(lambda (a, b): a * b, zip(row, mu)))
			deltaj = (numerator / denominator) - 0.5
			delta.append(deltaj)
		else:
			continue

	return delta

# def plot_vs(filename):
# 	data = parse_preprocessed_data(filename)
# 	for sample, numClusters, lengths, RDRs, meanBAFs, classifications in zip(*data):
# 		diploidInd = classifications.index('DIPLOID')
# 		diploidRDR = RDRs[diploidInd]
# 		diploidBAF = meanBAFs[diploidInd]

# 		leftx = diploidRDR * 0.5
# 		lefty = 0.5

# 		m0 = (diploidBAF - lefty) / (diploidRDR - leftx)
# 		b0 = diploidBAF - (m0 * diploidRDR)
# 		m1 = -(m0**-1)
# 		scores = []
# 		for RDR, BAF, classification in zip(RDRs, meanBAFs, classifications):
# 			if classification != "DIPLOID" and classification != "AMP":
# 				scores.append(float('inf'))
# 				continue
# 			b1 = BAF - (m1 * RDR)
# 			contactx = (b1 - b0) / (m0 - m1)
# 			contacty = (m0 * contactx) + b0
# 			distToContact = euclidean([RDR, BAF], [contactx, contacty])
# 			score = distToContact + log(BAF)
# 			scores.append(score)

# 		diploidInd = np.argmin(scores)
# 		diploidRDR = RDRs[diploidInd]
# 		diploidBAF = meanBAFs[diploidInd]

# 		singleCopyParamInds, zeroCopyParamInds, ampParamInds, diploidInd = classify_clusters_given_diploid(zip(RDRs, meanBAFs), diploidInd)
# 		for i in range(len(classifications)):
# 			if i == diploidInd:
# 				classifications[i] = "DIPLOID"
# 			elif i in singleCopyParamInds:
# 				classifications[i] = "SINGLE"
# 			elif i in zeroCopyParamInds:
# 				classifications[i] = "ZERO"
# 			else:
# 				classifications[i] = "AMP"

# 		leftx = diploidRDR * 0.5
# 		lefty = 0.5

# 		fig = plt.figure()
# 		ax = fig.add_subplot(111)
# 		ax.plot(RDRs, meanBAFs, 'x', color='orange', mew=5, markersize=10, zorder=3)
# 		ax.plot([leftx, diploidRDR], [lefty, diploidBAF], 'black', zorder=1)
# 		ax.set_xlim([0.0, 5.0])
# 		ax.set_ylim([0.0, 0.5])

# 		m0 = (diploidBAF - lefty) / (diploidRDR - leftx)
# 		b0 = diploidBAF - (m0 * diploidRDR)
# 		m1 = -(m0**-1)
# 		scores = []
# 		for RDR, BAF, classification in zip(RDRs, meanBAFs, classifications):
# 			if classification != "SINGLE" and classification != "ZERO":
# 				scores.append(float('inf'))
# 				continue
# 			b1 = BAF - (m1 * RDR)
# 			contactx = (b1 - b0) / (m0 - m1)
# 			contacty = (m0 * contactx) + b0
# 			ax.plot([RDR, contactx], [BAF, contacty], 'r')
# 			distToContact = euclidean([RDR, BAF], [contactx, contacty])
# 			distToIntercept = euclidean([RDR, BAF], [0.0, b0])
# 			score = distToContact + distToIntercept
# 			scores.append(score)

# 		clonalsingleCopyInd = np.argmin(scores)
# 		clonalsingleCopyRDR = RDRs[clonalsingleCopyInd]
# 		clonalsingleCopyBAF = meanBAFs[clonalsingleCopyInd]
# 		ax.plot([clonalsingleCopyRDR], [clonalsingleCopyBAF], 'rx', mew=5, markersize=10, zorder=4)

# 		x = clonalsingleCopyRDR
# 		stepSize = abs(diploidRDR - clonalsingleCopyRDR)
# 		i = 0
# 		while x < 5 and i < 100:
# 			ax.plot([x, x], [0, 0.5], 'black', zorder=2)
# 			x += stepSize
# 			i += 1

# 		intervalfilename = "/gpfs/main/research/compbio/projects/THetA/ICGC-PanCan/processed_data//pilot64/" + sample + "/" + sample + ".gamma.0.2.RD.BAF.intervals.txt"
# 		missingData, intervals = read_interval_RD_BAF_file(intervalfilename, byChrm=False)
# 		mus = map(lambda row: row[5:7], intervals)
# 		numPoints = map(lambda row: row[7] / 10, intervals)
# 		data = generate_data2(mus, numPoints)
# 		xs = [row[0] for row in data]
# 		ys = [row[1] for row in data]
# 		ax.plot(xs, ys, 'o', zorder=1)

# 		fig.savefig(sample + "_v.png")
# 		plt.close('all')
