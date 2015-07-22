#!/usr/bin/python

import os
path = os.path.dirname(os.path.realpath(__file__)) + "/bnpy/results/"
os.environ["BNPYOUTDIR"] = path
import bnpy
from FileIO import read_binned_file
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, ceil
from scipy.spatial.distance import euclidean

def get_data(binned, gene, chrm):
	npArray = np.array(map(lambda x: x[5:7], binned))
	Data = bnpy.data.XData(X=npArray)
	if chrm is None:
		Data.name = gene
	else:
		Data.name = gene + "_chrm_" + str(chrm)
	Data.summary = "Clustering data for " + gene
	if chrm is not None:
		Data.summary += ", chromosome " + str(chrm)

	return Data

def get_normal_ind(mus):
	distances = map(lambda (mux, muy): sqrt(((mux - 1.0) * (mux - 1.0)) + (muy * muy)), mus)
	normal = np.argmin(distances)
	return normal

def plot_gaussian(ax, mu, Sigma, color):

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
    # TODO: instead, choose radius by percentage of prob mass contained within
    for r in radiusLengths:
        Z = r * Zellipse + mu[:, np.newaxis]
        ax.plot(
            Z[0], Z[1], '.', markerfacecolor=color, markeredgecolor=color, zorder=2)

def classify_clusters_old(mus, sigmas):
	normalParamInds = []
	delParamInds = []
	ampParamInds = []
	unknownParamInds = []
	for i in range(len(mus)):
		muX = mus[i][0]
		muY = mus[i][1]
		sigmaX = sigmas[i][0][0]
		sigmaY = sigmas[i][1][1]

		if (muX > 0.9) and (muX < 1.3) and (muY <= 0.2):
			normalParamInds.append(i)
		elif muX < 0.9:
			delParamInds.append(i)
		else:
			ampParamInds.append(i)

	hetDelParamInds = []
	homDelParamInds = []
	if normalParamInds == []:
		hetDelParamInds = delParamInds
	else:
		avgNormMuX = sum(map(lambda x: mus[x][0], normalParamInds)) / len(normalParamInds)
		avgNormMuY = sum(map(lambda x: mus[x][1], normalParamInds)) / len(normalParamInds)

		for i in delParamInds:
			muX = mus[i][0]
			muY = mus[i][1]
			sigmaX = sigmas[i][0][0]
			sigmaY = sigmas[i][1][1]

			if muX < avgNormMuX - 0.2 and muY < avgNormMuY + 0.1:
				homDelParamInds.append(i)
			else:
				hetDelParamInds.append(i)

	return hetDelParamInds, homDelParamInds, ampParamInds, unknownParamInds, normalParamInds #normalInd

def classify_clusters(mus, lengths, clusterAssignments):
	metaLengths = [0 for i in range(len(mus))]
	for i, length in enumerate(lengths):
		metaLengths[clusterAssignments[i]] += length

	meanBAFs = map(lambda x: x[1], mus)
	filteredLengths = map(lambda (BAF, length): -float('inf') if BAF > 0.2 else length, zip(meanBAFs, metaLengths))
	normalInd = np.argmax(filteredLengths)

	normMuX = mus[normalInd][0]
	normMuY = mus[normalInd][1]

	delParamInds = []
	ampParamInds = []
	for i in range(len(mus)):
		if i == normalInd: continue
		
		if mus[i][0] < normMuX:
			delParamInds.append(i)
		else:
			ampParamInds.append(i)

	hetDelParamInds = []
	homDelParamInds = []
	for i in delParamInds:
		muX = mus[i][0]
		muY = mus[i][1]

		if muX < normMuX - 0.2 and muY < normMuY + 0.1:
			homDelParamInds.append(i)
		else:
			hetDelParamInds.append(i)

	return hetDelParamInds, homDelParamInds, ampParamInds, normalInd

def cluster(data, geneName, fig, ax, sf=0.1, chrm=None):
	Data = get_data(data, geneName, chrm)
	
	K = 15
	if Data.X.shape[0] < 15:
		K = Data.X.shape[0]

	hmodel, Info = bnpy.Run.run(Data, 'DPMixtureModel', 'DiagGauss', 'moVB', nLap=100, nTask=1, K=K, moves='birth,merge', targetMaxSize=500, ECovMat='eye', mergeStartLap=10, sF=sf, doWriteStdOut=False)

	observationModel = hmodel.obsModel
	numClusters = observationModel.K

	mus = [observationModel.get_mean_for_comp(k=i) for i in range(numClusters)]
	sigmas = [observationModel.get_covar_mat_for_comp(k=i) for i in range(numClusters)]

	hetDelParamInds, homDelParamInds, ampParamInds, unknownParamInds, normalParamInds = classify_clusters_old(mus, sigmas)
	LP = hmodel.calc_local_params(Data)
	clusterAssignments = np.argmax(LP['resp'], axis=1)

	def color_map(num):
		if num in normalParamInds:
			return 'green'
		elif num in hetDelParamInds:
			return 'red'
		elif num in homDelParamInds:
			return 'orange'
		elif num in ampParamInds:
			return 'blue'
		else:
			return 'black'

	numPoints = []
	for i in range(numClusters):
		currMu = mus[i]
		currSigma = sigmas[i]
		currColor = color_map(i)

		currX = np.array([Data.X[j] for j in range(len(Data.X)) if clusterAssignments[j] == i])
		numPoints.append(currX.shape[0])

		xvals = currX[:,0]
		yvals = currX[:,1]

		plot_gaussian(ax, currMu, currSigma, currColor)
		ax.plot(xvals, yvals, 'o', color=currColor, zorder=1)

	if chrm is not None:
		ax.set_title("Chromosome " + str(chrm))
	else:
		ax.set_title(geneName)
	ax.set_xlim([0, 5])
	ax.set_ylim([0, 0.5])

	return mus, sigmas, numPoints

def meta_cluster(data, gene, binned):
	npArray = np.array(data)
	Data = bnpy.data.XData(X=npArray)
	Data.name = gene + "_meta"
	Data.summary = "Meta clustering."

	hmodel, Info = bnpy.Run.run(Data, 'DPMixtureModel', 'DiagGauss', 'moVB', nLap=200, nTask=1, K=15, moves='birth,merge', targetMaxSize=500, ECovMat='eye', mergeStartLap=10, sF=0.01, doWriteStdOut=False)
	
	observationModel = hmodel.obsModel
	numClusters = observationModel.K

	mus = [observationModel.get_mean_for_comp(k=i) for i in range(numClusters)]
	sigmas = [observationModel.get_covar_mat_for_comp(k=i) for i in range(numClusters)]

	npArray = np.array(map(lambda row: row[5:7], binned))
	Data = bnpy.data.XData(X=npArray)

	LP = hmodel.calc_local_params(Data)
	clusterAssignments = np.argmax(LP['resp'], axis=1)

	return mus, sigmas, clusterAssignments, numClusters

def plot_classifications(mus, sigmas, binned, clusterAssignments, numClusters, gene, hetDelParamInds, homDelParamInds, ampParamInds, normalInd):
	fig = plt.figure()
	ax = fig.add_subplot(111)

	def color_map(num):
		if num == normalInd:
			return 'green'
		elif num in hetDelParamInds:
			return 'red'
		elif num in homDelParamInds:
			return 'orange'
		else:
			return 'blue'

	for i in range(numClusters):
		currMu = mus[i]
		currSigma = sigmas[i]
		currColor = color_map(i)

		currX = np.array([binned[j][5:7] for j in range(len(binned)) if clusterAssignments[j] == i])
		if list(currX) == []: continue
		xvals = currX[:,0]
		yvals = currX[:,1]

		plot_gaussian(ax, currMu, currSigma, currColor)
		ax.plot(xvals, yvals, 'o', color=currColor, zorder=1)

	ax.set_title(gene + " meta Clustering")
	ax.set_xlim([0, 5])
	ax.set_ylim([0, 0.5])
	fig.savefig(gene + "_classifications.png")

def plot_clusters(binned, clusterAssignments, numClusters, geneName, amp_upper, stepSize, normMuX, stepPointX):
	cmap = plt.get_cmap('gist_rainbow')
	colors = [cmap(i) for i in np.linspace(0, 1, numClusters)]

	fig = plt.figure()
	ax = fig.add_subplot(111)
	xs = map(lambda row: row[5], binned)
	ys = map(lambda row: row[6], binned)
	colorAssignment = map(lambda assignment: colors[assignment], clusterAssignments)

	ax.scatter(xs, ys, c=colorAssignment)
	ax.plot([stepPointX, stepPointX], [0.0, 0.5], color='red')
	ax.plot([normMuX, normMuX], [0.0, 0.5], color='green')
	if amp_upper != []:
		maxStep = int(max(amp_upper) - 1)
	else:
		maxStep = 1
	for scale in range(1, maxStep):
		barX = (scale * stepSize) + normMuX
		ax.plot([barX, barX], [0.0, 0.5], color='blue')
	ax.set_ylim([0, 0.5])
	ax.set_xlim([0, ((maxStep * stepSize) + normMuX)])
	fig.savefig(geneName + "_assignment.png")

def generate_data(data, sd=0.02):
	generatedData = []
	for chrm, start, end, tumorCounts, normalCounts, corrRatio, meanBAF, numSNPs in data:
		x = np.random.normal(corrRatio, sd, numSNPs / 10)
		y = np.random.normal(meanBAF, sd, numSNPs / 10)
		newRows = map(lambda (ratio, baf): [chrm, start, end, tumorCounts, normalCounts, ratio, baf, numSNPs], zip(x, y))
		generatedData.append(newRows)

	generatedData = [row for subData in generatedData for row in subData]
	return generatedData

def generate_data2(mus, numPoints, sd=0.05):
	generatedData = []
	for mu, num in zip(mus, numPoints):
		x = np.random.normal(mu[0], sd, num)
		y = np.random.normal(mu[1], sd, num)
		newRows = np.transpose([x, y])
		generatedData.append(newRows)

	generatedData = [row for subData in generatedData for row in subData]
	return generatedData

def group_to_meta_interval(lengths, tumorCounts, normalCounts, m, upper_bounds, lower_bounds, clusterAssignments, numClusters):
	metaLengths = [0 for i in range(numClusters)]
	metaTumorCounts = [0 for i in range(numClusters)]
	metaNormalCounts = [0 for i in range(numClusters)]
	meta_lower_bounds = [2 for i in range(numClusters)]
	meta_upper_bounds = [2 for i in range(numClusters)]
	intervalMap = [[] for i in range(numClusters)]
	for i in range(m):
		if upper_bounds == "X" or lower_bounds == "X": continue

		intervalMap[clusterAssignments[i]].append(i)
		metaLengths[clusterAssignments[i]] += lengths[i]
		metaTumorCounts[clusterAssignments[i]] += tumorCounts[i]
		metaNormalCounts[clusterAssignments[i]] += normalCounts[i]
		meta_lower_bounds[clusterAssignments[i]] = lower_bounds[i]
		meta_upper_bounds[clusterAssignments[i]] = upper_bounds[i]

	return intervalMap, metaLengths, metaTumorCounts, metaNormalCounts, meta_lower_bounds, meta_upper_bounds

def cluster_wrapper((binnedChrm, geneName, fig, currAx, chrm, generateData)):
	if generateData:
		binnedChrm = generate_data(binnedChrm)

	mu, sigmas, numPoints = cluster(binnedChrm, geneName, fig, currAx, chrm=chrm)
	metaDataRow = generate_data2(mu, numPoints)

	return mu, sigmas, numPoints, metaDataRow


def clustering_BAF(filename, byChrm=True, generateData=True):
	geneName = os.path.basename(filename).split(".")[0]
	missingData, binned = read_binned_file(filename, byChrm=byChrm)
	
	if byChrm:
		nrows = 6
		ncols = 4
	else:
		nrows = 1
		ncols = 1
	fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20,20))

	if byChrm:
		metaData = []
		for chrm, binnedChrm in enumerate(binned):
			if binnedChrm == []: continue
			chrm += 1
			print "Clustering " + str(chrm) + "/24"
			currAx = ax[chrm / 4][chrm % 4]

			if generateData:
				binnedChrm = generate_data(binnedChrm)

			mu, sigmas, numPoints = cluster(binnedChrm, geneName, fig, currAx, chrm=chrm)
			metaData.append(generate_data2(mu, numPoints))
		metaData = [row for subData in metaData for row in subData]
	else:
		if generateData:
			binned = generate_data(binned)
		mu, sigmas, numPoints = cluster(binned, geneName, fig, ax)
		metaData = generate_data2(mu, numPoints)

	if byChrm:
		binned = [row for subData in binned for row in subData]

	metaMu, metaSigma, clusterAssignments, numClusters = meta_cluster(metaData, geneName, binned)

	intervalLengths = map(lambda row: row[2] - row[1] + 2, binned)
	hetDelParamInds, homDelParamInds, ampParamInds, normalInd = classify_clusters(metaMu, intervalLengths, clusterAssignments)
	
	plot_classifications(metaMu, metaSigma, binned, clusterAssignments, numClusters, geneName, hetDelParamInds, homDelParamInds, ampParamInds, normalInd)
	fig.savefig(geneName + "_by_chromosome.png")

	normMu = metaMu[normalInd]
	normMuX = normMu[0]

	if hetDelParamInds != []:
		hetBAFs = map(lambda x: metaMu[x][1] if metaMu[x][0] < (normMuX - 0.1) else -float("inf"), hetDelParamInds)
		stepSizeInd = np.argmax(hetBAFs)
		if hetBAFs[stepSizeInd] == -float("inf"):
			hetRDRs = map(lambda x: metaMu[x][0], hetDelParamInds)
			stepSizeInd = np.argmin(hetRDRs)
		stepPoint = metaMu[hetDelParamInds[stepSizeInd]]
		stepPointX = stepPoint[0]
		stepSize = normMuX - stepPoint[0]
	else:
		stepPointX = 0.0
		stepSize = 0.5

	if ampParamInds != []:
		#amplification means
		ampMus = map(lambda x: metaMu[x], ampParamInds)
		ampDistances = map(lambda point: point[0] - normMuX, ampMus)
		amp_upper = map(lambda distance: ceil(distance / stepSize) + 2, ampDistances)
		amp_upper_map = dict(zip(ampParamInds, amp_upper))
	else:
		amp_upper = []
		ampDistances = []

	plot_clusters(binned, clusterAssignments, numClusters, geneName, amp_upper, stepSize, normMuX, stepPointX)

	m = len(binned) + len(missingData)
	lengths = range(m)
	tumorCounts = range(m)
	normalCounts = range(m)
	upper_bounds = range(m)
	lower_bounds = range(m)
	fullClusterAssignments = range(m)
	for row in missingData:
		lengths[row[-1]] = None
		tumorCounts[row[-1]] = None
		normalCounts[row[-1]] = None
		upper_bounds[row[-1]] = None
		lower_bounds[row[-1]] = None
		fullClusterAssignments[row[-1]] = None

	j = 0 #counter for data
	k = 0 #counter for missingData

	for i in range(m):
		if lengths[i] is None:
			row = missingData[k]

			lengths[i] = row[2] - row[1] + 1
			tumorCounts[i] = row[3]
			normalCounts[i] = row[4]
			upper_bounds[i] = "X"
			lower_bounds[i] = "X"
			fullClusterAssignments[i] = -1
			k += 1
		else:
			row = binned[j]
			
			length = row[2] - row[1] + 1
			lengths[i] = length

			tumorCounts[i] = row[3]
			normalCounts[i] = row[4]
			fullClusterAssignments[i] = clusterAssignments[j]

			if clusterAssignments[j] in ampParamInds:
				lower_bounds[i] = 2
				upper_bounds[i] = amp_upper_map[clusterAssignments[j]]
			else:
				upper_bounds[i] = 2
				if clusterAssignments[j] == normalInd:
					lower_bounds[i] = 2
				elif clusterAssignments[j] in hetDelParamInds:
					lower_bounds[i] = 1
				else:
					lower_bounds[i] = 0
			j += 1
	
	return lengths, tumorCounts, normalCounts, m, upper_bounds, lower_bounds, fullClusterAssignments, numClusters, metaMu

def write_clusters_for_all_samples(samplelist):
	f = open("all_sample_clusters.txt", 'w')
	f.write("#length\tmeanRD\tmeanBAF\n")
	for sample in samplelist:
		f.write(sample)
		filename = "/research/compbio/projects/THetA/ICGC-PanCan/processed_data/pilot64/" + sample + "/" + sample + ".gamma.0.2.RD.BAF.intervals.txt"
		try:
			lengths, tumorCounts, normalCounts, m, upper_bounds, lower_bounds, clusterAssignments, numClusters, clusterMeans = clustering_BAF(filename)
			intervalMap, metaLengths, metaTumorCounts, metaNormalCounts, meta_lower_bounds, meta_upper_bounds = group_to_meta_interval(lengths, tumorCounts, normalCounts, m, upper_bounds, lower_bounds, clusterAssignments, numClusters)
			f.write(sample + "\n")
			f.write(str(numClusters) + "\n")
			for i in range(numClusters):
				f.write(str(metaLengths[i]) + "\t" + str(clusterMeans[i][0]) + "\t" + str(clusterMeans[i][1]) + "\n")
			f.write("\n")
		except IOError:
			continue
	f.close()

def show_histogram(filename):
	sampleList = []
	numClustersList = []
	clusterLengths = []
	clusterMeanBAFs = []
	clusterRDRs = []
	with open(filename) as f:
		f.readline()
		while True:
			sample = f.readline().strip()
			if sample == "": break

			sampleList.append(sample)
			numClusters = int(f.readline().strip())
			numClustersList.append(numClusters)
			
			lengths = []
			meanBAFs = []
			RDRs = []
			for i in range(numClusters):
				length, RDR, BAF = f.readline().strip().split("\t")
				lengths.append(int(length))
				RDRs.append(float(RDR))
				meanBAFs.append(float(BAF))
			clusterLengths.append(lengths)
			clusterMeanBAFs.append(meanBAFs)
			clusterRDRs.append(RDRs)

			f.readline()

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

def plot_max_BAF(filename):
	clusterMeanBAFs = []
	clusterRDRs = []
	clusterLengths = []
	with open(filename) as f:
		f.readline()
		while True:
			sample = f.readline().strip()
			if sample == "": break

			numClusters = int(f.readline().strip())

			lengths = []
			meanBAFs = []
			RDRs = []
			for i in range(numClusters):
				length, RDR, BAF = f.readline().strip().split("\t")
				lengths.append(int(length))
				RDRs.append(float(RDR))
				meanBAFs.append(float(BAF))

			sortedlist = sorted(zip(lengths, RDRs, meanBAFs), key=lambda (a, b, c): a)
			lengths, RDRs, meanBAFs = zip(*sortedlist)
			lengths = list(lengths)[-2:]
			RDRs = list(RDRs)[-2:]
			meanBAFs = list(meanBAFs)[-2:]

			clusterLengths.append(lengths)
			clusterMeanBAFs.append(meanBAFs)
			clusterRDRs.append(RDRs)

			f.readline()

	scatterFig = plt.figure()
	scatterAx = scatterFig.add_subplot(111)
	for BAF, RDR in zip(clusterMeanBAFs, clusterRDRs):
		if BAF[-1] < 0.1: continue
		scatterAx.plot(RDR, BAF, 'black', zorder=1)
		scatterAx.scatter(RDR, BAF, c=['red', 'blue'], zorder=2)
	scatterFig.savefig('two_largest_mass_clusters.png')

	# secondRDR = map(lambda x: x[0], clusterRDRs)
	# secondBAF = map(lambda x: x[0], clusterMeanBAFs)
	# secondScatterFig = plt.figure()
	# secondScatterAx = secondScatterFig.add_subplot(111)
	# secondScatterAx.plot(secondRDR, secondBAF, 'o')
	# secondScatterFig.savefig('second_largest_mass_clusters.png')

	# scatterAx.plot(clusterRDRs, clusterMeanBAFs, 'o')
	# scatterFig.savefig('largest_mass_clusters.png')

	# histFig = plt.figure()
	# histAx = histFig.add_subplot(111)
	# histAx.hist(clusterMeanBAFs, bins=30)
	# histFig.savefig('largest_mass_histogram.png')

	# scaledMeanBAFs = []
	# for BAF, length in zip(clusterMeanBAFs, clusterLengths):
	# 	BAFlist = map(lambda x: BAF, range(length / 100000))
	# 	scaledMeanBAFs.append(BAFlist)
	# scaledMeanBAFs = [BAF for sublist in scaledMeanBAFs for BAF in sublist]

	# scaledHistFig = plt.figure()
	# scaledHistAx = scaledHistFig.add_subplot(111)
	# scaledHistAx.hist(scaledMeanBAFs, bins=30)
	# scaledHistFig.savefig('largest_mass_histogram_scaled.png')


if __name__ == "__main__":
	# plot_max_BAF('all_sample_clusters.txt')
	import os
	samplelist = os.listdir("/research/compbio/projects/THetA/ICGC-PanCan/processed_data/pilot64")[:-1]
	# genes = ['0c7af04b-e171-47c4-8be5-5db33f20148e',
	# 		'6847e993-1414-4e6f-a2af-39ebe218dd7c',
	#		'46f19b5c-3eba-4b23-a1ab-9748090ca4e5',
			# '29a00d78-b9bb-4c6b-b142-d5b8bfa63455',
			# '786fc3e4-e2bf-4914-9251-41c800ebb2fa',
	#genes =	['6aa00162-6294-4ce7-b6b7-0c3452e24cd6'] #,
			# '4853fd17-7214-4f0c-984b-1be0346ca4ab']
	#write_clusters_for_all_samples(samplelist)
	for gene in samplelist:
		clustering_BAF("/gpfs/main/research/compbio/projects/THetA/ICGC-PanCan/processed_data//pilot64/" + gene + "/" + gene + ".gamma.0.2.RD.BAF.intervals.txt")
	#print group_to_meta_interval(*clustering_BAF("/research/compbio/projects/THetA/ICGC-PanCan/processed_data/pilot64/" + gene + "/" + gene + ".gamma.0.2.RD.BAF.intervals.txt"))