#!/usr/bin/python

import os
path = os.path.dirname(os.path.realpath(__file__)) + "/bnpy/results/"
os.environ["BNPYOUTDIR"] = path
import bnpy
from FileIO import read_binned_file
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import pylab
from math import sqrt, ceil

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

#====NOTE: this was the old plot_gaussian I was using. current one is stolen from mike's code=====#
# def plot_gaussian(ax, mu, sigma, color):
# 	muX, muY = mu
# 	sigmaX = sigma[0][0]
# 	sigmaY = sigma[1][1]

# 	xmin = muX - (5 * sigmaX)
# 	xmax = muX + (5 * sigmaX)
# 	ymin = muY - (5 * sigmaY)
# 	ymax = muY + (5 * sigmaY)
	
# 	delta = 0.001
# 	x = np.arange(xmin, xmax, delta)
# 	y = np.arange(ymin, ymax, delta)
# 	X, Y = np.meshgrid(x, y)
# 	Z = mlab.bivariate_normal(X, Y, sigmaX, sigmaY, muX, muY)
# 	ax.contour(X, Y, Z, colors=color, zorder=100)

def classify_clusters(mus, sigmas):
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

def cluster(data, geneName, fig, ax, sf=0.1, chrm=None):
	Data = get_data(data, geneName, chrm)
	
	K = 15
	if Data.X.shape[0] < 15:
		K = Data.X.shape[0]

	hmodel, Info = bnpy.Run.run(Data, 'DPMixtureModel', 'DiagGauss', 'moVB', nLap=100, nTask=1, K=K, moves='birth,merge', targetMaxSize=500, ECovMat='eye', mergeStartLap=10, sF=sf, doWriteStdOut=False)

	observationModel = hmodel.obsModel
	numClusters = observationModel.K

	#xvals = Data.X[:,0]
	#yvals = Data.X[:,1]
	mus = [observationModel.get_mean_for_comp(k=i) for i in range(numClusters)]
	sigmas = [observationModel.get_covar_mat_for_comp(k=i) for i in range(numClusters)]

	#===NOTE: necessary for bnpy built-in plotting only===#
	# if chrm is not None:
	# 	resultsPath = path + geneName + "_chrm_" + str(chrm) + "/defaultjob/"
	# else:
	# 	resultsPath = path + geneName + "/defaultjob/"

	hetDelParamInds, homDelParamInds, ampParamInds, unknownParamInds, normalParamInds = classify_clusters(mus, sigmas)
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

	#bnpy built-in plotting
	# bnpy.viz.PlotComps.plotCompsForJob(resultsPath, figH=fig)

	if chrm is not None:
		ax.set_title("Chromosome " + str(chrm))
	else:
		ax.set_title(geneName)
	ax.set_xlim([0, 5])
	ax.set_ylim([0, 0.5])

	return mus, sigmas, numPoints

def plot_all_chrms(fullDataX, fullDataY, gene):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	for i in range(24):
		if fullDataX[i] is None and fullDataY[i] is None: continue
		ax.plot(fullDataX[i], fullDataY[i], 'o', zorder=1)
		chrm = i + 1
		resultsPath = path + gene + "_chrm_" + str(chrm) + "/defaultjob/"
		bnpy.viz.PlotComps.plotCompsForJob(resultsPath, figH=fig)
	plt.show()

def plot_all_means(muX, muY):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(muX, muY, 'o')
	plt.show()

def meta_cluster(data, gene):
	npArray = np.array(data)
	Data = bnpy.data.XData(X=npArray)
	Data.name = gene + "_meta"
	Data.summary = "Meta clustering."

	hmodel, Info = bnpy.Run.run(Data, 'DPMixtureModel', 'DiagGauss', 'moVB', nLap=200, nTask=1, K=15, moves='birth,merge', targetMaxSize=500, ECovMat='eye', mergeStartLap=10, sF=0.01, doWriteStdOut=False)
	
	observationModel = hmodel.obsModel
	numClusters = observationModel.K

	mus = [observationModel.get_mean_for_comp(k=i) for i in range(numClusters)]
	sigmas = [observationModel.get_covar_mat_for_comp(k=i) for i in range(numClusters)]

	hetDelParamInds, homDelParamInds, ampParamInds, unknownParamInds, normalParamInds = classify_clusters(mus, sigmas)
	LP = hmodel.calc_local_params(Data)
	clusterAssignments = np.argmax(LP['resp'], axis=1)

	fig = plt.figure()
	ax = fig.add_subplot(111)

	return mus, sigmas, hmodel

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

	ax.plot(Data.X[:,0], Data.X[:,1], 'o', zorder=1)
	ax.set_title(gene + " meta Clustering")
	ax.set_xlim([0, 5])
	ax.set_ylim([0, 0.5])
	fig.savefig(gene + "_meta.png")

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
	intervalMap = {}
	
	#layla - change to dictionary to handle cluster of -1
	for val in range(numClusters):
		intervalMap[val] = []
	intervalMap[-1] = []
	#intervalMap = [[] for i in range(numClusters)]
	for i in range(m):
		if upper_bounds[i] == "X" or lower_bounds[i] == "X" or clusterAssignments[i] == -1: 
			intervalMap[clusterAssignments[i]].append(i)
			continue
		else:
			intervalMap[clusterAssignments[i]].append(i)
			metaLengths[clusterAssignments[i]] += lengths[i]
			metaTumorCounts[clusterAssignments[i]] += tumorCounts[i]
			metaNormalCounts[clusterAssignments[i]] += normalCounts[i]
			meta_lower_bounds[clusterAssignments[i]] = lower_bounds[i]
			meta_upper_bounds[clusterAssignments[i]] = upper_bounds[i]

	return intervalMap, metaLengths, metaTumorCounts, metaNormalCounts, meta_lower_bounds, meta_upper_bounds

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
		fullDataX = []
		fullDataY = []
		muX = []
		muY = []
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

	metaMu, metaSigma, metaHModel = meta_cluster(metaData, geneName)

	if byChrm:
		binned = [row for subData in binned for row in subData]

	npArray = np.array(map(lambda row: row[5:7], binned))
	Data = bnpy.data.XData(X=npArray)

	LP = metaHModel.calc_local_params(Data)
	clusterAssignments = np.argmax(LP['resp'], axis=1)

	cmap = plt.get_cmap('gist_rainbow')
	numClusters = len(metaMu)
	colors = [cmap(i) for i in np.linspace(0, 1, numClusters)]

	clusterFig = plt.figure()
	clusterAx = clusterFig.add_subplot(111)
	xs = npArray[:,0]
	ys = npArray[:,1]
	colorAssignment = map(lambda assignment: colors[assignment], clusterAssignments)

	clusterAx.scatter(xs, ys, c=colorAssignment)
	clusterFig.savefig(geneName + "_meta_assignment.png")
	fig.savefig(geneName + "_clusters.png")

	hetDelParamInds, homDelParamInds, ampParamInds, unknownParamInds, normalParamInds = classify_clusters(metaMu, metaSigma)

	if ampParamInds != []:
		if normalParamInds == []:
			avgNormMuX = 1.0
			avgNormMuY = 0.0
		else:
			avgNormMuX = sum(map(lambda x: metaMu[x][0], normalParamInds)) / len(normalParamInds)
			avgNormMuY = sum(map(lambda x: metaMu[x][1], normalParamInds)) / len(normalParamInds)

		if hetDelParamInds != []:
			avgHetMuX = sum(map(lambda x: metaMu[x][0], hetDelParamInds)) / len(hetDelParamInds)
			avgHetMuY = sum(map(lambda x: metaMu[x][1], hetDelParamInds)) / len(hetDelParamInds)
		else:
			avgHetMuX = 0.0
			avgHetMuY = 1.0

		length = sqrt(((avgHetMuX - avgNormMuX)**2) + ((avgHetMuY - avgNormMuY)**2))

		ampMuXs = map(lambda x: metaMu[x][0], ampParamInds)
		ampMuYs = map(lambda x: metaMu[x][1], ampParamInds)

		squareDiffMuX = map(lambda x: (x - avgNormMuX)**2, ampMuXs)
		squareDiffMuY = map(lambda x: (x - avgNormMuY)**2, ampMuYs)

		distances = map(lambda (dX, dY): sqrt(dX + dY), zip(squareDiffMuX, squareDiffMuY))
		maxInd = np.argmax(distances)

		maxAmpMuX = metaMu[ampParamInds[maxInd]][0]
		maxAmpMuY = metaMu[ampParamInds[maxInd]][1]

		ampDist = sqrt((maxAmpMuX - avgNormMuX)**2 + (maxAmpMuY - avgNormMuY)**2)
		amp_upper = int(ceil(ampDist / length) + 2)
	else:
		amp_upper = 2

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
				upper_bounds[i] = amp_upper
			else:
				upper_bounds[i] = 2
				if clusterAssignments[j] in normalParamInds:
					lower_bounds[i] = 2
				elif clusterAssignments[j] in hetDelParamInds:
					lower_bounds[i] = 1
				else:
					lower_bounds[i] = 0
			j += 1
	
	return lengths, tumorCounts, normalCounts, m, upper_bounds, lower_bounds, fullClusterAssignments, numClusters, metaMu

#experiment in clustering across interval files
# def plot_intervals(genes):
# 	X = []
# 	Y = []
# 	for gene in genes:
# 		filename = "/research/compbio/projects/THetA/ICGC-PanCan/processed_data/pilot64/" + gene + "/" + gene + ".gamma.0.2.RD.BAF.intervals.txt"
# 		binned = read_binned_file(filename, byChrm=False)
# 		binned = generate_data(binned)
# 		x = map(lambda row: row[5], binned)
# 		y = map(lambda row: row[6], binned)
# 		X.append(x)
# 		Y.append(y)
	
# 	X = [element for sublist in X for element in sublist]
# 	Y = [element for sublist in Y for element in sublist]
# 	nparray = np.array([X, Y])
# 	Data = bnpy.data.XData(X=np.transpose(nparray))
# 	Data.name = "full"
# 	Data.summary = "Clustering data for full set" + gene

# 	hmodel, Info = bnpy.Run.run(Data, 'DPMixtureModel', 'DiagGauss', 'moVB', nLap=250, nTask=1, K=15, moves='birth,merge', targetMaxSize=500, ECovMat='eye', mergeStartLap=10, sF=0.1)
# 	resultsPath = path + "full" + "/defaultjob/"
# 	bnpy.viz.PlotComps.plotCompsForJob(resultsPath)
# 	plt.xlim([0, 5])
# 	plt.ylim([0, 0.5])
# 	plt.plot(X, Y, 'o', zorder=1)
# 	plt.show()

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


if __name__ == "__main__":
	import os
	samplelist = os.listdir("/research/compbio/projects/THetA/ICGC-PanCan/processed_data/pilot64")[:-1]
	# genes = ['0c7af04b-e171-47c4-8be5-5db33f20148e',
	# 			'6847e993-1414-4e6f-a2af-39ebe218dd7c',
	# 			'46f19b5c-3eba-4b23-a1ab-9748090ca4e5',
	# 			'29a00d78-b9bb-4c6b-b142-d5b8bfa63455',
	# 			'786fc3e4-e2bf-4914-9251-41c800ebb2fa',
	# 			'6aa00162-6294-4ce7-b6b7-0c3452e24cd6',
	# 			'4853fd17-7214-4f0c-984b-1be0346ca4ab']
	write_clusters_for_all_samples(samplelist)
	# for gene in genes:
	#clustering_BAF("/gpfs/main/research/compbio/projects/THetA/ICGC-PanCan/data/" + gene + "/" + gene + ".gamma.0.2.RD.BAF.intervals.txt")
	#print group_to_meta_interval(*clustering_BAF("/research/compbio/projects/THetA/ICGC-PanCan/processed_data/pilot64/" + gene + "/" + gene + ".gamma.0.2.RD.BAF.intervals.txt"))