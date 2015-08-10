#!/usr/bin/python

import os
os.environ["BNPYOUTDIR"] = "./"
import bnpy
from FileIO import *
from ClusterPlottingTools import *
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, ceil, log
from scipy.spatial.distance import euclidean
from multiprocessing import Pool

def clustering_BAF(filename, byChrm=True, generateData=True, prefix=None, outdir="./", numProcesses=1):
	if not outdir.endswith("/"):
		outdir += "/"

	geneName = os.path.basename(filename).split(".")[0]
	if prefix is None:
		prefix = geneName

	os.environ["BNPYOUTDIR"] = outdir + prefix + "_cluster_data/"

	missingData, intervals = read_interval_RD_BAF_file(filename, byChrm=byChrm)

	metaData = generate_meta_data(intervals, byChrm, numProcesses, geneName, generateData, outdir)

	if byChrm:
		intervals = [row for subData in intervals for row in subData]

	print "Begin meta clustering..."
	metaMu, metaSigma, clusterAssignments, numPoints, numClusters = cluster(metaData, geneName, sf=0.01, intervals=intervals)

	print "Classifying clusters..."
	intervalLengths = [row[2] - row[1] + 2 for row in intervals]
	hetDelParamInds, clonalHetDelParamInd, homDelParamInds, ampParamInds, normalInd = classify_clusters(metaMu, intervalLengths, clusterAssignments)
	
	print "Plotting classifications..."
	plot_classifications(metaMu, metaSigma, intervals, clusterAssignments, numClusters, geneName, hetDelParamInds, homDelParamInds, ampParamInds, normalInd, outdir)

	print "Determining copy number bounds..."
	normMu = metaMu[normalInd]
	normMuX = normMu[0]

	if hetDelParamInds != []:
		stepPointX = metaMu[clonalHetDelParamInd][0]
		stepSize = normMuX - stepPointX
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

	m = len(intervals) + len(missingData)
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
			row = intervals[j]
			
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
	
	print "Plotting clusters..."
	plot_clusters(intervals, clusterAssignments, numClusters, geneName, amp_upper, stepSize, normMuX, stepPointX, outdir)

	return lengths, tumorCounts, normalCounts, m, upper_bounds, lower_bounds, fullClusterAssignments, numClusters, metaMu, normalInd

def generate_meta_data(intervals, byChrm, numProcesses, geneName, generateData, outdir):
	if byChrm:
		print "First round of clustering..."
		metaData = []
		p = Pool(numProcesses)
		
		linearizedChrm = range(24)
		linearizedGeneName = [geneName for i in linearizedChrm]
		linearizedGenerateData = [generateData for i in linearizedChrm]

		results = p.map(cluster_wrapper, zip(intervals, linearizedGeneName, linearizedChrm, linearizedGenerateData))
		
		fig, ax = plt.subplots(nrows=6, ncols=4, figsize=(20,20))
		metaData = []
		for i in linearizedChrm:
			row = results[i]
			if row is None: continue

			generatedData, mus, sigmas, clusterAssignments, metaDataRow = row
			currAx = ax[i / 4][i % 4]
			plot_chromosome_clustering(generatedData, mus, sigmas, clusterAssignments, currAx)
			metaData += metaDataRow

		fig.savefig(outdir + geneName + "_by_chromosome.png")

	else:
		metaData = map(lambda row: row[5:7], intervals)
		if generateData:
			numPoints = map(lambda row: (row[2] - row[1] + 1) / 100000, intervals)
			metaData = generate_data(metaData, numPoints)

	return metaData

def cluster_wrapper((binnedChrm, geneName, chrm, generateData)):
	if binnedChrm == []:
		return None

	if generateData:
		means = [row[5:7] for row in binnedChrm]
		numPoints = [row[7] / 10 for row in binnedChrm]
		binnedChrm = generate_data(means, numPoints)
	else:
		binnedChrm = [row[5:7] for row in binnedChrm]

	mus, sigmas, clusterAssignments, numPoints, numClusters = cluster(binnedChrm, geneName, chrm=chrm)
	metaDataRow = generate_data(mus, numPoints)

	return binnedChrm, mus, sigmas, clusterAssignments, metaDataRow

def generate_data(mus, numPoints, sdx=0.1, sdy=0.02):
	generatedData = []
	for mu, num in zip(mus, numPoints):
		np.random.seed(seed=0)
		x = np.random.normal(mu[0], sdx, num)
		y = np.random.normal(mu[1], sdy, num)
		newRows = np.transpose([x, y])
		generatedData.append(newRows)

	generatedData = [row for subData in generatedData for row in subData]
	return generatedData

def cluster(data, geneName, sf=0.1, chrm=-1, intervals=None):
	Data = format_data(data, geneName, chrm)

	K = 15
	if Data.X.shape[0] < 15:
		K = Data.X.shape[0]

	hmodel, Info = bnpy.Run.run(Data, 'DPMixtureModel', 'DiagGauss', 'moVB', nLap=100, nTask=1, K=K, moves='birth,merge', targetMaxSize=500, ECovMat='eye', mergeStartLap=10, sF=sf, doWriteStdOut=False)

	observationModel = hmodel.obsModel
	numClusters = observationModel.K

	mus = [observationModel.get_mean_for_comp(k=i) for i in range(numClusters)]
	sigmas = [observationModel.get_covar_mat_for_comp(k=i) for i in range(numClusters)]

	if intervals is not None:
		points = [row[5:7] for row in intervals]
		Data = format_data(points, geneName, -1)

	LP = hmodel.calc_local_params(Data)
	clusterAssignments = np.argmax(LP['resp'], axis=1)

	numPoints = []
	for i in range(numClusters):
		currX = np.array([Data.X[j] for j in range(len(Data.X)) if clusterAssignments[j] == i])
		numPoints.append(currX.shape[0])

	return mus, sigmas, clusterAssignments, numPoints, numClusters

def format_data(data, gene, chrm):
	npArray = np.array(data)
	Data = bnpy.data.XData(X=npArray)
	if chrm == -1:
		Data.name = gene + "_meta"
		Data.summary = "Meta clustering for " + gene
	else:
		Data.name = gene + "_chrm_" + str(chrm)
		Data.summary = "Clustering data for " + gene + ", chromosome " + str(chrm)

	return Data

def classify_clusters(mus, lengths, clusterAssignments):
	metaLengths = [0 for i in range(len(mus))]
	for i, length in enumerate(lengths):
		if length is None: continue
		metaLengths[clusterAssignments[i]] += length

	meanBAFs = map(lambda x: x[1], mus)
	filteredLengths = map(lambda (BAF, length): -float('inf') if BAF > 0.2 else length, zip(meanBAFs, metaLengths))
	normalInd = np.argmax(filteredLengths)

	hetDelParamInds, homDelParamInds, ampParamInds, normalInd = classify_clusters_given_normal(mus, normalInd)

	normalRDR = mus[normalInd][0]
	normalBAF = mus[normalInd][1]

	leftx = normalRDR * 0.5
	lefty = 0.5

	m0 = (normalBAF - lefty) / (normalRDR - leftx)
	b0 = normalBAF - (m0 * normalRDR)
	m1 = -(m0**-1)
	scores = []
	for i in range(len(mus)):
		if i != normalInd and i not in ampParamInds:
			scores.append(float('inf'))
			continue

		RDR = mus[i][0]
		BAF = mus[i][1]
		b1 = BAF - (m1 * RDR)
		contactx = (b1 - b0) / (m0 - m1)
		contacty = (m0 * contactx) + b0
		distToContact = euclidean([RDR, BAF], [contactx, contacty])
		score = distToContact + log(BAF)
		scores.append(score)

	normalInd = np.argmin(scores)

	hetDelParamInds, homDelParamInds, ampParamInds, normalInd = classify_clusters_given_normal(mus, normalInd)

	normalRDR = mus[normalInd][0]
	normalBAF = mus[normalInd][1]
	leftx = normalRDR * 0.5
	lefty = 0.5

	m0 = (normalBAF - lefty) / (normalRDR - leftx)
	b0 = normalBAF - (m0 * normalRDR)
	m1 = -(m0**-1)
	scores = []
	for i in range(len(mus)):
		RDR = mus[i][0]
		BAF = mus[i][1]

		if i not in hetDelParamInds and i not in homDelParamInds:
			scores.append(float('inf'))
			continue
		b1 = BAF - (m1 * RDR)
		contactx = (b1 - b0) / (m0 - m1)
		contacty = (m0 * contactx) + b0
		distToContact = euclidean([RDR, BAF], [contactx, contacty])
		distToIntercept = euclidean([RDR, BAF], [0.0, b0])
		score = distToContact + distToIntercept
		scores.append(score)

	clonalHetDelParamInd = np.argmin(scores)
	
	return hetDelParamInds, clonalHetDelParamInd, homDelParamInds, ampParamInds, normalInd

def classify_clusters_given_normal(mus, normalInd):
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

def write_clusters_for_all_samples(samplelist):
	f = open("all_sample_clusters.txt", 'w')
	f.write("#length\tmeanRD\tmeanBAF\tclassification\n")
	for sample in samplelist:
		filename = "/research/compbio/projects/THetA/ICGC-PanCan/processed_data/pilot64/" + sample + "/" + sample + ".gamma.0.2.RD.BAF.intervals.txt"
		try:
			lengths, tumorCounts, normalCounts, m, upper_bounds, lower_bounds, clusterAssignments, numClusters, clusterMeans, normalInd = clustering_BAF(filename)
			hetDelParamInds, homDelParamInds, ampParamInds, normalInd = classify_clusters(clusterMeans, lengths, clusterAssignments)
			intervalMap, metaLengths, metaTumorCounts, metaNormalCounts, meta_lower_bounds, meta_upper_bounds = group_to_meta_interval(lengths, tumorCounts, normalCounts, m, upper_bounds, lower_bounds, clusterAssignments, numClusters)
			f.write(sample + "\n")
			f.write(str(numClusters) + "\n")
			for i in range(numClusters):
				f.write(str(metaLengths[i]) + "\t" + str(clusterMeans[i][0]) + "\t" + str(clusterMeans[i][1]) + "\t")
				if i == normalInd:
					f.write("NORMAL")
				elif i in hetDelParamInds:
					f.write("HETDEL")
				elif i in homDelParamInds:
					f.write("HOMDEL")
				else:
					f.write("AMP")
				f.write("\n")
			f.write("\n")
			f.flush()
		except IOError:
			continue
	f.close()