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
	"""
	Performs clustering on interval data and analyzes clusters to determine copy
	number bounds and which intervals belong to similar clusters.

	Arguments:
		filename (string): Location of the input interval file.
		byChrm (bool): Indicates if intra-chromosome interval clustering should be
						performed before full genome clustering. True to do intra-chromosome
						interval clustering, False otherwise.
		generateData (bool): Indicates if clustering should be performed on random data generated
								around intervals, where the number of data points is proportional
								to the interval length. True to generate data, False otherwise.
		prefix (string): Prefix used for output files.
		outdir (string): Target directory for output files.
		numProcesses (int): Number of processors used (only applies for intra-chromosome clustering).

	Returns:
		lengths (list of ints): Length of each interval
		tumorCounts (list of ints): Tumor count for each interval
		normalCounts (list of ints): Normal count for each interval
		m (int): The total number of intervals (including intervals that have been marked as invalid).
		upper_bounds (list of ints): The upper copy number bound assigned to each interval.
		lower_bounds (list of ints): The lower copy number bound assigned to each interval.
		clusterAssignments (list of ints): The assignment of each interval to a cluster, where an entry
											j at index i means the ith interval has been assigned to the
											jth meta-interval.
		numClusters (int): The number of clusters.
		metaMu (list of lists of floats): The means of each cluster.
		normalInd (int): The index of the cluster that has been called normal.
	"""

	#ensure / at the end of output directory
	if not outdir.endswith("/"):
		outdir += "/"


	sampleName = os.path.basename(filename).split(".")[0]
	if prefix is None:
		prefix = sampleName

	#setting bnpy outdir
	os.environ["BNPYOUTDIR"] = outdir + prefix + "_cluster_data/"

	missingData, intervals = read_interval_RD_BAF_file(filename, byChrm=byChrm)

	metaData = generate_meta_data(intervals, byChrm, numProcesses, sampleName, generateData, outdir)

	#flatten interval data if it has been grouped by chromosome
	if byChrm:
		intervals = [row for subData in intervals for row in subData]

	print "Begin meta clustering..."
	metaMu, metaSigma, clusterAssignments, numPoints, numClusters = cluster(metaData, sampleName, sf=0.01, intervals=intervals)

	intervalLengths = [row[2] - row[1] + 2 for row in intervals]
	hetDelParamInds, clonalHetDelParamInd, homDelParamInds, ampParamInds, normalInd = classify_clusters(metaMu, intervalLengths, clusterAssignments)
	
	plot_classifications(metaMu, metaSigma, intervals, clusterAssignments, numClusters, sampleName, hetDelParamInds, homDelParamInds, ampParamInds, normalInd, outdir)

	lengths, tumorCounts, normalCounts, upper_bounds, lower_bounds, clusterAssignments, m = process_classifications(intervals, missingData, metaMu, clusterAssignments, numClusters, normalInd, clonalHetDelParamInd, hetDelParamInds, ampParamInds, sampleName, outdir)

	return lengths, tumorCounts, normalCounts, m, upper_bounds, lower_bounds, clusterAssignments, numClusters, metaMu, normalInd

def generate_meta_data(intervals, byChrm, numProcesses, sampleName, generateData, outdir):
	"""
	Generates and/or packages data to be clustered from interval data.

	Arguments:
		intervals (list of lists): 2D list containing data related to each interval.
		byChrm (bool): Indicates if intra-chromosome interval clustering should be
						performed before full genome clustering. True to do intra-chromosome
						interval clustering, False otherwise.
		numProcesses (int): Number of processors to be used.
		sampleName (string): The name of the input sample.
		generateData (bool): Indicates if clustering should be performed on random data generated
								around intervals, where the number of data points is proportional
								to the interval length. True to generate data, False otherwise.
		outdir (string): Target directory for output files.

	Returns:
		metaData (list of lists of floats): Data points to be clustered.
	"""

	if byChrm:
		print "First round of clustering..."
		metaData = []
		p = Pool(numProcesses)
		
		linearizedChrm = range(24)
		linearizedSampleName = [sampleName for i in linearizedChrm]
		linearizedGenerateData = [generateData for i in linearizedChrm]

		results = p.map(cluster_wrapper, zip(intervals, linearizedSampleName, linearizedChrm, linearizedGenerateData))
		
		fig, ax = plt.subplots(nrows=6, ncols=4, figsize=(20,20))
		metaData = []
		for i in linearizedChrm:
			row = results[i]
			if row is None: continue

			generatedData, mus, sigmas, clusterAssignments, metaDataRow = row
			currAx = ax[i / 4][i % 4]
			plot_chromosome_clustering(generatedData, mus, sigmas, clusterAssignments, currAx)
			metaData += metaDataRow

		fig.savefig(outdir + sampleName + "_by_chromosome.png")
	else:
		metaData = [row[5:7] for row in intervals]
		if generateData:
			numPoints = [(row[2] - row[1] + 1) / 100000 for row in intervals]
			metaData = generate_data(metaData, numPoints)

	return metaData

def cluster_wrapper((binnedChrm, sampleName, chrm, generateData)):
	"""
	Wrapper function used for intra-chromosome clustering.

	Arguments:
		binnedChrm (list of lists): 2D list containing data related to each interval within a chromosome.
		sampleName (string): The name of the input sample.
		chrm (int): Number corresponding to the chromosome that is being clustered.
		generateData (bool): Indicates if clustering should be performed on random data generated
								around intervals, where the number of data points is proportional
								to the interval length. True to generate data, False otherwise.
	
	Returns:
		binnedChrm (list of lists of floats): 2D list of points that have been clustered.
		mus (list of lists of floats): Means of clusters.
		sigmas (list of 2D lists of floats): list of 2D lists, where each 2D list represents the
												covariance matrix of a cluster.
		clusterAssignments (list of ints): The assignment of each interval to a cluster, where an entry
											j at index i means the ith interval has been assigned to the
											jth meta-interval.
		metaDataRow (list of lists of floats): 2D list of points to be clustered in full-genome clustering
	"""

	if binnedChrm == []:
		return None

	if generateData:
		means = [row[5:7] for row in binnedChrm]
		numPoints = [row[7] / 10 for row in binnedChrm]
		binnedChrm = generate_data(means, numPoints)
	else:
		binnedChrm = [row[5:7] for row in binnedChrm]

	mus, sigmas, clusterAssignments, numPoints, numClusters = cluster(binnedChrm, sampleName, chrm=chrm)
	metaDataRow = generate_data(mus, numPoints)

	return binnedChrm, mus, sigmas, clusterAssignments, metaDataRow

def generate_data(mus, numPoints, sdx=0.1, sdy=0.02):
	"""
	Randomly generates normally distributed data points about several means.

	Arguments:
		mus (list of list of floats): Means about which data points are drawn.
		numPoints (list of ints): Number of points to draw from each distribution.
		sdx (float): x standard deviation for each distribution.
		sdy (float): y standard deviation for each distribution

	Returns:
		generatedData (list of lists of floats): Data points generated from the distributions.
	"""

	generatedData = []
	for mu, num in zip(mus, numPoints):
		np.random.seed(seed=0)
		x = np.random.normal(mu[0], sdx, num)
		y = np.random.normal(mu[1], sdy, num)
		newRows = np.transpose([x, y])
		generatedData.append(newRows)

	generatedData = [row for subData in generatedData for row in subData]
	return generatedData

def cluster(data, sampleName, sf=0.1, chrm=None, intervals=None):
	"""
	Clusters a set of data points lying in an arbitrary number of clusters.

	Arguments:
		data (list of lists of floats): list of data points to be clustered.
		sampleName (string): The name of the input sample.
		sf (float): Tuning parameter for clustering; used to determine initial size of
					distribution covariances. Small sf indicates a belief that clusters
					are of small size.
		chrm (int): Number corresponding to the chromosome that is being clustered; None if 
					clustering is not intra-chromosome.
		intervals (list of lists of floats): List of original interval points. Set to None
											if generated data is given clusterAssignments
											instead of original intervals.
	
	Returns:
		mus (list of lists of floats): List of cluster means.
		sigmas (list of 2D lists of floats): List of cluster covariances.
		clusterAssignments (list of ints): The assignment of each interval to a cluster, where an entry
											j at index i means the ith interval has been assigned to the
											jth meta-interval.
		numPoints (list of ints): Number of points assigned to each cluster
		numClusters (int): The number of clusters.
	"""

	Data = format_data(data, sampleName, chrm)

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
		Data = format_data(points, sampleName, -1)

	LP = hmodel.calc_local_params(Data)
	clusterAssignments = np.argmax(LP['resp'], axis=1)

	numPoints = []
	for i in range(numClusters):
		currX = np.array([Data.X[j] for j in range(len(Data.X)) if clusterAssignments[j] == i])
		numPoints.append(currX.shape[0])

	return mus, sigmas, clusterAssignments, numPoints, numClusters

def format_data(data, sampleName, chrm):
	"""
	Creates XData object to be clustered by bnpy.

	Arguments:
		data (list of lists of floats): List of data points to be formatted for clustering.
		sampleName (string): The name of the input sample.
		chrm (int): Number corresponding to the chromosome that is being formatted; None if 
					clustering is not intra-chromosome.

	Returns:
		Data (XData object): input object for bnpy clustering.
	"""

	npArray = np.array(data)
	Data = bnpy.data.XData(X=npArray)
	if chrm is None:
		Data.name = sampleName + "_meta"
		Data.summary = "Meta clustering for " + sampleName
	else:
		Data.name = sampleName + "_chrm_" + str(chrm)
		Data.summary = "Clustering data for " + sampleName + ", chromosome " + str(chrm)

	return Data

def classify_clusters(mus, lengths, clusterAssignments):
	"""
	Assigns each cluster created by clustering a classification.

	Arguments:
		mus (list of lists of floats): List of cluster means
		lengths (list of ints): Length of each input interval.
		clusterAssignments (list of ints): The assignment of each interval to a cluster, where an entry
											j at index i means the ith interval has been assigned to the
											jth meta-interval.

	Returns:
		hetDelParamInds (list of ints): List of indices of clusters that have been classified as 
										heterozygous deletions.
		clonalHetDelParamInd (int): Index of the cluster that has been classified as the clonal 
									heterozygous deletion.
		homDelParamInds (list of ints): List of indices of clusters that have been classified as 
										homozygous deletions.
		ampParamInds (list of ints): List of indices of clusters that have been classified as
									 amplifications.
		normalInd (int): Index of the cluster that has been classified as the normal cluster.
	"""

	print "Classifying clusters..."

	#generate lengths for meta-intervals
	metaLengths = [0 for i in range(len(mus))]
	for length, assignment in zip(lengths, clusterAssignments):
		if length is not None: metaLengths[assignment] += length

	#guess the normal cluster to be the largest cluster that has BAF less than 0.2
	meanBAFs = [x[1] for x in mus]
	filteredLengths = map(lambda (BAF, length): -float('inf') if BAF > 0.2 else length, zip(meanBAFs, metaLengths))
	normalInd = np.argmax(filteredLengths)

	#classify clusters given naive guess
	hetDelParamInds, homDelParamInds, ampParamInds = classify_clusters_given_normal(mus, normalInd)
	#revise decision of normal cluster
	normalInd = revise_normal_ind(mus, normalInd, ampParamInds)
	#reclassify clusters
	hetDelParamInds, homDelParamInds, ampParamInds = classify_clusters_given_normal(mus, normalInd)

	clonalHetDelParamInd = determine_clonal_heterozygous_deletion(mus, normalInd, hetDelParamInds, homDelParamInds)
	
	return hetDelParamInds, clonalHetDelParamInd, homDelParamInds, ampParamInds, normalInd

def revise_normal_ind(mus, normalInd, ampParamInds):
	"""
	Revises the decision of what has been called the normal cluster by checking to see if another cluster
	lies on the line of heterozygous deletions that is more likely the normal.

	Arguments:
		mus (list of lists of floats): List of cluster means
		normalInd (int): Index of the cluster that has been guessed to be the normal cluster.
		ampParamInds (list of ints): List of indices of clusters that have been guessed to be
									 amplifications.

	Returns:
		normalInd (int): Index of the cluster that has been classified as the normal cluster.
	"""

	normalRDR = mus[normalInd][0]
	normalBAF = mus[normalInd][1]

	leftx = normalRDR * 0.5
	lefty = 0.5

	#slope of guessed line of heterozygous deletions
	m0 = (normalBAF - lefty) / (normalRDR - leftx)
	#y-intercept of guessed line of heterozygous deletions
	b0 = normalBAF - (m0 * normalRDR)
	#slope of lines perpendicular to line of heterozygous deletions
	m1 = -(m0**-1)

	def score_for_normal(mu, i):
		"""
		Scoring function. The score assigned to a point is the sum of log BAF and its 
		distance to the guessed line of heterozygous deletions.
		"""

		#clusters that are not guessed to be normal or amplifications are disqualified
		if i != normalInd and i not in ampParamInds:
			return float('inf')

		RDR = mu[0]
		BAF = mu[1]
		#y-intercept of point to be scored
		b1 = BAF - (m1 * RDR)
		#x coordinate of point on the line of heterozygous deletions closest to the point being scored
		contactx = (b1 - b0) / (m0 - m1)
		#y coordinate of point on the line of heterozygous deletions closest to the point being scored
		contacty = (m0 * contactx) + b0
		#distance from point being scored to the line of heterozygous deletions
		distToContact = euclidean([RDR, BAF], [contactx, contacty])
		score = distToContact + log(BAF)
		return score

	scores = [score_for_normal(mu, i) for (i, mu) in enumerate(mus)]
	#new normal cluster is that with the lowest score
	normalInd = np.argmin(scores)

	return normalInd

def determine_clonal_heterozygous_deletion(mus, normalInd, hetDelParamInds, homDelParamInds):
	"""
	Determines which cluster is the clonal heterozygous deletion.

	Arguments:
		mus (list of lists of floats): List of cluster means
		normalInd (int): Index of the cluster that has been classified as the normal cluster.
		hetDelParamInds (list of ints): List of indices of clusters that have been classified as 
										heterozygous deletions.
		homDelParamInds (list of ints): List of indices of clusters that have been classified as 
										homozygous deletions.
	"""

	normalRDR = mus[normalInd][0]
	normalBAF = mus[normalInd][1]
	leftx = normalRDR * 0.5
	lefty = 0.5

	#slope of line of heterozygous deletions
	m0 = (normalBAF - lefty) / (normalRDR - leftx)
	#y-intercept of line of heterozygous deletions
	b0 = normalBAF - (m0 * normalRDR)
	#slope of lines perpendicular to the line of heterozygous deletions
	m1 = -(m0**-1)

	def score_for_clonal_het_del(mu, i):
		"""
		Scoring function. The score assigned to a point is the sum of 
		the distance from the point to the line of heterozygous deletions and
		the distance from the point to the y-intercept of the line of
		heterozygous deletions.
		"""

		if i not in hetDelParamInds and i not in homDelParamInds:
			return float('inf')

		RDR = mu[0]
		BAF = mu[1]
		#y-intercept of the point to be scored
		b1 = BAF - (m1 * RDR)
		#x coordinate of point on the line of heterozygous deletions closest to the point being scored
		contactx = (b1 - b0) / (m0 - m1)
		#y coordinate of point on the line of heterozygous deletions closest to the point being scored
		contacty = (m0 * contactx) + b0
		#distance from point being scored to the line of heterozygous deletions
		distToContact = euclidean([RDR, BAF], [contactx, contacty])
		#distance from point being scored to the y-intercept of the line of heterozygous deletions.
		distToIntercept = euclidean([RDR, BAF], [0.0, b0])
		score = distToContact + distToIntercept
		return score

	scores = [score_for_clonal_het_del(mu, i) for (i, mu) in enumerate(mus)]
	#clonal heterozygous deletion cluster is that with the lowest score
	clonalHetDelParamInd = np.argmin(scores)
	return clonalHetDelParamInd

def classify_clusters_given_normal(mus, normalInd):
	"""
	Classifies clusters given that one cluster has been classified as normal.

	Arguments:
		mus (list of lists of floats): List of cluster means
		normalInd (int): Index of the cluster that has been classified as the normal cluster.

	Returns:
		hetDelParamInds (list of ints): List of indices of clusters that have been classified as 
										heterozygous deletions.
		homDelParamInds (list of ints): List of indices of clusters that have been classified as 
										homozygous deletions.
		ampParamInds (list of ints): List of indices of clusters that have been classified as
									 amplifications.
	"""

	normMuX = mus[normalInd][0]
	normMuY = mus[normalInd][1]

	delParamInds = []
	ampParamInds = []
	for i in range(len(mus)):
		if i == normalInd: continue
		
		#deletions have lower RDR than the normal; amplifications have greater RDR than the normal
		if mus[i][0] < normMuX:
			delParamInds.append(i)
		else:
			ampParamInds.append(i)

	hetDelParamInds = []
	homDelParamInds = []
	for i in delParamInds:
		muX = mus[i][0]
		muY = mus[i][1]

		#condition for deciding if a deletion is homozygous (note: needs to be revised)
		if muX < normMuX - 0.2 and muY < normMuY + 0.1:
			homDelParamInds.append(i)
		else:
			hetDelParamInds.append(i)

	return hetDelParamInds, homDelParamInds, ampParamInds

def process_classifications(intervals, missingData, metaMu, clusterAssignments, numClusters, normalInd, clonalHetDelParamInd, hetDelParamInds, ampParamInds, sampleName, outdir):
	"""
	Processes intervals using clustering information by assigning intervals to meta-intervals,
	determining meta-interval length, and deciding copy number bounds for each meta-interval.

	Arguments:
		intervals (list of lists): 2D list containing data related to each interval.
		missingData (list of lists): 2D list containing data related to intervals that have been considered invalid.
		metaMu (list of lists of floats): The means of each cluster.
		clusterAssignments (list of ints): The assignment of each interval to a cluster, where an entry
											j at index i means the ith interval has been assigned to the
											jth meta-interval.
		numClusters (int): The number of clusters.
		normalInd (int): The index of the cluster that has been called normal.
		clonalHetDelParamInd (int): Index of the cluster that has been classified as the clonal 
									heterozygous deletion.
		hetDelParamInds (list of ints): List of indices of clusters that have been classified as 
										heterozygous deletions.
		ampParamInds (list of ints): List of indices of clusters that have been classified as
									 amplifications.
		sampleName (string): The name of the input sample.
		outdir (string): Target directory for output files.

	Returns:
		lengths (list of ints): Length of each interval
		tumorCounts (list of ints): Tumor count for each interval
		normalCounts (list of ints): Normal count for each interval
		upper_bounds (list of ints): The upper copy number bound assigned to each interval.
		lower_bounds (list of ints): The lower copy number bound assigned to each interval.
		fullClusterAssignments (list of ints): The assignment of each interval to a cluster, where an entry
											j at index i means the ith interval has been assigned to the
											jth meta-interval. Includes entries for intervals marked as invalid.
		m (int): The total number of intervals (including intervals that have been marked as invalid).
	"""

	print "Determining copy number bounds..."
	normMu = metaMu[normalInd]
	normRDR = normMu[0]

	if hetDelParamInds != []:
		clonalHetDelRDR = metaMu[clonalHetDelParamInd][0]
		stepSize = normRDR - clonalHetDelRDR
	else:
		clonalHetDelRDR = 0.0
		stepSize = 0.5

	if ampParamInds != []:
		ampMus = [metaMu[ind] for ind in ampParamInds] #amplification means
		ampDistances = [mu[0] - normRDR for mu in ampMus] #difference between amplification RDRs and the normal RDR
		amp_upper = [ceil(distance / stepSize) + 2 for distance in ampDistances] #upper bound on copy numbers for amplifications
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

	plot_clusters(intervals, clusterAssignments, numClusters, sampleName, amp_upper, stepSize, normRDR, clonalHetDelRDR, outdir)

	return lengths, tumorCounts, normalCounts, upper_bounds, lower_bounds, fullClusterAssignments, m


def group_to_meta_interval(lengths, tumorCounts, normalCounts, m, upper_bounds, lower_bounds, clusterAssignments, numClusters):
	"""
	Produces meta-intervals given intervals and their clustering.

	Arguments:
		lengths (list of ints): Length of each interval
		tumorCounts (list of ints): Tumor count for each interval
		normalCounts (list of ints): Normal count for each interval
		m (int): The total number of intervals (including intervals that have been marked as invalid).
		upper_bounds (list of ints): The upper copy number bound assigned to each interval.
		lower_bounds (list of ints): The lower copy number bound assigned to each interval.
		clusterAssignments (list of ints): The assignment of each interval to a cluster, where an entry
											j at index i means the ith interval has been assigned to the
											jth meta-interval. Includes entries for intervals marked as invalid.
		numClusters (int): The number of clusters.
	
	Returns:
		intervalMap (dictionary (int, list of ints)): Maps from a meta-interval index to a list of indices of intervals
														that have been assigned to that meta-interval
		metaLengths (list of ints): "Length" of each meta-interval, where the length is the sum of the lengths
								of the intervals assigned to that meta-interval.
		metaTumorCounts (list of ints): "Tumor counts" of each meta-interval, where the count is the sum of
									the tumor counts of the intervals assigned to that meta-interval.
		metaNormalCounts (list of ints): "Normal counts" of each meta-interval, where the count is the sum of
										the tumor counts of the intervals assigned to that meta-interval.
		meta_lower_bounds (list of ints): The lower copy number bound assigned to each meta-interval.
		meta_upper_bounds (list of ints): The upper copy number bound assigned to each meta-interval.
	"""

	metaLengths = [0 for i in range(numClusters)]
	metaTumorCounts = [0 for i in range(numClusters)]
	metaNormalCounts = [0 for i in range(numClusters)]
	meta_lower_bounds = [2 for i in range(numClusters)]
	meta_upper_bounds = [2 for i in range(numClusters)]
	intervalMap = {}
	
	for val in range(numClusters):
		intervalMap[val] = []
	intervalMap[-1] = []

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
	"""
	Writes a summary file of clustering for a large set of samples.

	Arguments:
		sampleList (list of strings): A list of samples to cluster.
	"""

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