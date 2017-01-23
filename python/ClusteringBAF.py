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

 # @author Layla Oesper, Ahmad Mahmoody, Benjamin J. Raphael, Gryte Satas, and Alex Ashery
 ###

import os
os.environ["BNPYOUTDIR"] = "./"
import bnpy
from FileIO import *
from shutil import rmtree
from ClusterPlottingTools import *
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, ceil, log
from scipy.spatial.distance import euclidean
from multiprocessing import Pool

def clustering_BAF(n, intervals=None, missingData=None, filename=None, byChrm=True, generateData=True, prefix=None, outdir="./", numProcesses=1):
	"""
	Performs clustering on interval data and analyzes clusters to determine copy
	number bounds and which intervals belong to similar clusters.

	Arguments:
		n (int): The number of components to consider.
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
		normalCounts (list of ints): normal count for each interval
		m (int): The total number of intervals (including intervals that have been marked as invalid).
		upper_bounds (list of ints): The upper copy number bound assigned to each interval.
		lower_bounds (list of ints): The lower copy number bound assigned to each interval.
		clusterAssignments (list of ints): The assignment of each interval to a cluster, where an entry
											j at index i means the ith interval has been assigned to the
											jth meta-interval.
		numClusters (int): The number of clusters.
		metaMu (list of lists of floats): The means of each cluster.
		diploidInd (int): The index of the cluster that has been called diploid.
	"""

	#ensure / at the end of output directory
	#if not outdir.endswith("/"):
	#	outdir += "/"

	if prefix is None:
		sampleName = os.path.basename(filename).split(".")[0]
		prefix = sampleName
	else:
		sampleName = prefix

	#setting bnpy outdir
	#bnpyoutdir = outdir + prefix + "_cluster_data/"
	bnpyoutdir = os.path.join(outdir,prefix+"_"+str(n)+"_cluster_data/")
	os.environ["BNPYOUTDIR"] = bnpyoutdir

	if intervals is None and missingData is None:
		missingData, intervals = read_interval_RD_BAF_file(filename, byChrm=byChrm)

	metaData = generate_meta_data(intervals, byChrm, numProcesses, sampleName, generateData, outdir)

	#flatten interval data if it has been grouped by chromosome
	if byChrm:
		intervals = [row for subData in intervals for row in subData]

	print "Begin meta clustering..."
	metaMu, metaSigma, clusterAssignments, numPoints, numClusters = cluster(metaData, sampleName, sf=0.01, intervals=intervals)

	intervalLengths = [row[2] - row[1] + 2 for row in intervals]
	singleCopyParamInds, clonalsingleCopyParamInd, zeroCopyParamInds, ampParamInds, diploidInd = classify_clusters(metaMu, intervalLengths, clusterAssignments)

	plot_classifications(metaMu, metaSigma, intervals, clusterAssignments, numClusters, sampleName, singleCopyParamInds, zeroCopyParamInds, ampParamInds, diploidInd, outdir)

	lengths, tumorCounts, normalCounts, upper_bounds, lower_bounds, clusterAssignments, m = process_classifications(intervals, missingData, metaMu, clusterAssignments, numClusters, diploidInd, clonalsingleCopyParamInd, singleCopyParamInds, ampParamInds, sampleName, outdir)

	#Layla 9-15-15 -- Commenting out this section fixes the crash problem...need to investigate further
	try:
		rmtree(bnpyoutdir)
	except OSError:
		print "WARNING: Was unable to remove bnpy output. This can be manually removed after THetA has completed. bnpy output has been stored in " + bnpyoutdir

	return lengths, tumorCounts, normalCounts, m, upper_bounds, lower_bounds, clusterAssignments, numClusters, metaMu, diploidInd

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

		fig.savefig(os.path.join(outdir,sampleName + "_by_chromosome.png"))
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
		numPoints = [(row[2] - row[1] + 1) / 100000 for row in binnedChrm]
		binnedChrm = generate_data(means, numPoints, sdx=0.02, sdy=0.02)
	else:
		binnedChrm = [row[5:7] for row in binnedChrm]

	mus, sigmas, clusterAssignments, numPoints, numClusters = cluster(binnedChrm, sampleName, chrm=chrm)
	metaDataRow = generate_data(mus, numPoints)

	return binnedChrm, mus, sigmas, clusterAssignments, metaDataRow

def generate_data(mus, numPoints, sdx=0.05, sdy=0.05):
	"""
	Randomly generates diploidly distributed data points about several means.

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
		singleCopyParamInds (list of ints): List of indices of clusters that have been classified as
										single copy states.
		clonalsingleCopyParamInd (int): Index of the cluster that has been classified as the clonal
									single copy state.
		zeroCopyParamInds (list of ints): List of indices of clusters that have been classified as
										zero copy states.
		ampParamInds (list of ints): List of indices of clusters that have been classified as
									 amplifications.
		diploidInd (int): Index of the cluster that has been classified as the diploid cluster.
	"""

	print "Classifying clusters..."

	#generate lengths for meta-intervals
	metaLengths = [0 for i in range(len(mus))]
	for length, assignment in zip(lengths, clusterAssignments):
		if length is not None: metaLengths[assignment] += length

	#guess the diploid cluster to be the largest cluster that has BAF less than 0.2
	meanBAFs = [x[1] for x in mus]
	filteredLengths = map(lambda (BAF, length): -float('inf') if BAF > 0.2 else length, zip(meanBAFs, metaLengths))
	diploidInd = np.argmax(filteredLengths)

	#classify clusters given naive guess
	singleCopyParamInds, zeroCopyParamInds, ampParamInds = classify_clusters_given_diploid(mus, diploidInd)
	#revise decision of diploid cluster
	diploidInd = revise_diploid_ind(mus, diploidInd, ampParamInds)
	#reclassify clusters
	singleCopyParamInds, zeroCopyParamInds, ampParamInds = classify_clusters_given_diploid(mus, diploidInd)

	clonalsingleCopyParamInd = determine_clonal_single_copy_state(mus, diploidInd, singleCopyParamInds, zeroCopyParamInds)

	return singleCopyParamInds, clonalsingleCopyParamInd, zeroCopyParamInds, ampParamInds, diploidInd

def revise_diploid_ind(mus, diploidInd, ampParamInds):
	"""
	Revises the decision of what has been called the diploid cluster by checking to see if another cluster
	lies on the line of single copy states that is more likely the diploid.

	Arguments:
		mus (list of lists of floats): List of cluster means
		diploidInd (int): Index of the cluster that has been guessed to be the diploid cluster.
		ampParamInds (list of ints): List of indices of clusters that have been guessed to be
									 amplifications.

	Returns:
		diploidInd (int): Index of the cluster that has been classified as the diploid cluster.
	"""

	diploidRDR = mus[diploidInd][0]
	diploidBAF = mus[diploidInd][1]

	leftx = diploidRDR * 0.5
	lefty = 0.5

	#slope of guessed line of single copy states
	m0 = (diploidBAF - lefty) / (diploidRDR - leftx)
	#y-intercept of guessed line of single copy states
	b0 = diploidBAF - (m0 * diploidRDR)
	#slope of lines perpendicular to line of single copy states
	m1 = -(m0**-1)

	def score_for_diploid(mu, i):
		"""
		Scoring function. The score assigned to a point is the sum of log BAF and its
		distance to the guessed line of single copy states.
		"""

		#clusters that are not guessed to be diploid or amplifications are disqualified
		if i != diploidInd and i not in ampParamInds:
			return float('inf')

		RDR = mu[0]
		BAF = mu[1]
		#y-intercept of point to be scored
		b1 = BAF - (m1 * RDR)
		#x coordinate of point on the line of single copy states closest to the point being scored
		contactx = (b1 - b0) / (m0 - m1)
		#y coordinate of point on the line of single copy states closest to the point being scored
		contacty = (m0 * contactx) + b0
		#distance from point being scored to the line of single copy states
		distToContact = euclidean([RDR, BAF], [contactx, contacty])
		score = distToContact + log(BAF)
		return score

	scores = [score_for_diploid(mu, i) for (i, mu) in enumerate(mus)]
	#new diploid cluster is that with the lowest score
	diploidInd = np.argmin(scores)

	return diploidInd

def determine_clonal_single_copy_state(mus, diploidInd, singleCopyParamInds, zeroCopyParamInds):
	"""
	Determines which cluster is the clonal single copy state.

	Arguments:
		mus (list of lists of floats): List of cluster means
		diploidInd (int): Index of the cluster that has been classified as the diploid cluster.
		singleCopyParamInds (list of ints): List of indices of clusters that have been classified as
										single copy states.
		zeroCopyParamInds (list of ints): List of indices of clusters that have been classified as
										zero copy states.
	"""

	diploidRDR = mus[diploidInd][0]
	diploidBAF = mus[diploidInd][1]
	leftx = diploidRDR * 0.5
	lefty = 0.5

	#slope of line of single copy states
	m0 = (diploidBAF - lefty) / (diploidRDR - leftx)
	#y-intercept of line of single copy states
	b0 = diploidBAF - (m0 * diploidRDR)
	#slope of lines perpendicular to the line of single copy states
	m1 = -(m0**-1)

	def score_for_clonal_single_copy(mu, i):
		"""
		Scoring function. The score assigned to a point is the sum of
		the distance from the point to the line of single copy states and
		the distance from the point to the y-intercept of the line of
		single copy states.
		"""

		if i not in singleCopyParamInds and i not in zeroCopyParamInds:
			return float('inf')

		RDR = mu[0]
		BAF = mu[1]
		#y-intercept of the point to be scored
		b1 = BAF - (m1 * RDR)
		#x coordinate of point on the line of single copy states closest to the point being scored
		contactx = (b1 - b0) / (m0 - m1)
		#y coordinate of point on the line of single copy states closest to the point being scored
		contacty = (m0 * contactx) + b0
		#distance from point being scored to the line of single copy states
		distToContact = euclidean([RDR, BAF], [contactx, contacty])
		#distance from point being scored to the y-intercept of the line of single copy states.
		distToIntercept = euclidean([RDR, BAF], [0.0, b0])
		score = distToContact + distToIntercept
		return score

	scores = [score_for_clonal_single_copy(mu, i) for (i, mu) in enumerate(mus)]
	#clonal single copy state cluster is that with the lowest score
	clonalsingleCopyParamInd = np.argmin(scores)
	return clonalsingleCopyParamInd

def classify_clusters_given_diploid(mus, diploidInd):
	"""
	Classifies clusters given that one cluster has been classified as diploid.

	Arguments:
		mus (list of lists of floats): List of cluster means
		diploidInd (int): Index of the cluster that has been classified as the diploid cluster.

	Returns:
		singleCopyParamInds (list of ints): List of indices of clusters that have been classified as
										single copy states.
		zeroCopyParamInds (list of ints): List of indices of clusters that have been classified as
										zero copy states.
		ampParamInds (list of ints): List of indices of clusters that have been classified as
									 amplifications.
	"""

	diploidMuX = mus[diploidInd][0]
	diploidMuY = mus[diploidInd][1]

	delParamInds = []
	ampParamInds = []
	for i in range(len(mus)):
		if i == diploidInd: continue

		#deletions have lower RDR than the diploid; amplifications have greater RDR than the diploid
		if mus[i][0] < diploidMuX:
			delParamInds.append(i)
		else:
			ampParamInds.append(i)

	singleCopyParamInds = []
	zeroCopyParamInds = []
	for i in delParamInds:
		muX = mus[i][0]
		muY = mus[i][1]

		#condition for deciding if a deletion is a zero copy state (note: needs to be revised)
		if muX < diploidMuX - 0.2 and muY < diploidMuY + 0.1:
			zeroCopyParamInds.append(i)
		else:
			singleCopyParamInds.append(i)

	return singleCopyParamInds, zeroCopyParamInds, ampParamInds

def process_classifications(intervals, missingData, metaMu, clusterAssignments, numClusters, diploidInd, clonalsingleCopyParamInd, singleCopyParamInds, ampParamInds, sampleName, outdir):
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
		diploidInd (int): The index of the cluster that has been called diploid.
		clonalsingleCopyParamInd (int): Index of the cluster that has been classified as the clonal
									single copy state.
		singleCopyParamInds (list of ints): List of indices of clusters that have been classified as
										single copy states.
		ampParamInds (list of ints): List of indices of clusters that have been classified as
									 amplifications.
		sampleName (string): The name of the input sample.
		outdir (string): Target directory for output files.

	Returns:
		lengths (list of ints): Length of each interval
		tumorCounts (list of ints): Tumor count for each interval
		normalCounts (list of ints): normal count for each interval
		upper_bounds (list of ints): The upper copy number bound assigned to each interval.
		lower_bounds (list of ints): The lower copy number bound assigned to each interval.
		fullClusterAssignments (list of ints): The assignment of each interval to a cluster, where an entry
											j at index i means the ith interval has been assigned to the
											jth meta-interval. Includes entries for intervals marked as invalid.
		m (int): The total number of intervals (including intervals that have been marked as invalid).
	"""

	print "Determining copy number bounds..."
	diploidMu = metaMu[diploidInd]
	diploidRDR = diploidMu[0]

	if singleCopyParamInds != []:
		clonalsingleCopyRDR = metaMu[clonalsingleCopyParamInd][0]
		stepSize = diploidRDR - clonalsingleCopyRDR
	else:
		clonalsingleCopyRDR = 0.0
		stepSize = 0.5

	if ampParamInds != []:
		ampMus = [metaMu[ind] for ind in ampParamInds] #amplification means
		ampDistances = [mu[0] - diploidRDR for mu in ampMus] #difference between amplification RDRs and the diploid RDR
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
				if clusterAssignments[j] == diploidInd:
					lower_bounds[i] = 2
				elif clusterAssignments[j] in singleCopyParamInds:
					lower_bounds[i] = 1
				else:
					lower_bounds[i] = 0
			j += 1

	plot_clusters(intervals, clusterAssignments, numClusters, sampleName, amp_upper, stepSize, diploidRDR, clonalsingleCopyRDR, outdir)

	return lengths, tumorCounts, normalCounts, upper_bounds, lower_bounds, fullClusterAssignments, m


def group_to_meta_interval(lengths, tumorCounts, normalCounts, m, upper_bounds, lower_bounds, clusterAssignments, numClusters):
	"""
	Produces meta-intervals given intervals and their clustering.

	Arguments:
		lengths (list of ints): Length of each interval
		tumorCounts (list of ints): Tumor count for each interval
		normalCounts (list of ints): normal count for each interval
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
		metaNormalCounts (list of ints): "normal counts" of each meta-interval, where the count is the sum of
										the normal counts of the intervals assigned to that meta-interval.
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
			lengths, tumorCounts, normalCounts, m, upper_bounds, lower_bounds, clusterAssignments, numClusters, clusterMeans, diploidInd = clustering_BAF(filename)
			singleCopyParamInds, zeroCopyParamInds, ampParamInds, diploidInd = classify_clusters(clusterMeans, lengths, clusterAssignments)
			intervalMap, metaLengths, metaTumorCounts, metaNormalCounts, meta_lower_bounds, meta_upper_bounds = group_to_meta_interval(lengths, tumorCounts, normalCounts, m, upper_bounds, lower_bounds, clusterAssignments, numClusters)
			f.write(sample + "\n")
			f.write(str(numClusters) + "\n")
			for i in range(numClusters):
				f.write(str(metaLengths[i]) + "\t" + str(clusterMeans[i][0]) + "\t" + str(clusterMeans[i][1]) + "\t")
				if i == diploidInd:
					f.write("DIPLOID")
				elif i in singleCopyParamInds:
					f.write("SINGLE")
				elif i in zeroCopyParamInds:
					f.write("ZERO")
				else:
					f.write("AMP")
				f.write("\n")
			f.write("\n")
			f.flush()
		except IOError:
			continue
	f.close()
