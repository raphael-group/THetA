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
from math import *
from os.path import basename, join
import matplotlib.pyplot as plt
import sys
from numpy import linspace
from scipy.stats import norm, beta
from multiprocessing import Pool

def run_BAF_model(resultsFile, tumor=None, normal=None, tumorBAF=None, normalBAF=None, chrmsToUse=None, intervals=None, tumorSNP=None, normalSNP=None, intervalFile=None, prefix=None, directory="./", plotOption="best", model="gaussian", width=12.0, height=12.0, gamma=0.05, numProcesses=1):
	"""
	Runs the BAF model on SNP and interval data.

	Arguments:
		tumorSNP (string): location of tumor SNP data
		normalSNP (string): location of normal SNP data
		intervalFile (string): location of original input file for THetA
		resultsFile (string): location of original output file for THetA
		prefix (string): prefix used for output file
		directory (string): directory where the output file is saved
		plotOption (string): option for plotting; "all" to plot all results,
								"best" to plot the optimal result only
		model (string): distribution used for the BAF model (note: only "gaussian"
						is currently supported)
		width (float): width of output plot in inches
		height (float): height of output plot in inches
		gamma (float): parameter for determining heterozygosity of an allele
	"""

	minSize = 2000000 #lower bound on interval size
	minSNP = 10 #lower bound on number of reads

	if tumor is None and normal is None and tumorBAF is None and normalBAF is None and chrmsToUse is None:
		#reading in tumor BAF data
		tumor = read_snp_file(tumorSNP)

		#reading in normal snp data, interval data, and THetA results
		normal = read_snp_file(normalSNP)
		chrmsToUse, intervals = read_interval_file_BAF(intervalFile)

		#calculate the BAFs and remove irrelevant data
		tumorBAF, normalBAF, tumor, normal = calculate_BAF(tumor, normal, chrmsToUse, minSNP, gamma, numProcesses)

	results = read_results_file_full(resultsFile)
	
	#possible parameters estimated by THetA
	k = results['k']
	C = results['C']
	mu = results['mu']

	#arrays for storing results
	BAFVec = []
	meansVec = []
	posVec = []
	chrmVec = []
	NLLVec = []

	#run model on each THetA output
	for i in range(k):
		print "Calculating NLL for model " + str(i + 1)

		#unpack model parameters
		currC = C[i]
		currMu = mu[i]

		#filter out intervals that are too small or have invalid entries in C matrix
		filteredData = filter(lambda ((chrm, start, end), cj): (end - start + 1) >= minSize and -1 not in cj, zip(intervals, currC))
		currIntervals, currC = zip(*filteredData)
		currIntervals = list(currIntervals)
		currC = list(currC)

		#pi maps from chromosome and position to interval number
		pi = generate_pi(currIntervals)

		#calculate the NLL for the given model
		if model == "gaussian":
			currBAF, currMeans, currPos, currChrmVec, currNLL = get_gaussian_NLL(tumor, tumorBAF, normal, normalBAF, currC, currMu, pi, numProcesses)
		else:
			print model + " is not a supported model. Exiting program..."
			exit(1)

		#storing results
		BAFVec.append(currBAF)
		meansVec.append(currMeans)
		posVec.append(currPos)
		chrmVec.append(currChrmVec)
		NLLVec.append(currNLL)

	plotDim = (width, height)

	if prefix == None:
		prefix = ".".join(basename(resultsFile).split(".")[0:2])

	#generate results
	plot_results(BAFVec, meansVec, posVec, chrmVec, NLLVec, chrmsToUse, plotOption, directory, prefix, plotDim)
	results['BAF_NLL'] = NLLVec
	write_out_NLL_result(directory, prefix, results)

def plot_single_result(BAF, means, pos, chrm, NLL, chrmsToUse, numberResults, fig, colors, plotNum=0):
	"""
	Creates a subplot from the BAF model calculations for one THetA output.

	Arguments:
		BAF (list of floats): BAFs calculated by the model
		means (list of floats): means calculated by the model
		pos (list of ints): interval number for the SNP associated with each BAF
		chrm (list of ints): chromosome number for the SNP associated with each BAF
		NLL (float): negative log likelihood calculated by the model
		chrmsToUse (list of ints): a list of chromosomes on which the BAF model was run
		numberResults (int): number of results that theta output
		fig (matplotlib figure): figure used for plotting the data
		colors (matplotlib colors): colors used for distinguishing plot points by chromsome number]
		i (int): result number, used to choose which subplot the plot is placed on

	Returns:
		fig with a new subplot added for result i
	"""

	print "Plotting model " + str(plotNum)

	#create new subplot
	ax = fig.add_subplot(numberResults, 1, plotNum)
	
	#parameter for plotting chromosome separations
	mag = 6

	dataArray = zip(BAF, means, pos, chrm)
	dataDict = dict(zip(chrmsToUse, map(lambda x: [], chrmsToUse)))
	
	#sorting by chromosome number
	for row in dataArray:
		dataDict[row[3]].append(row[:3])
		
	#sets offset for each chromosome
	offset = 0

	#stores location for setting xticks to chromosome numbers
	xlabelPoints = []

	#ensuring chromosomes are in sorted order
	chrmsToUse = sorted(chrmsToUse)

	#plotting
	for chrm in chrmsToUse:
		relevantData = dataDict[chrm]
		
		xs = []
		ys = []
		mus = []

		color = colors[chrm - 1]

		#determines offset for next chromosome
		maxPos = offset
		for row in relevantData:
			BAF = row[0]
			mean = row[1]
			pos = row[2]
			x = pos + offset

			xs.append(x)
			ys.append(BAF)
			mus.append(mean)

			maxPos = x if x > maxPos else maxPos

		#plotting chromosome points
		xlabelPoints.append((offset + maxPos) / 2.0)
		offset = maxPos + (2 * (10**mag))
		ax.plot(xs, ys, 'o', color=color, ms=2, markeredgecolor='none', zorder=1, solid_capstyle="butt")
		ax.plot(xs, mus, 's', color='black', ms=2, zorder=2, solid_capstyle="butt")
		ax.plot([maxPos + (10**mag), maxPos + (10**mag)], [0, 1], color='black', zorder=3, linewidth=2)
			
	#formatting plot
	ax.set_title('BAF Model NLL: ' + str(NLL))
	ax.set_xticks(xlabelPoints)
	ax.set_xticklabels(chrmsToUse)
	ax.set_xlabel('Chromosome')
	ax.set_ylabel('BAF')
	ax.set_xlim([0, maxPos])

	return fig

def plot_results(BAFVec, meansVec, posVec, chrmVec, NLLVec, chrmsToUse, plotOption, directory, prefix, plotDim):
	"""
	Plots results from BAF model.

	Arguments:
		BAFVec (list of lists of floats): BAFs calculated for each THetA output
		meansVec (list of lists of floats): means calculated for each THetA output
		posVec (list of lists of ints): interval numbers for the SNPs associated with the BAFs in each THetA output
		chrmVec (list of lists of ints): chromsome numbers for the SNPs associated with the BAFs in each THetA output
		NLLVec (list of ints): negative log likelihood calculated by the BAF model for each THetA output
		chrmsToUse (list of ints): a list of chromosomes on which the BAF model was run
		plotOption (string): option for plotting; "all" to plot all results,
									"best" to plot the optimal result only
		directory (string): directory where the output file is saved
		prefix (string): prefix used for output file
		plotDim (2-tuple of floats): dimensions of output plot, stored as (width, height)
	"""

	numberChrms = len(chrmsToUse)
	numberResults = len(NLLVec)
	
	#calculating colors used for distinguishing data by chromosome number
	cmap = plt.get_cmap('gist_rainbow')
	colors = [cmap(i) for i in linspace(0, 1, numberChrms)]

	fig = plt.figure(figsize=plotDim)

	if plotOption == "all":
		#plot all results
		for i in range(numberResults):
			fig = plot_single_result(BAFVec[i], meansVec[i], posVec[i], 
									 chrmVec[i], NLLVec[i], chrmsToUse,
									 numberResults, fig, colors, plotNum=i)
	elif plotOption == "best":
		#determine the optimal result
		currBestInd = 0
		currBest = NLLVec[0]
		val, idx = min((val, idx) for (idx, val) in enumerate(NLLVec))

		#plot best result
		fig = plot_single_result(BAFVec[idx], meansVec[idx], posVec[idx], 
									 chrmVec[idx], NLLVec[idx], chrmsToUse,
									 1, fig, colors)
	else:
		print "Plot option not recognized. Exiting..."
		exit(1)

	#formatting and saving plot
	fig.tight_layout()

	#layla - fix plotting location
	fig_file=os.path.join(directory,prefix+".BAF.plot." + plotOption +".png")
	plt.savefig(fig_file)
	#plt.savefig(directory + prefix + ".BAF.plot." + plotOption +".png")

def is_heterozygous((n_a, n_b, gamma)):
	"""
	Determines if an allele should be considered heterozygous.

	Arguments:
		n_a (int): number of a alleles counted. Used as alpha parameter for the beta distribution
		n_b (int): number of b alleles counted. Used as beta parameter for the beta distribution
		gamma (float): parameter used for deciding heterozygosity; determined via a beta distribution
						with 1 - gamma confidence

	Returns:
		A boolean indicating whether or not the allele should be considered heterozygous.
	"""

	if n_a == -1 or n_b == -1: return False

	p_lower = gamma / 2.0
	p_upper = 1 - p_lower

	[c_lower, c_upper] = beta.ppf([p_lower, p_upper], n_a + 1, n_b + 1)
	return c_lower <= 0.5 and c_upper >= 0.5

def calculate_BAF(tumorData, normalData, chrmsToUse, minSNP, gamma, numProcesses):
	"""
	Calculates the BAF for tumor SNP data and normal SNP data. Also filters for useful data using a variety of
	heuristics.

	Arguments:
		tumorData (2D list): 2D list representing the tumor data
		normalData (2D list): 2D list representing the normal data
		chrmsToUse (list of ints): a list of chromosomes on which the BAF model should be run
		minSNP (int): lower bound on number of SNP reads, used to filter out data that doesn't have enough information
		gamma (float): parameter for determining heterozygosity of an allele

	Returns:
		tumorBAF (list of floats): BAF calculated for each tumor SNP
		normalBAF (list of floats): BAF calcualted for each normal SNP
		newTumorData (2D list): tumorData filtered for relevant rows
		newNormalData (2D list): normalData filtered for relevant rows
	"""

	#function to select columns from a 2D list
	select_col = lambda array, colNum: map(lambda x: x[colNum], array)

	#vectors of tumor data
	tumorMutCount = select_col(tumorData, 3)
	tumorRefCount = select_col(tumorData, 2)

	#vectors of normal data
	normalMutCount = select_col(normalData, 3)
	normalRefCount = select_col(normalData, 2)

	#denominators for BAFs
	tumorDenom = map(sum, zip(tumorMutCount, tumorRefCount))
	normalDenom = map(sum, zip(normalMutCount, normalRefCount))

	tumorBAF = []
	normalBAF = []
	newTumorData = []
	newNormalData = []
	print "Determining heterozygosity."
	p = Pool(numProcesses)
	repGamma = [gamma for i in range(len(tumorData))]
	isHet = p.map(is_heterozygous, zip(normalRefCount, normalMutCount, repGamma))
	print "Calculating BAFs."
	for i in range(len(tumorData)):
		chrm = tumorData[i][0]
		#filter out data that uses irrelevant chromosomes
		if chrm not in chrmsToUse: continue
		#filter out data where there aren't enough tumor reads
		if tumorMutCount[i] + tumorRefCount[i] < minSNP: continue
		#filter out data where there aren't enough normal reads
		if normalMutCount[i] + normalRefCount[i] < minSNP: continue

		currTumorNum = tumorMutCount[i]
		currTumorDenom = tumorDenom[i]

		currNormalNum = normalMutCount[i]
		currNormalDenom = normalDenom[i]

		#filter out points with denominators that are 0 to prevent zero division error
		if currTumorDenom == 0 or currNormalDenom == 0: continue
		else:
			tumorBAFj = currTumorNum / currTumorDenom
			normalBAFj = currNormalNum / currNormalDenom
			#filter out data where normal BAFs do not fit in bounds correctly
			if isHet[i]:
				tumorBAF.append(tumorBAFj)
				normalBAF.append(normalBAFj)
				newTumorData.append(tumorData[i])
				newNormalData.append(normalData[i])
			else:
				continue

	return tumorBAF, normalBAF, newTumorData, newNormalData

def generate_delta(C, mu):
	"""
	Generates delta for the BAF gaussian model, which corresponds to the
	distribution mean. See the THetA2 supplement for an explanation of delta 
	and the helper method phi.

	Arguments:
		C (list of lists of ints): copy number count matrix calculated using THetA
		mu (list of floats): tumor population distribution calculated using THetA

	returns:
		delta (list of floats): delta parameters calculated for each interval
	"""

	def phi(a):
		if a == 0:
			return 0.0
		elif a == 3:
			return 2.0
		else:
			return 1.0

	delta = []
	for row in C:
		numerator = sum(map(lambda (a, b): phi(a) * b, zip(row, mu)))
		denominator = sum(map(lambda (a, b): a * b, zip(row, mu)))
		deltaj = (numerator / denominator) - 0.5
		delta.append(deltaj)
	return delta

def generate_pi(intervals):
	"""
	Generates pi for the BAF gaussian model, which maps from SNP location to
	interval number. In this implementation, pi maps from a chromsome number to
	ranges of positions and the interval number associated with that range. See
	the THetA2 supplement for a more in-depth explanation

	Arguments:
		intervals (2D list): 2D list of data taken from the intervals file

	Returns:
		pi (dict of lists of ints): pi map
	"""

	pi = {}
	j = 0
	for chrm, start_pos, end_pos in intervals:
		try:
			pi[chrm].append((start_pos, end_pos, j))
		except KeyError:
			pi[chrm] = [(start_pos, end_pos, j)]
		j += 1

	return pi

def calculate_interval(pi, chrm, pos):
	"""
	Determines how a position is mapped to an interval number using pi.

	Arguments:
		pi (dict of lists of ints): the pi map
		chrm (int): chromosome number
		pos (int): SNP location
	Returns:
		The interval associated with pos on chrm if the mapping is successful, otherwise None
	"""

	try:
		chrmArray = pi[chrm]
		for start, end, ind in chrmArray:
			if start <= pos and pos <= end:
				return ind
		#position was not found in intervals
		return None
	except KeyError:
		#chromosome was not found in pi
		return None

def generate_sigma(normal, normalBAF, pi, m):
	"""
	Generates sigma for the BAF gaussian model, which corresponds to the
	distribution's standard deviation.

	Arguments:
		normal (2D list): 2D list of normal snp data
		normalBAF (list of floats): BAFs calculated from the normal snp data
		pi (dict of lists of ints): pi map
		m: number of intervals

	Returns:
		sigma (list of floats): sigma parameter calculated for each interval
	"""

	numerator = [0] * m
	denominator = [0] * m
	for row, BAF in zip(normal, normalBAF):
		chrm = row[0]
		pos = row[1]
		j = calculate_interval(pi, chrm, pos)
		#ignore data whose interval can't be found
		if j is None: continue

		numerator[j] += (BAF - 0.5)**2
		denominator[j] += 1

	sigma = map(lambda (n, d): n / d if d != 0 else None, zip(numerator, denominator))
	return sigma

def normal_BAF_pdf(x, delta, sigma):
	"""
	Calculates the normal distribution values used in the BAF model.

	Arguments:
		x (float): input for the distribution
		delta (float): delta parameter used to calculate mu
		sigma (float): standard deviation of the distribution

	Returns:
		mu (float): the mean of the distribution
		p (float): the probability of seeing x in the distribution
	"""
	
	#casting values as floats for safety
	x = float(x)
	delta = float(delta)
	sigma = sqrt(float(sigma))

	#sgn is used to calculate mu. See the THetA2 supplement for a full explanation.
	sgn = lambda y: 1.0 if y >= 0 else -1.0
	mu = 0.5 + (sgn(x - 0.5) * delta)
	p = norm(mu, sigma).pdf(x)
	return mu, p

def get_gaussian_NLL(tumor, tumorBAF, normal, normalBAF, C, mu, pi, numProcesses):
	"""
	Calculates the negative log likelihood of seeing a result from THetA using a gaussian distribution.

	Arguments:
		tumor (2D list): A 2D list of tumor snp data
		tumorBAF (list of floats): a list of the BAFs for the tumor sample
		normal (2D list): A 2D list of normal snp data
		normalBAF (list of floats): a list of the BAFs for the normal sample
		chrmsToUse (list of ints): a list of chromosomes on which the BAF model should be run
		C (list of lists of ints): copy number count matrix calculated using THetA
		mu (list of floats): tumor population distribution calculated using THetA
		pi (dict of lists of ints): map for determining the interval that a SNP sits on
	Returns:
		tumorBAF (list of floats): The BAFs calculated for relevant tumor SNPs
		means (list of floats): The mean calculated for each relevant tumor SNP
		posVec (list of ints): The positions corresponding to each tumor SNP that is used
		chrmVec (list of ints): The chromosome number corresponding to each tumor SNP that is used
		NLL (float): the negative log likelihood of seeing C and mu assuming that the BAF in a sample
					 is generated by a gaussian distribution
	"""

	#calculating distribution parameters
	delta = generate_delta(C, mu)
	sigma = generate_sigma(normal, normalBAF, pi, len(C)) #not sure about len(C)

	NLL = 0
	means = []
	posVec = []
	chrmVec = []
	for i in range(len(tumorBAF)):
		chrm = tumor[i][0]
		pos = tumor[i][1]
		j = calculate_interval(pi, chrm, pos)
		#ignore snps whose intervals cannot be calculated, or have a 0 standard deviation
		if j is None or sigma[j] is None or sigma[j] == 0: continue
		else:
			mean, likelihood = normal_BAF_pdf(tumorBAF[i], delta[j], sigma[j])
			NLL -= log(likelihood)
			means.append(mean)
			posVec.append(pos)
			chrmVec.append(chrm)

	return tumorBAF, means, posVec, chrmVec, NLL

def main():
	kwargs = parse_BAF_arguments()
	run_BAF_model(**kwargs)

if __name__ == "__main__":
	main()
