# This file contains different methods for setting bounds on intervals.
import numpy as np
import bisect
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema

def find_lt(a, x):
    'Find the index of the bin where x would be inserted into a'
    i = bisect.bisect_left(a, x)
    if i:
        return i-1
    print "\n"
    raise ValueError



def set_new_bounds(new_bounds_file):
	"""
	Takes in a RD and BAF file.

	Returns:
	upper_bounds, lower_bounds
	"""
	upper_bounds = []
	lower_bounds = []

	#load in file
	cols=(1,2,5,6)
	X=np.loadtxt(new_bounds_file, usecols=cols)
	goodRows = np.logical_and(X[:,2] != -1, X[:,3] != -1)

	#Create weighted histogram of points
	Y=None
	cov_mat = [[0.002,0],[0,0.002]]
	for row in X[goodRows]:
		length = row[1] - row[0] +1
		if length < 1000000 or row[2] > 3: continue
		num_points = int(round(length/100000))
		if num_points == 0: num_points =1

		mean = [row[2],row[3]]
		pts = np.random.multivariate_normal(mean, cov_mat, num_points)

		if Y is None:
			Y = pts
		else:
			Y = np.concatenate((Y, pts), axis=0)


	# Use KDE to create bins
	rd_vals = Y[:,0]
	x_grid = np.linspace(0,3,1000)
	kde = gaussian_kde(rd_vals)
	pdf = kde.evaluate(x_grid)

	min_pts = argrelextrema(pdf,np.less)
	#max_pts = argrelextrema(pdf, np.greater)
	max_idx = np.argmax(pdf)
	max_x_pt = x_grid[max_idx]

	bins = x_grid[min_pts].tolist()
	bins.append(0)
	bins.append(3)
	bins.sort()

	cluster_assignment = []

	for row in X:
		if row[2] == -1 or row[3] == -1 or (row[1] - row[0] + 1) < 1000000 or row[2] > 3:
			cluster_assignment.append(-1)
		else:
			cluster_assignment.append(find_lt(bins,row[2]))

	norm_cluster = find_lt(bins,max_x_pt)
	cluster_bounds = get_cluster_bounds(cluster_assignment,bins,norm_cluster)

	for v in cluster_assignment:
		(lb,ub) = cluster_bounds[v]
		lower_bounds.append(lb)
		upper_bounds.append(ub)
	#upper_bounds = None
	#lower_bounds = None

	cluster_properties = get_cluster_rd_baf(cluster_assignment, X)
	numClusters=len(bins)-1
	
	#for key in cluster_properties.keys(): print cluster_properties[key]

	return upper_bounds, lower_bounds, cluster_assignment, numClusters


def get_cluster_bounds(cluster_assignment,bins,norm_cluster):
	"""
	Returns a dictonary of cluster ids to bounds.
	"""
	cluster_bounds = dict()
	num_bins = len(bins) - 1

	for i in xrange(num_bins):

		if i < norm_cluster:
			cluster_bounds[i]=(1,2)

		elif i == norm_cluster:
			cluster_bounds[i]=(2,2)

		else:
			cluster_bounds[i]=(2,3)

	cluster_bounds[-1]=('X','X')

	return cluster_bounds


def get_cluster_rd_baf(cluster_assignments, X):
	"""
	Returns the average read depth and baf for all intervals
	assigned to the different clusters.
	"""

	cluster_ids = set(cluster_assignments)
	cluster_properties = dict()
	for id in cluster_ids:
		cluster_properties[id]=(0,0,0)

	for i,row in enumerate(cluster_assignments):
		rd,baf,count=cluster_properties[row]
		rd+=X[i,2]
		baf+=X[i,3]
		count+=1
		cluster_properties[row]=(rd, baf, count)

	for key in cluster_properties.keys():
		rd,baf, count = cluster_properties[key]
		if count != 0:
			mean_rd = rd/float(count)
			mean_baf = baf/float(count)
			cluster_properties[key] = (mean_rd,mean_baf)
		else:
			cluster_properties[key] = (-1,-1)

	return cluster_properties

	

