def isClose(v1, v2, margin = 10e-4):
	"""
	Check if every pair of floats (v1[i],v2[i] are within the margin value of eachother)
	
	Args:
		v1, v2 (list of floats): lists with len(v1) = len(v2)
		margin (float): margin distance to judge "closeness"
	"""
	for (a,b) in zip(v1,v2):
		if abs(a-b) > margin: 
			return False
	return True

def inRange(v1, minVal = 0, maxVal = 1):
	"""
	Check if every value in v1 is in the range (minVal, maxVal)
	"""
	for v in v1: 
		if v < minVal or v > maxVal: 
			return False

	return True
