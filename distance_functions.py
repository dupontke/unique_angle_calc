#!/mnt/lustre_fs/users/mjmcc/apps/python2.7/bin/python
##!/Library/Frameworks/Python.framework/Versions/2.7/bin/python

# USAGE:
# from distance_functions import *

# PREAMBLE:

import numpy as np

sqrt = np.sqrt
sums = np.sum
square = np.square
zeros = np.zeros

# SUBROUTINES:

def RMSD(x,y,n):
	""" Calculates the Root Mean Squared Distance between two arrays of the same size

	Usage: rmsd = RMSD(x,y,n)

	Arguments:
	x, y: numpy arrays with the same shape (n X 3)
	n: number of particles being summed over; ex: number of atoms in the atom selection being analyzed;
		if n = 1, this function calculates the distance between x and y arrays

	"""
	
	return sqrt(sums(square(x-y))/n)

def MSD(x,y,n):
	""" Calculates the Mean Squared Distance between two arrays of the same size

	Usage: msd = MSD(x,y,n)

	Arguments:
	x, y: numpy arrays with the same shape
	n: number of particles being summed over; ex: number of atoms in the atom selection being analyzed;
		if n = 1, this function calculates the distance squared between x and y arrays

	"""

	return sums(square(x-y))/n

def wrapping(x,dim):
	""" Calculates the translation matrix needed to wrap a particle back into the original periodic box
	
	Usage: t = wrapping(x,dim)

	Arguments:
	x: a numpy array of size (3) that corresponds to the xyz coordinates of an ATOM/COM/COG of a residue
	dim: a numpy array of size (3) that holds the xyz dimensions of the periodic box at that timestep

	"""
	
	t = zeros(3)
	dim2 = dim/2.
	for i in range(3):
		if (x[i]<-dim2[i]) or (x[i]>dim2[i]):
			t[i] = -dim[i]*round(x[i]/dim[i])
	return t

def euclid_dist(x,y):
	""" Calculates the Euclidian Distance between two arrays of the same size
	Usage: dist,dist2 = euclid_dist(x,y)

	Arguments:
	x, y: numpy arrays with the same size
	"""

	dist2 = sums(square(x-y))
	dist = sqrt(dist2)
	return dist, dist2

def computePbcDist2(r1,r2,box):
	""" compute the distance between two points taking into account periodic boundary conditions
	Usage: dist = computePbcDist2(r1,r2,box):

	Arguments:
	r1, r2: two points that are defined
	box: dimensions of the box containing protein and solvent
	"""

	dist2 = 0

	for j in range(0,3):
		temp = r1[j]-r2[j]
		if temp < -box[j]/2.0:
			temp += box[j]
		elif temp > box[j]/2.0:
			temp -= box[j]
		dist2 += temp*temp

	dist2 = math.sqrt(dist2)
	return dist2;
