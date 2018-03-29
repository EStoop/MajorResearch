import cPickle
import numpy as np
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u
import scipy.optimize as s


def root(varpi, e_varpi):
	'''
	Calculates the real roots of the equation defined by r^3/L - 2r^2 + (varpi/e_varpi^2)*r - e_varpi^-2 = 0 (Bailer-Jones 2015)
	
	Parameters
	----------
	varpi	The parallax in arcsec
	e_varpi	The uncertainty on the parallax in arcsec
	
	Output
	------
	r		The real roots in varpisec
	'''
	L = (1.00 * u.kpc).to(u.pc) # (Bailer-Jones 2015)
	
	# The coefficience gives p[0]*r^0 + p[1]*r^1+...+ p[n]
	p = [(1./L).value, -2., (varpi/e_varpi**2), -1 * (e_varpi**(-2))] 
	# Calculates the roots
	root = np.roots(p)
	# Only return the real varpit of the returned values
	r = root.real[root.imag == 0.0] 
	return r


def rmode(varpi, e_varpi):
	'''
	Determines the r_mode from the returned values of root. If there is one real root than this is the mode. If there are three roots then: if varpi>0 the smallest root is taken. Otherwise the postive value is r_mode. 
	
	Parameters
	----------
	varpi	The parallax in arcsec
	e_varpi	The uncertainty on the parallax in arcsec

	Output
	------
	rmode	The mode in parsec
	'''
	r = root(varpi, e_varpi)
	if len(r) > 1:
		if varpi > 0.0: rmode = np.min(r)
		else: rmode = r[r > 0.0]
	else: 
		rmode = r[0]
	return rmode

def main(varpi, e_varpi):
	'''
	Returnes the distance belonging to varpi and e_varpi
	
	Parameters
	----------
	varpi	An array with parallax in arcsec
	e_varpi	An array containing uncertainty on the parallax in arcsec

	Output
	------
	dist	An array containing the predicted distances in parsec
	
	'''
	dist = np.zeros_like(varpi.value)
	for i in range(len(varpi)):
		dist[i] = rmode((varpi[i]).value, (e_varpi[i]).value)
	return dist




