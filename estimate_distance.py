import cPickle
import numpy as np
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u


# Read the gog sample file
dir_gog = '/net/blekekil/data2/gog/random_samples/'
nSamples = 1837070.

print "\nStart reading data"
data = cPickle.load( open(dir_gog + 'gog_parameters_' + str(int(nSamples)) + '.pkl' ))
edata = cPickle.load( open(dir_gog + 'gog_errors_' + str(int(nSamples)) + '.pkl' ))
print "Data loaded"

parallax = np.array(data['parallax']) * u.mas 			# mas
e_parallax = np.array(edata['parallax_err']) * u.mas 	# mas

ra  = np.array(data['alpha']) * u.rad # rad 
dec = np.array(data['delta']) * u.rad # rad


# Remove parallaxes with relative errors (f_obs) > 1.0
ind = np.where(np.abs(e_parallax/parallax) > 1.)[0]
par, e_par = np.delete(par, ind), np.delete(e_par, ind)

# Likelihood
def Likelihood(r, std, varpi):
	'''
	Gives the likelihood on the parallax, varpi, given the distance r and standard deviation std. 
	----------------
	Input Parameters
	----------------
	r		Distance in parsec.
	std		The observed uncertainty on the parallax in arcseconds
	varpi	The observed parallax in arcseconds. 
	
	------
	Output
	------
	P		The likelihood on varpi given the uncertainty and the distance	
	'''
	
	return (1.0 / np.sqrt(2.0*np.pi)*std) * np.exp( (-1.0/(2.0* std**2)) * (varpi - r**(-1))**2 )

# The stellar number density
def rho_MW(r, l, b):
	
	'''
	Calculates the stellar number density of the Milky Way. The number of stars per volumne element dV. 
	----------------
	Input Parameters
	----------------
	r		Distance in parsec.
	l		The longitude coordiante
	b		The latitude coordinate
	
	------
	Output
	------
	rho_MW	The stellar number density of the Milky Way consisting of the sum of the number density of 
			the buldge, disc and halo.
	'''
	
	
	# Convert l,b to Galactocentric distance and height in cylindrical coordinates
	c_gal = coord.SkyCoord(l=l, b=b, frame='galactic')
	c_galacto = c_gal.transform_to(coord.Galactocentric,  representation_type='cylindrical')
	z, R = c_galacto.z, c_galacto.rho
	# Transform to spherical Galactocentric coordinates
	c_galacto_sph = c_gal.transform_to(coord.Galactocentric,  representation_type='spherical')
	r_G = c_galacto_sph.distance
	
	
	# Constants used
	# Bulge (Binney & Tremaine 2008)
	rho_b0	= 1.722 * u.solMass / u.pc**3
	alpha_b = 1.8
	a_b 	= 1.0 * u.kpc
	r_b 	= 1.9 * u.kpc
	q_b 	= 0.6
	qmin 	= 1e-2 * u.kpc
	# Disc (Bland-Hawthron & Gerhard 2016)
	sig_t	= 970.294 * u.solMass / u.pc**2
	sig_T	= 268.648 * u.solMass / u.pc**2
	R_t		= 2.6 * u.kpc
	R_T		= 2.0 * u.kpc
	z_t		= 0.3 * u.kpc
	z_T		= 0.9 * u.kpc
	# Halo (Kafle et al 2014)
	rho_h0 	= 5.075e-6 * u.solMass / u.pc**-3
	alpha_h1 = 2.4
	alpha_h2 = 4.5
	r_hb	= 17.2 * u.kpc
	r_ht	= 97.7 * u.kpc
	del_h	= 7.1 * u.kpc
	ep_h	= (r_ht/del_h) - 4.5
	r_hmin	= 0.5 * u.kpc
		
	
	# The bulge model as described by Binney & Tremaine (2008)
	q = np.sqrt(R**2 + (z/q_b)**2)
	if q > qmin:
		rho_b = rho_b0 * (q/a_b)**(-alpha_b) * np.exp(-q**2/r_b**2)
	else:
		rho_b = rho_b0 * (qmin/a_b)**(-alpha_b) * np.exp(-q**2/r_b**2)

	
	# The disc model described by Binney & Tremaine (2008) consisting of both the thin and thick disc
	rho_d = (sig_t/(2*z_t))*np.exp(-(R/R_t)-(np.abs(z)/z_t)) + (sig_T/(2*z_T)) * np.exp(-(R/R_T)-(np.abs(z)/z_T))

	
	# The halo model as an isotropic double-power-law function with an exponential decline beyond a certain distance (Kafle et al. 2014)
	for r_G < r_hmin:
		rho_h = rho_h0 * (r_hmin/r_hb)**(-alpha_h1)
	for r_hmin <= r_G < r_hb:
		rho_h = rho_h0 * (r_G/r_hb)**(-alpha_h1)
	for r_hb <= r_G < rht:
		rho_h = rho_h0 * (r_G/r_hb)**(-alpha_h2)
	for r_G >= r_ht:
		rho_h = rho_h0 * (r_ht/r_hb)**(-alpha_h2) * (r_G/r_ht)**ep_h * np.exp(-(r_G - r_ht)/del_h)


	return rho_b + rho_d + rho_h



def Prior(r, l, b)
	'''
	The Milky Way prior as described in Astraatmadja & Bailer-Jones (2016).
	----------------
	Input Parameters
	----------------
	r		Distance r in pc
	l,b 	The location on the sky in Galactic Coordinates
	
	------
	Output
	------
	P		The prior. The probability on the distance r.	
	'''
	return r**2 * rho_MW(r,l,b) * p_obs(r,l,b)


# Probability on r
def P(r, std, varpi):
	'''
	The probability distribution on the distance r, given the parallax and uncertainty on the parallax. Using the prior and the likelihood on the parallax.
	----------------
	Input Parameters
	----------------
	r		Distance in parsec. If not, it will be converged to parsec units
	std		The observed uncertainty on the parallax in arcseconds
	varpi	The observed parallax in arcseconds. 
	
	------
	Output
	------
	P		The probability distribution of the distance	
	'''
	
	# Change to the right units
	r = r.to(u.pc)
	std = std.to(u.arcsec)
	varpi = varpi.to(u.arcsec)	
	
	
	P = 1.0/norm * Likelihood(r, std, varpi) * Prior(r, l, b)
	return P




