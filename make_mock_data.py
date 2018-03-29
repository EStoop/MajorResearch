import cPickle
import numpy as np
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u
import estimate_distance_v2 as est_dist
import make_fig

# Constants used
vsun = 220.0 * u.km/u.s
vesc = 521.0 * u.km/u.s # From Williams et al. 2017
vthick  = 180.0 * u.km/u.s
vthin = 70.0 * u.km/u.s

# Make the boundaries

def Escape(V):
	return np.sqrt(vesc.value**2 - (V+vsun.value)**2)

def ThickDisk(V):
	return np.sqrt(vthick.value**2 - V**2)

def ThinDisk(V):
	return np.sqrt(vthin.value**2 - V**2)

def HaloStars(V, UW, index):
	istars = []
	for i in index:
		if UW[i] < Escape(V[i]):
			if (V[i] > (-1*vthick.value)) and (V[i] < vthick.value) and  (UW[i] > ThickDisk(V[i])):
				istars.append(i)
			if ((V[i] < (-1*vthick.value)) or (V[i] > vthick.value)) and (UW[i] > 0.0):
				istars.append(i)
	return istars

def DiskStars(V, UW, index):
	istars = []
	for i in index:
		if UW[i] < ThickDisk(V[i]):
			if (V[i] > (-1*vthick.value)) and (V[i] < vthick.value):
				istars.append(i)
	return istars

def ThinDiskStars(V, UW, index):
	istars = []
	for i in index:
		if UW[i] < ThinDisk(V[i]):
			if (V[i] > (-1*vthin.value)) and (V[i] < vthin.value):
				istars.append(i)
	return istars

def ThickDiskStars(V, UW, index):
	istars = []
	for i in index:
		if UW[i] < ThickDisk(V[i]):
			if ((-1*vthick.value) < V[i] < (-1*vthin.value)) or (vthin.value < V[i] < vthick.value):
				istars.append(i)
			if ((V[i] > -1*vthin.value) or (V[i] < vthin.value)) and UW[i] > ThinDisk(V[i]):
				istars.append(i)
	return istars

def EscapedStars(V,UW,index):
	istars = []
	for i in index:
		if (V[i] > (-1*vesc.value - vsun.value)) and (V[i] < (vesc.value - vsun.value)) and (UW[i] > Escape(V[i])):
			istars.append(i)
		if ((V[i] < (-1*vesc.value - vsun.value)) or (V[i] > vesc.value - vsun.value)) and (UW[i] > 0.0):
			istars.append(i)
	return istars


def RemoveRows(ind, arr):
	arr_new = []
	for a in arr:
		arr_new.append(np.delete(a, ind))
	return arr_new

def v_esc(r):
	R0 = 8.5*u.kpc
	alpha = 0.37
	v0 = 521.0*u.km/u.s
	
	v = v0 * (r / R0)**(-0.5*alpha)
	return v #km/s
#--------------------------------------------------------------------------------------------------------------------------


# Read the gog sample file
dir_gog = '/net/blekekil/data2/gog/random_samples/'
nSamples = 1837070.
#nSamples = 114816.
print "\nStart reading data..."
data = cPickle.load( open(dir_gog + 'gog_parameters_' + str(int(nSamples)) + '.pkl' ))
edata = cPickle.load( open(dir_gog + 'gog_errors_' + str(int(nSamples)) + '.pkl' ))
print "Data loaded"

# Split the data
ra  = np.array(data['alpha']) 			# rad 
dec = np.array(data['delta']) 			# rad
parallax = np.array(data['parallax']) 	# mas
mu_ra = np.array(data['mu_alpha_star']) # mas/yr - proper motion in ra direction
mu_dec = np.array(data['mu_delta']) 	# mas/yr - proper motion in dec direction
v_r = np.array(data['v_rad']) 			# km/s
e_parallax = np.array(edata['parallax_err']) # mas

# Remove stars with abs(e_parallax / parallax) > 1
ind1 = np.where(np.abs(e_parallax/parallax) > 1.)[0]
[ra, dec, parallax, mu_ra, mu_dec, v_r, e_parallax] = RemoveRows(ind1, [ra, dec, parallax, mu_ra, mu_dec, v_r, e_parallax])

# Remove nans v_r
ind2 = np.where(np.isnan(v_r))[0]
[ra, dec, parallax, mu_ra, mu_dec, v_r, e_parallax] = RemoveRows(ind2, [ra, dec, parallax, mu_ra, mu_dec, v_r, e_parallax])
print "Stars {0} removed".format(len(ind1) + len(ind2))

# Calculate the distance
print "Calculating distances..."
dist = est_dist.main((parallax*u.mas).to(u.arcsec), (e_parallax*u.mas).to(u.arcsec))

# Give all stars an index value
index = range(len(ra))

# Apply LSR
c1 = coord.ICRS(ra = ra*u.rad, dec = dec*u.rad, distance = dist*u.pc, pm_ra_cosdec = mu_ra*u.mas/u.yr, pm_dec = mu_dec*u.mas/u.yr, radial_velocity = v_r*u.km/u.s) 		# ICRS coordinates
gc1 = c1.transform_to(coord.Galactocentric) # transform the coordinates to LSR frame

U, V, W = gc1.v_x, gc1.v_y - vsun, gc1.v_z   # km/s remove the solar movement of V coordinate
UW = np.sqrt(U**2 + W**2) # km/s

# Define stars within boundaries
halo = HaloStars(V.value, UW.value, index)
escaped = EscapedStars(V.value,UW.value,index)
disk = DiskStars(V.value, UW.value, index)
thin = ThinDiskStars(V.value, UW.value, index)
thick = ThickDiskStars(V.value, UW.value, index)


print 'Found {0}/{1} halo stars'.format(len(haloStars), len(ra))
print 'Found {0}/{1} escaped stars'.format(len(escapedStars),len(ra))
print 'Found {0}/{1} disk stars'.format(len(diskStars),len(ra))



make_fig.subplot_space_image(gc1, thin, thick, halo)

