import cPickle
import numpy as np
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u
import estimate_distance_v2 as est_dist
from matplotlib.colors import LogNorm
import make_fig
import seaborn as sns
import astropy.io.fits as fits

# Constants used
vsun = 220.0 * u.km/u.s
vesc = 521.0 * u.km/u.s # From Williams et al. 2017
vthick  = 180.0 * u.km/u.s
vthin = 70.0 * u.km/u.s

# Make the boundaries

def Escape(V):
	return np.sqrt(vesc.value**2 - (V+vsun.value)**2)

def ThickDisk(V, vth = vthick.value):
	return np.sqrt(vth**2 - V**2)

def ThinDisk(V, vth = vthin.value):
	return np.sqrt(vth**2 - V**2)

def HaloStars(V, UW, index, vth=vthick.value):
	istars = []
	for i in index:
		if UW[i] < Escape(V[i]):
			if (V[i] > (-1*vth)) and (V[i] < vth) and  (UW[i] > ThickDisk(V[i], vth)):
				istars.append(i)
			if ((V[i] < (-1*vth)) or (V[i] > vth)) and (UW[i] > 0.0):
				istars.append(i)
	return istars

def DiskStars(V, UW, index):
	istars = []
	for i in index:
		if UW[i] < ThickDisk(V[i]):
			if (V[i] > (-1*vthick.value)) and (V[i] < vthick.value):
				istars.append(i)
	return istars

def ThinDiskStars(V, UW, index, vth = vthin.value):
	istars = []
	for i in index:
		if UW[i] < ThinDisk(V[i], vth):
			if (V[i] > (-1*vth)) and (V[i] < vth):
				istars.append(i)
	return istars

def ThickDiskStars(V, UW, index, vth1 = vthick.value, vth2 = vthin.value):
	istars = []
	for i in index:
		if UW[i] < ThickDisk(V[i], vth1):
			if ((-1*vth1) < V[i] < (-1*vth2)) or (vth2 < V[i] < vth1):
				istars.append(i)
			if ((V[i] > -1*vth2) or (V[i] < vth2)) and UW[i] > ThinDisk(V[i], vth2):
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

'''
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
'''
data_file = "/net/boorne/data2/marchetti/GAIA_DR2/NN/halo_mock_catalogue.fits"
hdu = fits.open(data_file)
data = hdu[1].data[0:100000]

ra  = data['ra']
dec = data['dec']
parallax   = data['varpi']
e_parallax = data['e_varpi']
mu_ra   = data['pmra']
mu_dec   = data['pmdec']
v_r = data['vrad']

# Remove stars with abs(e_parallax / parallax) > 1
ind1 = np.where(np.abs(e_parallax/parallax) > 1.)[0]
[ra, dec, parallax, mu_ra, mu_dec, v_r, e_parallax] = RemoveRows(ind1, [ra, dec, parallax, mu_ra, mu_dec, v_r, e_parallax])

# Remove nans v_r
ind2 = np.where(np.isnan(v_r))[0]
[ra, dec, parallax, mu_ra, mu_dec, v_r, e_parallax] = RemoveRows(ind2, [ra, dec, parallax, mu_ra, mu_dec, v_r, e_parallax])
print "Stars {0} removed".format(len(ind1) + len(ind2))

print "Number of stars {}".format(len(ra))
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


print '\nFound {0}/{1} halo stars'.format(len(halo), len(ra))
print 'Found {0}/{1} escaped stars'.format(len(escaped),len(ra))	
print 'Found {0}/{1} disk stars'.format(len(disk),len(ra))


make_fig.Toombre(V, UW)


#### Test ####
v_thin, v_thick = 50, 140
halo = HaloStars(V.value, UW.value, index, v_thick)
escaped = EscapedStars(V.value,UW.value,index)
disk = DiskStars(V.value, UW.value, index)
thin = ThinDiskStars(V.value, UW.value, index, v_thin)
thick = ThickDiskStars(V.value, UW.value, index, v_thick, v_thin)

make_fig.Toombre(V, UW, v_thin, v_thick)
# make_fig.subplot_space_image(gc1, thin, thick, halo, v_thin, v_thick)


"""
###### Start on the hao/disk distinction ######
for v_thin in range(70-25, 70+25, 5):
	for v_thick in range(180-39, 180+39, 5):
		print '\nThin disc = {0}, Thick disc={1}'.format(v_thin, v_thick)
		halo = HaloStars(V.value, UW.value, index, vth = v_thick)
		escaped = EscapedStars(V.value,UW.value,index)
		thin = ThinDiskStars(V.value, UW.value, index, vth = v_thick)
		thick = ThickDiskStars(V.value, UW.value, index, vth1 = v_thick, vth2 = v_thin)
		make_fig.subplot_space_image(gc1, thin, thick, halo, v_thin, v_thick, savefig=True)
"""

######### Make Relative Error Space image ###############

gc1.representation = 'cylindrical'
R = gc1.rho # distance from galactic centre in pc
z = gc1.z # height in pc

err = e_parallax/parallax
'''
plt.scatter(R,z, c=err, s=0.2, cmap = 'jet', alpha=0.5)
plt.colorbar()
plt.show()
'''
# plt.pcolor(X=R,Y=z,C = err,norm=LogNorm())
# plt.colorbar()
# plt.show()


###### Distribution of disk (z ~ 0.0) #############

ind = np.where((-0.0005*u.pc >= z) & (z <= 0.0005*u.pc))[0]
plt.hist(R.value[ind]/1000., label = 'Total', bins = 100, histtype='step', color = '#000000')
plt.hist(R.value[np.intersect1d(ind, thin)]/1000., label = 'Thin', bins = 100, histtype='step', color = 'blue')
plt.hist(R.value[np.intersect1d(ind,thick)]/1000., label = 'Thick', bins = 100, histtype='step', color = 'orange')
plt.hist(R.value[np.intersect1d(ind,halo)]/1000., label = 'Halo + bulge', bins = 100, histtype='step', color = 'red')
plt.yscale('log')
plt.title("Distance Distribution")
plt.xlabel("Distance [kpc]")
plt.legend()
plt.show()


