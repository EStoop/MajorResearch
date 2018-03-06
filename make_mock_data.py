import cPickle
import numpy as np
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u

# Constants used
vsun = 220.0 * u.km/u.s
vesc = 521.0 * u.km/u.s # From Williams et al. 2017
vth  = 180.0 * u.km/u.s
vthin = 70.0 * u.km/u.s
# Make the boundaries

def Escape(V):
	return np.sqrt(vesc.value**2 - (V+vsun.value)**2)

def ThickDisk(V):
	return np.sqrt(vth.value**2 - V**2)

def ThinDisk(V):
	return np.sqrt(vthin.value**2 - V**2)

def HaloStars(V, UW, index):
	istars = []
	for i in index:
		if UW[i] < Escape(V[i]):
			if (V[i] > (-1*vth.value)) and (V[i] < vth.value) and  (UW[i] > ThickDisk(V[i])):
				istars.append(i)
			if ((V[i] < (-1*vth.value)) or (V[i] > vth.value)) and (UW[i] > 0.0):
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
#--------------------------------------------------------------------------------------------------------------------------


# Read the gog sample file
dir_gog = '/net/blekekil/data2/gog/random_samples/'
nSamples = 1837070.
#nSamples = 114816.
print "\nStart reading data"
data = cPickle.load( open(dir_gog + 'gog_parameters_' + str(int(nSamples)) + '.pkl' ))
edata = cPickle.load( open(dir_gog + 'gog_errors_' + str(int(nSamples)) + '.pkl' ))
print "Data loaded"

# Split the data
ra  = np.array(data['alpha']) # rad 
dec = np.array(data['delta']) # rad
parallax = np.array(data['parallax']) # mas
mu_ra = np.array(data['mu_alpha_star']) # mas/yr - proper motion in ra direction
mu_dec = np.array(data['mu_delta']) # mas/yr - proper motion in dec direction
v_r = np.array(data['v_rad']) # km/s

e_parallax = np.array(edata['parallax_err']) # mas


# Remove stars with abs(e_parallax / parallax) > 1
ind1 = np.where(np.abs(e_parallax/parallax) > 1.)[0]
[ra, dec, parallax, mu_ra, mu_dec, v_r, e_parallax] = RemoveRows(ind1, [ra, dec, parallax, mu_ra, mu_dec, v_r, e_parallax])

# Remove nans v_r
ind2 = np.where(np.isnan(v_r))[0]
[ra, dec, parallax, mu_ra, mu_dec, v_r, e_parallax] = RemoveRows(ind2, [ra, dec, parallax, mu_ra, mu_dec, v_r, e_parallax])

print "Stars {0} removed".format(len(ind1) + len(ind2))

# Give all stars an index value
index = range(len(ra))
# Apply LSR
c1 = coord.ICRS(ra = ra * u.rad, dec = dec * u.rad, distance = (parallax* u.mas).to(u.pc, u.parallax()), pm_ra_cosdec = mu_ra * u.mas/u.yr, pm_dec = mu_dec * u.mas/u.yr, radial_velocity = v_r * u.km/u.s) # ICRS coordinates
gc1 = c1.transform_to(coord.Galactocentric) # transform the coordinates to LSR frame
U, V, W = gc1.v_x, gc1.v_y - vsun, gc1.v_z   # km/s remove the solar movement of V coordinate
UW = np.sqrt(U**2 + W**2) # km/s

# Plot the stars and boundaries
x = np.arange(-800, 800)
plt.scatter(V, UW, s=0.1, color = 'blue')
plt.plot(x,Escape(x), c='red')
plt.plot(x,ThickDisk(x), c = 'red')
plt.plot(x,ThinDisk(x), c='red')
plt.xlim(-1000,550)
plt.ylim(0,1600)
plt.show()

# Define stars within boundaries
haloStars = HaloStars(V.value, UW.value, index)
print 'Found {0}/{1} halo stars'.format(len(haloStars), len(ra))
escapedStars = EscapedStars(V.value,UW.value,index)
print 'Found {0}/{1} escaped stars'.format(len(escapedStars),len(ra))

# Relative error on parallax
err = e_parallax[haloStars]/parallax[haloStars]
ind_err = np.where((err<0.2) & (err > 0.0))[0]

err_esc = e_parallax[escapedStars]/parallax[escapedStars]
ind_err_esc = np.where((err_esc < 0.2) & (err_esc > 0.0))[0] 

err_all = e_parallax/parallax
ind_err_all = np.where((err_all<0.2) & (err_all >0.0))[0]

# Plot the Halo Stars 
plt.scatter(V, UW, s=0.2, alpha=0.5, label= 'Relative error > 0.2')
plt.scatter(np.array(V)[ind_err_all], np.array(UW)[ind_err_all], color = 'orange', s=0.8, label= 'Relative error < 0.2')
plt.plot(x,Escape(x), c='red')
plt.plot(x,ThickDisk(x), c = 'red')
plt.plot(x,ThinDisk(x), c='red')
plt.xlim(-1000,1000)
plt.ylim(0,1600)
plt.legend()
plt.show()

'''
# Plot the Halo Stars
plt.scatter(V, UW, s=0.2, alpha=0.5)
plt.scatter(np.array(V)[haloStars], np.array(UW)[haloStars], color = 'green', s=0.5, alpha=0.5)
plt.scatter(np.array(V)[escapedStars], np.array(UW)[escapedStars], color = 'red', s=0.5, alpha=0.5)
plt.plot(x,Escape(x), c='red')
plt.plot(x,ThickDisk(x), c = 'red')
plt.plot(x,ThinDisk(x), c='red')
plt.xlim(-1000,1000)
plt.ylim(0,1600)
plt.show()
'''

# Make a velocity distribution of halo stars
plt.hist(np.array(V)[haloStars], bins=100, alpha=0.5, label = 'V')
plt.hist(np.array(U)[haloStars],bins=100, alpha=0.5, label = 'U')
plt.hist(np.array(W)[haloStars],bins=100, alpha=0.5, label = 'W')
plt.hist(np.array(UW)[haloStars],bins=100, alpha=0.5, label='UW')
plt.legend()
plt.show()

print 'Statistics halo stars:\nmu, std'
print 'V:', np.nanmean(V[haloStars]), np.nanstd(V[haloStars])
print 'U:', np.nanmean(U[haloStars]), np.nanstd(U[haloStars])
print 'W:',np.nanmean(W[haloStars]), np.nanstd(W[haloStars])

# Make error distribution of halo stars on parallax
plt.suptitle(r'Relative error distribution halo stars ')
plt.title('$\mu=${0:.3f} , $\sigma=${1:.3f}'.format(np.mean(err), np.std(err)))
plt.xlabel(r'$\sigma_{\varpi}/\varpi$')
plt.hist(err, bins=100, alpha=0.5)
plt.axvline(x=np.mean(err),color = 'red')
plt.show()

# Make error distribution of escaped stars on parallax
plt.suptitle(r'Relative error distribution stars with escape velocities')
plt.title('$\mu=${0:.3f} , $\sigma=${1:.3f}'.format(np.mean(err_esc), np.std(err_esc)))
plt.xlabel(r'$\sigma_{\varpi}/\varpi$')
plt.hist(err_esc, bins=50, alpha=0.5)
plt.axvline(x=np.mean(err_esc),color = 'red')
plt.show()


