import numpy as np
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u
from matplotlib.colors import LogNorm
from scipy import interpolate
'''
Units:
------
Distances  	[pc]
Time		[Gyr]	
'''

# Constants
R_sun = 8300. * u.pc # Distance to Galactic Centre from the Sun


class Distribution(object):

	def __init__(self, max_distance = 15.0e3, max_scale_height = 3.e3, max_phi = 360., max_age=14., n_stars = int(1e4), norm = False, counts = False):
		self.n_stars = n_stars
		self.norm = norm
		self.counts = counts
		
		if max_distance != 0.0 and not norm: self.R = np.random.randint(0, max_distance, self.n_stars) * u.pc
		if norm: self.R = max_distance*u.pc

		if max_scale_height != 0.0: self.z = np.random.randint(-1*max_scale_height, max_scale_height, self.n_stars) * u.pc
		if norm: self.z = 27.*u.pc
		if max_scale_height == 0.0: self.z = np.full(self.n_stars, 0.0) * u.pc

		if max_phi != 0.0: self.phi = np.random.randint(0.0, max_phi, self.n_stars) * u.deg
		if norm: self.phi = 0.0*u.deg
		if max_phi == 0.0: self.phi = np.full(self.n_stars, 0.0)* u.deg

		self.coords = coord.SkyCoord(frame = 'galactocentric',representation='cylindrical', rho=self.R, phi=self.phi, z=self.z)
		self.x = self.coords.cartesian.x
		self.y = self.coords.cartesian.y
		
		#if norm: self.age = max_age * u.Gyr
		#else: self.age = np.random.uniform(0, max_age, self.n_stars) * u.Gyr
		self.age = np.random.randint(0, max_age, self.n_stars) * u.Gyr
		
		self.thin_dens = self.ThinDisk()
		self.thick_dens = self.ThickDisk()
		self.bulge_dens = self.Bulge()
		self.halo_dens = self.Halo()

	def Flare(self):
		''' Calculates the flaring of the disk at a given galactocentric coordinate. Unitless. '''
		
		g_flare = 0.545e-6 * u.pc**(-1)
		R_flare = 9500 * u.pc 
		return g_flare.value * (self.R - R_flare).value + 1.
	
	def Rho0_Epsilon_Thin(self):
		'''Returns the local density in [M_sun/pc3] and axis ratio epsilon depending on the age of the stars '''
		
		rho0 = np.zeros(self.n_stars)
		epsilon = np.zeros(self.n_stars)
		
		i1 = np.where(self.age <=0.15 * u.Gyr)[0]
		i2 = np.where((.15*u.Gyr < self.age) & (self.age <= 1.*u.Gyr))[0]
		i3 = np.where((1.*u.Gyr < self.age) & (self.age <= 2.*u.Gyr))[0]
		i4 = np.where((2.*u.Gyr < self.age) & (self.age <= 3.*u.Gyr))[0]
		i5 = np.where((3.*u.Gyr < self.age) & (self.age <= 5.*u.Gyr))[0]
		i6 = np.where((5.*u.Gyr < self.age) & (self.age <= 7.*u.Gyr))[0]
		i7 = np.where((7.*u.Gyr < self.age) & (self.age <= 10.*u.Gyr))[0]
		
		rho0[i1], epsilon[i1] = 4.0e-3, 0.0140
		rho0[i2], epsilon[i2] = 7.9e-3, 0.0268
		rho0[i3], epsilon[i3] = 6.2e-3, 0.0375
		rho0[i4], epsilon[i4] = 4.0e-3, 0.0551
		rho0[i5], epsilon[i5] = 4.0e-3, 0.0696
		rho0[i6], epsilon[i6] = 4.9e-3, 0.0785
		rho0[i7], epsilon[i7] = 6.6e-3, 0.0791
		
		return rho0*0.5 / u.pc**3, epsilon # 0.5 is correction factor to get stars/pc3

	def ThinDisk(self, d0 = 1.):
		'''Returns the density distribution of the thin disk population in [M_sun/pc3] '''
		
		output = np.zeros(self.n_stars)
		k_flare = self.Flare()
		rho0, epsilon = self.Rho0_Epsilon_Thin()
		a = np.sqrt( self.R**2 + (self.z/epsilon)**2 )
		h_rp_y = 5000. * u.pc # h_r+ [pc]
		h_rm_y = 3000. * u.pc # h_r- [pc]
		h_rp_o = 2530. * u.pc # h_r+ [pc]
		h_rm_o = 1320. * u.pc # h_r- [pc]

		if self.n_stars == 1:
			if self.age < 0.15*u.Gyr:
				output = (rho0/d0/k_flare) * ( np.exp(-(a/h_rp_y)**2) - np.exp(-(a/h_rm_y)**2) )
				#if counts: output *= np.pi * (self.R.value**2 - 3000.**2) * self.z.value
				if self.counts: output *= np.pi * (5000.-3000.)**2 * 330.
			else: 
				output = (rho0/d0/k_flare) * ( np.exp(-(0.25+(a/h_rp_o)**2)**0.5) - np.exp(-(0.25 + (a/h_rm_o)**2)**0.5) )
				#if counts: output *= np.pi * (self.R.value**2 - 1320.**2) * self.z.value
				if self.counts: output *= np.pi * (2530.-1320.)**2 * 330.
		else:
			# Stars with ages < 0.15 Gyr
			ind_y = np.where(self.age <= 0.15*u.Gyr)[0]  # indices of young stars
			output[ind_y] = (rho0[ind_y]/d0/k_flare[ind_y]) * ( np.exp(-(a[ind_y]/h_rp_y)**2) - np.exp(-(a[ind_y]/h_rm_y)**2) )

			# Stars with ages > 0.15 Gyr
			ind_o = np.where(self.age > 0.15 * u.Gyr)[0] # stars > 0.15 Gyr
			output[ind_o] = (rho0[ind_o]/d0/k_flare[ind_o]) * ( np.exp(-(0.25+(a[ind_o]/h_rp_o)**2)**0.5) - np.exp(-(0.25 + (a[ind_o]/h_rm_o)**2)**0.5) )
		
			# Multiply by the volume to get the counts. Take the gap between bulge and disk into account
			if self.counts:
				#output[ind_y] *= np.pi * (self.R.value[ind_y]**2 - 3000.**2) * self.z.value[ind_y]
				#output[ind_o] *= np.pi * (self.R.value[ind_o]**2 - 1320.**2) * self.z.value[ind_o]
				output[ind_y] *= np.pi * (5000.-3000.)**2 * 330.
				output[ind_o] *= np.pi * (2530.-1320.)**2 * 330.
			
		return output

	def ThickDisk(self, d0 = 1.):
		''' Returns the density distribution of the thick disk population in [M_sun/pc3] '''
		rho0 = 2.83e-3/u.pc**3  # [stars pc-3]
		output = np.zeros(self.n_stars)
		x_l = 72. * u.pc
		h_R = 4000. * u.pc
		h_z = 1200. * u.pc
		k_flare = self.Flare()
		

		if self.n_stars == 1:
			if np.abs(self.z) <= x_l: 
				output = (rho0/d0/k_flare) * np.exp((self.R - R_sun)/h_R) * (1. - ((1/h_z)/(x_l * (2.+ (x_l/h_z)))) *self.z**2)
				#if counts: output *= np.pi * self.R.value**2 * self.z.value
				if self.counts: output *= np.pi * 4000.**2 * 72.
			else: 
				output = rho0 * np.exp(-(self.R - R_sun)/(h_R)) * ((np.exp(x_l/h_z))/(1+(x_l/(2*h_z)))) *np.exp(-(np.abs(self.z)/h_z))
				#if counts: output *= np.pi * self.R.value**2 * (self.z.value - x_l.value)
				if self.counts: output *= np.pi * 4000.**2 * (1200.* 72.)
		
		else:
			# Stars with a scale height < 72 pc
			ind_less = np.where(np.abs(self.z) <= x_l)[0]
			output[ind_less] = (rho0/d0/k_flare[ind_less]) * np.exp((self.R[ind_less] - R_sun)/h_R) * (1. - ((1/h_z)/(x_l * (2.+ (x_l/h_z)))) * self.z[ind_less]**2)

			# Stars with a scale height > 72 pc
			ind_more = np.where(np.abs(self.z) >  x_l)[0]
			output[ind_more] = rho0 * np.exp(-(self.R[ind_more] - R_sun)/(h_R)) * ((np.exp(x_l/h_z))/(1.+(x_l/(2*h_z)))) * np.exp(-(np.abs(self.z[ind_more])/h_z))
		
			if self.counts:
				#output[ind_less] *= np.pi * self.R.value[ind_less]**2 * self.z.value[ind_less]
				#output[ind_more] *= np.pi * self.R.value[ind_more]**2 * (self.z.value[ind_more] - x_l.value)
				output[ind_less] *= np.pi * 4000.**2 * 72.
				output[ind_more] *= np.pi * 4000.**2 * (1200.* 72.)
			
		return output

	def Halo(self, d0=1.):
		'''Returns the density distribution of the halo population consideringa spheroid in [M_sun/pc3] '''
		rho0 = 2.185e-5 / u.pc**3 # [stars pc-3]
		output = np.zeros(self.n_stars)
		epsilon = 0.76
		a = np.sqrt( self.R**2 + (self.z/epsilon)**2 ) 
		
		if self.n_stars == 1: 
			if a <= 500. * u.pc: output = rho0/d0 * (500.*u.pc/R_sun)**(-2.44)
			else: output = rho0 * (a/R_sun)**(-2.44)
		
		else:
			# Stars with a's smaller than 500 pc
			ind_less = np.where(a <= 500. * u.pc)[0]
			output[ind_less] = rho0/d0 * (500.*u.pc/R_sun)**(-2.44)

			# Stars with a's greater than 500 pc
			ind_more = np.where(a > 500. * u.pc)[0]
			output[ind_more] = rho0 * (a[ind_more]/R_sun)**(-2.44)
			
		if self.counts: 
			#output *= (4./3.) * np.pi * self.R.value**3 
			output *= (4./3.) * np.pi * 2100.**3 
		return output

	def Bulge(self):
		'''Returns the density distribution of the Bulge population considering a balk in [pc-3] '''

		output = np.zeros(self.n_stars)
		N = 13.7*u.pc**(-3) 							# [stars pc-3]
		x0, y0, z0 = 1590.*u.pc, 424.*u.pc, 424.*u.pc 	# normalizing factors [pc]
		R_c = 2540.*u.pc								# [pc]
		rs2 = np.sqrt( ( (self.x/x0)**2 + (self.y/y0)**2 )**2 + (self.z/z0)**4 ) 
		
		if self.n_stars == 1:
			if np.sqrt(self.x**2 + self.y**2) < R_c: output = N * np.exp(-0.5*rs2)
			else: output = N * np.exp(-0.5*rs2.value) * np.exp(-0.5*(np.sqrt(self.x.to(u.kpc)**2 + self.y.to(u.kpc)**2)).value**2)
		
		else:# Stars with a distance < R_c = 2540 pc
			ind_less = np.where(np.sqrt(self.x**2 + self.y**2) < R_c)[0]
			output[ind_less] = N * np.exp(-0.5*rs2[ind_less])
			
			# Stars with a distance > R_c = 2540 pc
			ind_more = np.where(np.sqrt(self.x**2 + self.y**2) > R_c)[0]
			output[ind_more] = N * np.exp(-0.5*rs2.value[ind_more]) * np.exp(-0.5 *(np.sqrt(self.x.to(u.kpc)[ind_more]**2 + self.y.to(u.kpc)[ind_more]**2)).value**2)
		
		if self.counts: 
			#output *= (4./3.) * np.pi * self.x.value * self.y.value * self.z.value
			output *= (4./3.) * np.pi * x0.value * y0.value * z0.value
		
		return output

	def Bin(self,x, y, n_bins = 75):
		'''Returns the binned array of a given array containing of n_bins'''
		
		binned = np.zeros(self.n_stars)
		if x.unit == 'mas': bins_x = np.linspace(x.min(), 1.0*u.mas, n_bins)
		else: bins_x = np.linspace(x.min(), x.max(), n_bins)
		
		for j in range(0,n_bins-1):
			inds = np.where((x > bins_x[j]) & (x <= bins_x[j+1]))[0]
			binned[inds] = np.mean(y[inds])
		return binned
		

	def Varpi(self):
		''' '''
		dist_to_sun = self.coords.icrs.distance
		varpi = ((1./dist_to_sun).value * u.arcsec).to(u.mas)
		return varpi
	
	def Figure(self, xlabel="Distance[pc]", ylabel="Density [pc^-3]", binned=True, varpi=False, xlimit = None):
		''' Makes a figure of the data '''
		
		if varpi: 
			x = self.Varpi()
			xlabel = "Parallax [mas]"
			plt.ylim(1,1e10)
			xlimit = [0,1]
		else: 
			x = np.abs(self.R)
		
		# Reorder the arrays
		inds = np.argsort(x)
		x = np.sort(x)
		y1 = self.thin_dens[inds]
		y2 = self.thick_dens[inds]
		y3 = self.bulge_dens[inds]
		y4 = self.halo_dens[inds]
		
		# Bin the data
		if binned:
			plt.plot(x, self.Bin(x=x, y=y1), label = 'Thin Disk', color = 'blue')
			plt.plot(x, self.Bin(x=x, y=y2), label = 'Thick Disk', color = 'green')
			plt.plot(x, self.Bin(x=x, y=y3), label = 'Bulge', color = '#5FE4F3')
			plt.plot(x, self.Bin(x=x, y=y4), label = 'Halo', color = 'red')
		else:
			plt.plot(x, y1, label = 'Thin Disk', color = 'blue')
			plt.plot(x, y2, label = 'Thick Disk', color = 'green')
			plt.plot(x, y3, label = 'Bulge', color = '#5FE4F3')
			plt.plot(x, y4, label = 'Halo', color = 'red')
		plt.legend()
		plt.yscale('log')
		plt.xlabel(xlabel)
		if self.counts: ylabel = 'Counts' 
		plt.ylabel(ylabel)
		if xlimit is not None: plt.xlim(xlimit[0],xlimit[-1])
		plt.show()

	def Normalise(self):
		self.__init__(max_distance = 8300., max_age = 4.5, norm = True, n_stars = 1, counts = True)
		print "thin dens: ", self.thin_dens
		d0_thin, d0_thick, d0_halo = self.thin_dens, self.thick_dens, self.halo_dens
		print "Thin: {0}\nThick: {1}\nHalo: {2} ".format(d0_thin, d0_thick, d0_halo)
		
		self.__init__()
		self.thin_dens/= d0_thin
		self.thick_dens/= d0_thick
		self.halo_dens/= d0_halo
		return

	def Classify(self):
		'''Classify the stars according to the heighest probability 
		0 = thin disk, 1 = thick disk, 2 = halo, 3 = bulge star'''
		print self.thin_dens[0], self.thick_dens[0], self.halo_dens[0], self.bulge_dens[0]
		densities = np.stack((self.thin_dens, self.thick_dens, self.halo_dens, self.bulge_dens)).T
		print densities.shape
		class_star = np.argmax(densities, axis=1)
		print class_star[0]
		return class_star

if __name__ == '__main__':
	
	a = Distribution()
	a.Figure()
	
	a.__init__(max_distance = 8300., max_age = 4.5, max_scale_height = 27., norm = True, n_stars = 1, counts = True)
	d0_thin, d0_thick, d0_halo = a.thin_dens, a.thick_dens, a.halo_dens
	print "Thin: {0}\nThick: {1}\nHalo: {2} ".format(d0_thin, d0_thick, d0_halo)
	a.__init__()
	a.thin_dens/= d0_thin.value
	a.thick_dens/= d0_thick.value
	a.halo_dens/= d0_halo.value
	
	classes = a.Classify()

'''
inds = np.argsort(a.R)
a.R = np.sort(a.R)
a.thin_dens = a.thin_dens[inds]
a.thick_dens = a.thick_dens[inds]
a.bulge_dens = a.bulge_dens[inds]
a.halo_dens = a.halo_dens[inds]


plt.plot(a.R, a.thin_dens, '.', ms=0.6),plt.show()
plt.plot(a.R, a.thick_dens, '.', ms=0.6),plt.show()
plt.plot(a.R, a.halo_dens, '.', ms=0.6),plt.show()
plt.plot(a.R, a.bulge_dens, '.', ms=0.6),plt.show()
'''
