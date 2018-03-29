import numpy as np
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u

################################
#      Toombre Diagram
################################
def Toombre(V,UW, over_plot = None):
	
	x = np.arange(-800, 800)
	plt.hist2d(V, UW, bins=500, norm=LogNorm())

	if over_plot != None:
		plt.scatter(extra_stars[0], extra_stars[1], color = 'orange', s=0.8)

	plt.plot(x,Escape(x), c='red')
	plt.plot(x,ThickDisk(x), c = 'red')
	plt.plot(x,ThinDisk(x), c='red')
	plt.xlim(-1000,550)
	plt.ylim(0,1600)
	plt.show()


################################
#    Velocity distributions
################################

def velocity_distr(halo, esc, disk):
# Make a velocity distribution of halo stars
	fig, (ax1, ax2, ax3)= plt.subplots(3,1, sharex=True)

	ax1.hist(np.array(halo[0]), bins=100, alpha=0.5, label = 'V')
	ax1.hist(np.array(halo[1]), bins=100, alpha=0.5, label = 'U')
	ax1.hist(np.array(halo[2]), bins=100, alpha=0.5, label = 'W')
	#ax1.hist(np.array(UWhalo),bins=100, alpha=0.5, label='UW')
	ax1.set_title('Halo stars')
	ax1.legend()

	ax2.hist(np.array(esc[0]), bins=100, alpha=0.5, label = 'V')
	ax2.hist(np.array(esc[1]), bins=100, alpha=0.5, label = 'U')
	ax2.hist(np.array(esc[2]), bins=100, alpha=0.5, label = 'W')
	#ax2.hist(np.array(UWesc),bins=100, alpha=0.5, label='UW')
	ax2.set_title('Escaped stars')

	ax3.hist(np.array(disk[0]), bins=100, alpha=0.5, label = 'V')
	ax3.hist(np.array(disk[1]), bins=100, alpha=0.5, label = 'U')
	ax3.hist(np.array(disk[2]), bins=100, alpha=0.5, label = 'W')
	#ax3.hist(np.array(UWdisk),bins=100, alpha=0.5, label='UW')
	ax3.set_title('Disk stars')

	plt.xlabel('Velocity [km/s]')
	plt.xlim(-550, 550)
	plt.show()




################################
#     Vesc vs. Radius
################################

def escape_velocity(U, V, W, over_plot = None):
	
	r = np.linspace(0.0*u.kpc, 50.0*u.kpc, 100)
	v = v_esc(r)

	gal = c1.transform_to(coord.Galactic())
	vtot = np.sqrt((V+vsun)**2 + U**2 + W**2)
	distance = np.sqrt(gal.cartesian.x**2 + gal.cartesian.y**2 + gal.cartesian.z**2)

	plt.hist2d(distance.to(u.kpc), vtot, bins = 100)
	if over_plot != None:
		plt.scatter(over_plot[0], over_plot[1], s=0.5, color='orange')

	plt.plot(r, v, label=r'$v_{esc}(r)$')
	plt.xlim(0., 25.)
	plt.xlabel('Distance from center [kpc]')
	plt.ylabel('Velociy [km/s]')
	plt.legend()
	plt.show()


################################
#      Make a space image
################################

def space_image(gc1, thin_lim, thick_lim):
	# Make a space image 
	gc1.representation = 'cylindrical'
	R = gc1.rho # distance from galactic centre in pc
	z = gc1.z # height in pc

	print 'Thin disk stars outside boundaries  = {}'.format(len(np.where(np.abs(z.value[thin]) > 330)[0]))
	print 'Thick disk stars outside boundaries = {}'.format(len(np.where(np.abs(z.value[thick]) > 1000)[0]))

	# General image:
	from matplotlib.colors import LogNorm
	plt.hist2d(R, z, bins = 500, norm = LogNorm())
	plt.colorbar()
	plt.axhline(thin_lim, color = 'black', ls='--', label = 'Thin disk limit') # (Chen et al., 2001)
	plt.axhline(-1*thin_lim, color = 'black', ls='--')
	plt.axhline(thick_lim, color = 'black', label = 'Thick disk limit') # (Bovy et al., 2012)
	plt.axhline(-1*thick_lim, color = 'black')
	plt.xlabel('Distance for Galactic Centre [pc]')
	plt.ylabel('Scale Height [pc]')
	plt.title('Space image')
	plt.ylim(-1e4, 1e4)
	plt.show()



# Subplot of all three types 
def subplot_space_image(gc1, thin, thick, halo):
	# Make a space image 
	gc1.representation = 'cylindrical'
	R = gc1.rho # distance from galactic centre in pc
	z = gc1.z # height in pc

	f, ax = plt.subplots(3,1, figsize = (7, 10), sharey = True, sharex = True)
		
	#ax[0].hist2d(R, z, bins = 500, norm = LogNorm(), cmap='gray')
	ax[0].hist2d(R[thin], z[thin], bins =500, norm=LogNorm())
	ax[0].axhline(300, color = 'black', ls='--', lw=1) # (Chen et al., 2001)
	ax[0].axhline(-300, color = 'black', ls='--', lw=1)
	ax[0].axhline(1000, color = 'black', lw=1) # (Bovy et al., 2012)
	ax[0].axhline(-1000, color = 'black', lw=1)
	ax[0].set_title('Thin disk')

	#ax[1].hist2d(R, z, bins = 500, norm = LogNorm(), cmap='gray')
	ax[1].hist2d(R[thick], z[thick], bins =500, norm=LogNorm())
	ax[1].axhline(300, color = 'black', ls='--', lw=1) # (Chen et al., 2001)
	ax[1].axhline(-300, color = 'black', ls='--', lw=1)
	ax[1].axhline(1000, color = 'black', lw=1) # (Bovy et al., 2012)
	ax[1].axhline(-1000, color = 'black', lw=1)
	ax[1].set_title('Thick disk')

	#ax[2].hist2d(R, z, bins = 500, norm = LogNorm(), cmap='gray')
	ax[2].hist2d(R[haloStars], z[haloStars], bins =500, norm=LogNorm())
	ax[2].axhline(300, color = 'black', ls='--', lw=1) # (Chen et al., 2001)
	ax[2].axhline(-300, color = 'black', ls='--', lw=1)
	ax[2].axhline(1000, color = 'black', lw=1) # (Bovy et al., 2012)
	ax[2].axhline(-1000, color = 'black', lw=1)
	ax[2].set_title('Halo')

	ax[2].set_xlabel('Distance for Galactic Centre [pc]')
	ax[1].set_ylabel('Scale Height [pc]')
	ax[0].set_ylim(-1e4, 1e4)
	plt.tight_layout()

	#plt.savefig('/home/stoop/disks/Major/Plots/Space_im_diff_types.png')
	plt.show()

