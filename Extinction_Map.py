import numpy as np
import astropy.coordinates as coord
import astropy.units as u
import pidly

# Calculate extinction
def Ag(l,b,r):
	'''
	Calculates the extinction in the G band using the extinction map of Drimmel et al. 2003 (ftp://ftp.to.astro.it/astrometria/extinction/ and https://github.com/jkrick/idlbin/tree/master/COBE/quidl) and the ratio Ag/Av = 0.695.
	------
	Input:
	------
	l	galactic longitude in degrees
	b	galactic latitude in degrees
	r	distance in kpc

	------
	Output:
	------
	Ag_rescaled		Extinction in G band with rescaling
	Ag				Extinction without rescaling
	'''
	
	l = l.to(u.deg).value
	b = b.to(u.deg).value
	r = r.to(u.kpc).value
	
	idl = pidly.IDL()
	idl('CD, "/disks/strw9/stoop/Major/Codes/ExtinctionMapFiles/idlbin-master/COBE/quidl/"')
	idl('findext, {0},{1},{2},abs,amod'.format(l,b,r))
	Av_rescaled, Av = idl.abs, idl.amod
	
	return Av_rescaled*0.695, Av*0.695


def Mg_faint(l,b,r, rescale = True):
	'''
	Calculates the faintest magnitude visible in the G band for a certain location and distance.
	------
	Input:
	------
	l	galactic longitude
	b	galactic latitude
	r	distance

	------
	Output:
	------
	Mg_faint	The faintest distance possible to observe in the G band for a certain location and distance. 
	'''
	
	Ag_res, Ag_mod = Ag(l,b,r)
	
	if rescale:
		Mg_faint = 20.0 - 5.0*np.log10(r.to(u.pc).value) + 5 - Ag_res
	else:
		Mg_faint = 20.0 - 5.0*np.log10(r.to(u.pc).value) + 5 - Ag_mod
	
	return Mg_faint




