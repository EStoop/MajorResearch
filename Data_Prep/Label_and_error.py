# -*- coding: utf-8 -*-
import numpy as np
import astropy.units as u
from astropy.io import votable
from astropy.io import fits
import glob
import random
import estimate_distance_v2 as est_dist


# Variables
n = 2
N_STARS = 333333 * n

def label_data(age):
	label = np.zeros_like(age)
	
	# Labels disk = 0, halo = 1, bulge = disk
	ind = np.where( age < 11)[0]
	label[ind] = 0
	ind = np.where(age == 14)[0]
	label[ind] = 1
	
	return label

# ------------------------------------------------------------------------------------------

# Read in data Halo stars
print "Read in halo stars"
file_robin = '/net/boorne/data2/marchetti/GAIA_DR2/NN/halo_mock_catalogue.fits'
hdu_halo = fits.open(file_robin)
Ra       = (hdu_halo[1].data['ra'] * u.rad).to(u.deg)      # deg
Dec      = (hdu_halo[1].data['dec'] * u.rad).to(u.deg)     # deg
Varpi    = hdu_halo[1].data['varpi'] * u.mas   # mas
E_varpi  = hdu_halo[1].data['e_varpi'] * u.mas # mas
Pmra     = hdu_halo[1].data['pmra']    # mas/yr
Pmdec    = hdu_halo[1].data['pmdec']   # mas/yr
#Vrad     = hdu_halo[1].data['vrad']    # km/s
Label    = np.full(len(Ra), 1)

print"Varpi    E_varpi"
print len(Varpi), len(E_varpi)


# Read in hyper velocity stars
print "Read in hyper velocity stars"
data_file2 = "/home/stoop/disks/Major/Dataset/HypVelStars/hyper_velocity_stars.npy"
data2 = np.load(data_file2)
data2 = data2.T
# Choose which ones to delete to get N_STARS of this class
n_hvs = len(data2)
delete = np.random.choice(range(n_hvs), n_hvs - N_STARS, replace=False) # sample of which stars to delete from halo stars
data2 = np.delete(obj=delete, arr=data2, axis = 0)

print "Hypervelocity stars left {}".format(len(data2))
Ra    = np.append(Ra,    (data2[:,0] * u.rad).to(u.deg))  # deg
Dec   = np.append(Dec,   (data2[:,1] * u.rad).to(u.deg)) # deg
Varpi = np.append(Varpi, data2[:,2]) # mas
E_varpi = np.append(E_varpi, np.full(len(data2), 1e-200))
Pmra  = np.append(Pmra,  data2[:,3]) # mas/yr
Pmdec = np.append(Pmdec, data2[:,4]) # mas/yr
Label = np.append(Label, np.full(len(data2), 2))

print len(Varpi), len(E_varpi)

# ------------------------------------------------------------------------------------------

# Read in the rest of the data
file_path = '/home/stoop/disks/Major/Dataset/GUMS/'

for i, file_name in enumerate(glob.glob(file_path + 'gums*')):
	
	# Read in the table
	table = votable.parse_single_table(file_name)
	print 'Has read in table {}'.format(file_name)
	
	age = table.array['Age'].data  # Gyr
	r   = table.array['r'].data    # pc
	G   = table.array['Gmag'].data
	V_I = table.array['V-I'].data
	pop = table.array['Pop'].data

	# Determine the label of the given star
	label_new = label_data(age)
	print 'Data is labelled'

	# Calculate the uncertainties as described in https://www.cosmos.esa.int/web/gaia/science-performance
	zG      = 10**(0.4*(G-15))
	ind     = np.where( 10**(0.4*(12.09-15)) >  zG)[0]
	zG[ind] = 10**(0.4*(12.09-15))
	e_par   = ((-1.631 + 680.766 * zG + 32.732 * zG**2)**(0.5) * (0.986 + (1 - 0.986) * V_I)) * 10**(-6)

	# Create potentially negative parallaxes
	varpi     = (1./r * u.arcsec).to(u.mas) # mas
	e_par     = (e_par* u.arcsec).to(u.mas) # mas
	new_varpi = [np.random.uniform(varpi.value[i] - e_par.value[i], varpi.value[i] + e_par.value[i], 1)[0] for i in range(len(varpi))] # mas
	new_varpi = new_varpi * u.mas
	
	# Combine the data
	Ra       = np.append(Ra,      table.array['_RAJ2000'].data)
	Dec      = np.append(Dec,     table.array['_DEJ2000'].data)
	Varpi    = np.append(Varpi,   new_varpi)
	E_varpi  = np.append(E_varpi, e_par)
	Pmra     = np.append(Pmra,    table.array['pmRA'].data)
	Pmdec    = np.append(Pmdec,   table.array['pmDE'].data)
	#Vrad     = np.append(Vrad,    table.array['RV'].data)
	Label    = np.append(Label,   label_new)
	
	print len(Varpi), len(E_varpi)
	
	print 'Data is combined'

np.save('All_data.npy', np.array([Ra, Dec, Varpi, E_varpi, Pmra, Pmdec, Label]))

#Ra, Dec, Varpi, E_varpi, Pmra, Pmdec, Label = np.load('All_data.npy')

# Remove the stars where the absolute relative error is bigger than 1
delete = np.where(np.abs(E_varpi/Varpi) > 1)[0]
print "Stars with big absolute relative error {}".format(len(delete))
Ra, Dec, Varpi, E_varpi, Pmra, Pmdec, Label = np.delete(arr=Ra, obj=delete, axis=0), np.delete(arr=Dec, obj=delete, axis=0), np.delete(arr=Varpi, obj=delete, axis=0), np.delete(arr=E_varpi, obj=delete, axis=0), np.delete(arr=Pmra, obj=delete, axis=0), np.delete(arr=Pmdec, obj=delete, axis=0), np.delete(arr=Label, obj=delete, axis=0)# delete rows from file

# Randomly sample Robin or gums halo stars to get N_STARS halo stars
ind_halo = np.where(Label == 1)[0]
print "Originally {} halo stars".format(len(ind_halo))
n_delete = len(ind_halo) - N_STARS # number of stars to delete
delete = np.random.choice(ind_halo, n_delete, replace=False) # sample of which stars to delete from halo stars
Ra, Dec, Varpi, E_varpi, Pmra, Pmdec, Label = np.delete(arr=Ra, obj=delete, axis=0), np.delete(arr=Dec, obj=delete, axis=0), np.delete(arr=Varpi, obj=delete, axis=0), np.delete(arr=E_varpi, obj=delete, axis=0), np.delete(arr=Pmra, obj=delete, axis=0), np.delete(arr=Pmdec, obj=delete, axis=0), np.delete(arr=Label, obj=delete, axis=0)# delete rows from file

# Randomly sample disk stars to get N_STARS halo stars
ind_disk = np.where(Label == 0)[0]
print "Originally {} disk stars".format(len(ind_disk))
n_delete = len(ind_disk) - N_STARS
delete = np.random.choice(ind_disk, n_delete, replace=False) # sample of which stars to delete from halo stars
Ra, Dec, Varpi, E_varpi, Pmra, Pmdec, Label = np.delete(arr=Ra, obj=delete, axis=0), np.delete(arr=Dec, obj=delete, axis=0), np.delete(arr=Varpi, obj=delete, axis=0), np.delete(arr=E_varpi, obj=delete, axis=0), np.delete(arr=Pmra, obj=delete, axis=0), np.delete(arr=Pmdec, obj=delete, axis=0), np.delete(arr=Label, obj=delete, axis=0)# delete rows from file

print 'Number of disk stars = {0}\nNumber of halo stars = {1}'.format( len(np.where(Label == 0)[0]), len(np.where(Label == 1)[0]) )

# Save the data
if (len(np.where(Label == 0)[0]) == N_STARS) and ( len(np.where(Label == 1)[0]) == N_STARS):
	
	# Save the data 
	print 'Saving data'
	file_name_out = 'GUMS_labelled_data_{}M.fits'.format(n)
	#c1 = fits.Column(name="ID",      unit = '',       array = )
	c2 = fits.Column(name="Ra",      unit = 'deg',    format = 'E', array = Ra)
	c3 = fits.Column(name="Dec",     unit = 'deg',    format = 'E', array = Dec)
	c4 = fits.Column(name="pmRa",    unit = 'mas/yr', format = 'E', array = Pmra)
	c5 = fits.Column(name="pmDec",   unit = 'mas/yr', format = 'E', array = Pmdec)
	c6 = fits.Column(name="varpi",   unit = 'mas',    format = 'E', array = Varpi)
	c7 = fits.Column(name="label",    unit = '',       format = 'E', array = Label)

	#t = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8, c9])
	t = fits.BinTableHDU.from_columns([c2,c3,c4,c5,c6,c7])
	t.writeto(file_name_out)

