import numpy as np
import sklearn as s
from astropy.io import fits
from astropy.table import Table
from sklearn import preprocessing

folder_name = '/home/stoop/disks/Major/Dataset/'
def read_data(file_name = 'Labelled_data_1M.fits'):
	''' Reads the created files containing the necessary data and splits the data in 
	the train values x'''
	print "read data"
	hdu = fits.open(folder_name + file_name)
	data = np.array(hdu[1].data.tolist())
	return data, hdu

def normalise(data):
	''' Normalises the data true axis 0 '''
	print "Mean: ", np.mean(data, axis=0)
	print "Std:  ", np.std(data, axis=0)
	
	data_norm = (data - np.mean(data, axis=0))/np.std(data, axis=0)
	return np.mean(data, axis=0), np.std(data, axis=0), data_norm

def write_out(mean, std, file_name = 'Mean_std.txt'):
	parameters = ['ra    ', 'dec   ', 'pm_ra ', 'pm_dec', 'varpi ']
	f = open(folder_name + file_name,'w')
	f.write('       Mean     Std\n')
	for i in range(len(mean)):
		f.write('{0} {1} {2}\n'.format(parameters[i], mean[i],std[i]))
	f.close()
	return

def main():
	x,hdul = read_data()
	mu, sig, normalised_data = normalise(x[:,0:5])
	write_out(mu, sig)
	t = Table(np.concatenate((normalised_data, np.array([x[:,5]]).T),axis=1), names = ['ra','dec','pmra','pmdec','parallax','label'])
	t.write(folder_name + 'Labelled_normed_data_1M.fits', format='fits', overwrite=True)
