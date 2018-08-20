import numpy as np
from astropy.io import fits

folder_name = '/home/stoop/disks/Major/Dataset/'
def data(normed = False):
	''' Create the data and split it in train, validation and test set 
	in the ratios 80:10:10 percent'''
	
	# Variables
	if normed: file_name = 'Labelled_normed_data_1M.fits'
	else: file_name = 'Labelled_data_1M.fits'
	
	# Read in the data
	hdu = fits.open(folder_name + file_name)
	data = np.array(hdu[1].data.tolist())
	
	# Split data in three groups
	np.random.shuffle(data)
	train, validation, test = np.split(data, [int(.8 * len(data)), int(.9 * len(data))])
	
	train_x, train_y = train[:, 0:-1], train[:,-1]
	validation_x, validation_y = validation[:, 0:-1], validation[:,-1]
	test_x, test_y= test[:, 0:-1], test[:,-1]
	
	# Change data to the right format
	train_y_vect = np.zeros((len(train_y), 3))
	for row, column in enumerate(train_y):
		train_y_vect[row, int(column)] = 1
	
	validation_y_vect = np.zeros((len(validation_y), 3))
	for row, column in enumerate(validation_y):
		validation_y_vect[row, int(column)] = 1
	
	test_y_vect = np.zeros((len(test_y), 3))
	for row, column in enumerate(test_y):
		test_y_vect[row, int(column)] = 1
	
	return train_x, train_y_vect, validation_x, validation_y_vect, test_x, test_y_vect

def write_out(normed = False):
	X_train, Y_train, X_val, Y_val, X_test, Y_test = data(normed)
	if normed:
		np.save(folder_name + 'normed_train_set.npy', X_train)
		np.save(folder_name + 'normed_validation_set.npy', X_val)
		np.save(folder_name + 'normed_test_set.npy', X_test)
	else:
		np.save(folder_name + 'train_set.npy', X_train)
		np.save(folder_name + 'validation_set.npy', X_val)
		np.save(folder_name + 'test_set.npy', X_test)
	
	np.save(folder_name + 'validation_set_label.npy', Y_val)
	np.save(folder_name + 'train_set_label.npy', Y_train)
	np.save(folder_name + 'test_set_label.npy', Y_test)
	return

def main(norm=False):
	write_out(norm)
