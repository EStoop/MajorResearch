import numpy as np
 
ra_mu,    ra_std    = 3.67819310709,   1.63856335551
dec_mu,   dec_std   = -0.172488118382, 0.65928924887
varpi_mu, varpi_std = 2.34303642327,   13.9298059115
pmra_mu,  pmra_std  = -2.65923048221,  4622.001842
pmdec_mu, pmdec_std = 236.687514308,   3231.67911015

folder_name = '/home/stoop/disks/Major/Dataset/Old_set/'
files_data  = ['CrossVal_Set.npy', 'Training_Set.npy', 'Test_Set.npy']
files_label = ['CrossVal_Set_output.npy', 'Training_Set_output.npy', 'Test_Set_output.npy']

ra,dec,varpi,pmra,pmdec = np.array([]),np.array([]),np.array([]),np.array([]),np.array([])

for i in range(3):
	x = np.load(folder_name + files_data[i])
	y = np.load(folder_name + files_label[i])
	
	ra     = np.append(ra,    (x[:, 1][y==1] * ra_std)    + ra_mu)      #rad
	dec    = np.append(dec,   (x[:, 2][y==1] * dec_std)   + dec_mu)     #rad
	varpi  = np.append(varpi, (x[:, 3][y==1] * varpi_std) + varpi_mu)   #mas
	pmra   = np.append(pmra,  (x[:, 4][y==1] * pmra_std)  + pmra_mu)    #mas/yr
	pmdec  = np.append(pmdec, (x[:, 5][y==1] * pmdec_std) + pmdec_mu)   #mas/yr

np.save('hyper_velocity_stars.npy', [ra,dec,varpi, pmra, pmdec])
