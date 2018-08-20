import numpy as np

total = 2143475885
n = [1, int(total/4.), int(total/2.), int((3./4.)*total), total]

for i in range(1, 5):
	random_values = np.random.randint(n[i-1], n[i], int(5e5))
	file_name = 'sample_{}.txt'.format(str(i))
	with open(file_name, 'w') as f:
		for item in random_values:
			f.write("{}\n".format(item))
	f.close


