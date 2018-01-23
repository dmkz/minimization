# Автор: Бураханова А.

import matplotlib.pyplot as plt 
import fileinput 
import re
import os

def rf(f):
	x = []
	y = []    
	line = f.readline()
	while line:
		a, b = line.split(" ")
		x.append(a)
		y.append(b)
		line = f.readline()

	return x,y

def find_files(path, regex):
    ret = []
    pattern = re.compile(regex)
    for root, dir, files in os.walk(path):
        for name in files:
            if pattern.match(name):
                ret.append(name)
    return ret

dir_path = os.path.dirname(os.path.realpath(__file__))
file_list = find_files(dir_path, ".*projection_")
print(file_list)



for file in file_list:
	f = open(file, 'r')
	x, y = rf(f)
	plt.xlim(0,1) 
	plt.ylim(0,1) 
	plt.xlabel('x') 
	plt.ylabel('y')
	plt.plot(x,y, 'rs', markersize = 1)
	plt.savefig(file+'.png')
	plt.clf()
	f.close()



