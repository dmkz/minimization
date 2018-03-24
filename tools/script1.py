# Автор: Бураханова А.

import matplotlib.pyplot as plt
import fileinput
import re
import os
import string

def rf(f):
	x = []
	y = []    
	line = f.readline()
	while line:
		a, b = line.split(" ")
		print(a, b)
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
file_list = find_files(dir_path, ".*subcube_integration_64_32768")
print(file_list)

i = 0
markers = ['rs','b^','g--']
sizes = [6,5.5,1]
for file in file_list:
	f = open(file, 'r')
	x, y = rf(f)
	s = file.find('_')
	l = file[0:s]
	plt.plot(x,y,markers[i],label = l, markersize = sizes[i])
	i = i + 1
	f.close()

plt.legend(loc= 4)
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('integration_test_result.png')
