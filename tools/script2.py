import matplotlib.pyplot as plt 
import fileinput 

f = open("testn.txt", 'r') 
x = [] 
y = [] 
line = f.readline() 
while line: 
	a, b = line.split(" ") 
	x.append(a) 
	y.append(b) 
	line = f.readline() 

plt.plot(x,y, 'rs', markersize = 1) 

f.close() 

plt.xlim(0,1) 
plt.ylim(0,1) 
plt.xlabel('x') 
plt.ylabel('y')
plt.show() 

