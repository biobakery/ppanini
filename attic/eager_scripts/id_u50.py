import os
import sys

if __name__ == '__main__':
	foo1 = open(sys.argv[1]) #to be excluded
	foo2 = open(sys.argv[2]) 
	
	foo1 = foo1.readlines()
	foo2 = foo2.readlines()
	
	foo = open(sys.argv[3], 'w')
	for i in foo2:
		if not i in foo1:
			foo.writelines([i])
	foo.close()	
