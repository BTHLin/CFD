import numpy as np
import matplotlib.pyplot as plt

def loadfile(mesh):
	if mesh == 0:
		v  = np.genfromtxt('0V.txt',dtype=float)
		E  = np.genfromtxt('0E.txt',dtype=int)
		BE = np.genfromtxt('0BE.txt',dtype=int)
		IE = np.genfromtxt('0IE.txt',dtype=int)
		V = [[0,0] for j in range(len(v)+1)]
		area = np.zeros(len(E))
		xcent = [np.zeros(2) for i in range(len(E))]
		for i in range(len(v)):
			V[i+1] = [v[i][0],v[i][1]]
		for i in range(len(E)):
		# find the 3 points for every element
			n1 = E[i][0]; n1x = V[n1][0]; n1y = V[n1][1]
			n2 = E[i][1]; n2x = V[n2][0]; n2y = V[n2][1]
			n3 = E[i][2]; n3x = V[n3][0]; n3y = V[n3][1]
		# compute the length of 3 edges on each element
			l1 = np.sqrt((n1x-n2x)**2+(n1y-n2y)**2)
			l2 = np.sqrt((n2x-n3x)**2+(n2y-n3y)**2)
			l3 = np.sqrt((n3x-n1x)**2+(n3y-n1y)**2)
		# compute the area and centroid for each element
			S = 0.5*(l1+l2+l3) # Heron's formula
			area[i] = np.sqrt(S*(S-l1)*(S-l2)*(S-l3))
			xcent[i] = [(n1x+n2x+n3x)/3,(n1y+n2y+n3y)/3]
	elif mesh == 1:
		v  = np.genfromtxt('1V.txt',dtype=float)
		E  = np.genfromtxt('1E.txt',dtype=int)
		BE = np.genfromtxt('1BE.txt',dtype=int)
		IE = np.genfromtxt('1IE.txt',dtype=int)
		V = [[0,0] for j in range(len(v)+1)]
		area = np.zeros(len(E))
		xcent = [np.zeros(2) for i in range(len(E))]
		for i in range(len(v)):
			V[i+1] = [v[i][0],v[i][1]]
		for i in range(len(E)):
		# find the 3 points for every element
			n1 = E[i][0]; n1x = V[n1][0]; n1y = V[n1][1]
			n2 = E[i][1]; n2x = V[n2][0]; n2y = V[n2][1]
			n3 = E[i][2]; n3x = V[n3][0]; n3y = V[n3][1]
		# compute the length of 3 edges on each element
			l1 = np.sqrt((n1x-n2x)**2+(n1y-n2y)**2)
			l2 = np.sqrt((n2x-n3x)**2+(n2y-n3y)**2)
			l3 = np.sqrt((n3x-n1x)**2+(n3y-n1y)**2)
		# compute the area and centroid for each element
			S = 0.5*(l1+l2+l3) # Heron's formula
			area[i] = np.sqrt(S*(S-l1)*(S-l2)*(S-l3))
			xcent[i] = [(n1x+n2x+n3x)/3,(n1y+n2y+n3y)/3]
	elif mesh == 2:
		v  = np.genfromtxt('2V.txt',dtype=float)
		E  = np.genfromtxt('2E.txt',dtype=int)
		BE = np.genfromtxt('2BE.txt',dtype=int)
		IE = np.genfromtxt('2IE.txt',dtype=int)
		V = [[0,0] for j in range(len(v)+1)]
		area = np.zeros(len(E))
		xcent = [np.zeros(2) for i in range(len(E))]
		for i in range(len(v)):
			V[i+1] = [v[i][0],v[i][1]]
		for i in range(len(E)):
		# find the 3 points for every element
			n1 = E[i][0]; n1x = V[n1][0]; n1y = V[n1][1]
			n2 = E[i][1]; n2x = V[n2][0]; n2y = V[n2][1]
			n3 = E[i][2]; n3x = V[n3][0]; n3y = V[n3][1]
		# compute the length of 3 edges on each element
			l1 = np.sqrt((n1x-n2x)**2+(n1y-n2y)**2)
			l2 = np.sqrt((n2x-n3x)**2+(n2y-n3y)**2)
			l3 = np.sqrt((n3x-n1x)**2+(n3y-n1y)**2)
		# compute the area and centroid for each element
			S = 0.5*(l1+l2+l3) # Heron's formula
			area[i] = np.sqrt(S*(S-l1)*(S-l2)*(S-l3))
			xcent[i] = [(n1x+n2x+n3x)/3,(n1y+n2y+n3y)/3]
	elif mesh == 3:
		v  = np.genfromtxt('3V.txt',dtype=float)
		E  = np.genfromtxt('3E.txt',dtype=int)
		BE = np.genfromtxt('3BE.txt',dtype=int)
		IE = np.genfromtxt('3IE.txt',dtype=int)
		V = [[0,0] for j in range(len(v)+1)]
		area  = np.zeros(len(E))
		xcent = [np.zeros(2) for i in range(len(E))]
		for i in range(len(v)):
			V[i+1] = [v[i][0],v[i][1]]
		for i in range(len(E)):
		# find the 3 points for every element
			n1 = E[i][0]; n1x = V[n1][0]; n1y = V[n1][1]
			n2 = E[i][1]; n2x = V[n2][0]; n2y = V[n2][1]
			n3 = E[i][2]; n3x = V[n3][0]; n3y = V[n3][1]
		# compute the length of 3 edges on each element
			l1 = np.sqrt((n1x-n2x)**2+(n1y-n2y)**2)
			l2 = np.sqrt((n2x-n3x)**2+(n2y-n3y)**2)
			l3 = np.sqrt((n3x-n1x)**2+(n3y-n1y)**2)
		# compute the area and centroid for each element
			S = 0.5*(l1+l2+l3) # Heron's formula
			area[i]  = np.sqrt(S*(S-l1)*(S-l2)*(S-l3))
			xcent[i] = [(n1x+n2x+n3x)/3,(n1y+n2y+n3y)/3]
		
	else:
		print('no file found')
	return v,V,E,IE,BE,area,xcent # v for plot, V for fv calculation

