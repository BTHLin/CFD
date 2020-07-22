'''
INSTRUCTION:
	This is the driver code for finite volume method. There are only two steps to start, as below.
	Step 1: set up which order you want to run, first(1) or second(2) for the first parameter.
	Step 2: set up which mesh you want to run, only accept integer numbers of 0 ~ 4 for the second parameter.
			mesh0(1149 elements), mesh1(2116 elements), mesh2(4124 elements), mesh3(8031 elements)
	Set 'fast_run = True' can use the pre-store data to plot contours.   
	Set 'LV_plot  = True' can plot the linear-variation Mach and pressure contours.
	Set 'airfoil_plot = True' can show the airfoil for this case.
'''
import time
import numpy as np
import matplotlib.pyplot as plt
from mesh import loadfile
from flux import FluxFunction,FluxBoundary
from fv2D import firstorderfv,secondorderfv
from post import postprocess, plotstate, airfoilplot,fastrun


def initialstate(parameter):
	rho  = parameter[0]; alfa  = parameter[1]*np.pi/180
	Minf = parameter[2]; gamma = parameter[3]
	u = Minf*np.cos(alfa)
	v = Minf*np.sin(alfa)
	E = 1/((gamma-1)*gamma) + Minf**2/2
	U0 = np.array([rho,rho*u,rho*v,rho*E])
	return U0

def main(order,mesh,fast_run,LV_plot,airfoil_plot):

	start = time.clock()
	maxiter = 2000; CFL = 0.9; tolerence = 10**-7
	parameter = [1, 23.45, .125, 1.4, 0.5588]   # [rho, alfa(deg), Minf, gamma, chord]
	U0 = initialstate(parameter)

	if fast_run == True:
		fastrun(order,mesh,LV_plot)
	elif fast_run == False:
		if order == 1:
			u,Cl,Cd,Cp,Cl_con,Res_con,Iter,midpoint,strlineIE,strlineBE = firstorderfv(maxiter,mesh,U0,CFL,tolerence,parameter)
			#np.savetxt('u%d.txt' %mesh, u, delimiter=',')
		elif order == 2:
			u,Cl,Cd,Cp,Cl_con,Res_con,Iter,midpoint,strlineIE,strlineBE = secondorderfv(maxiter,mesh,U0,CFL,tolerence,parameter)
			#np.savetxt('u%d_2nd.txt' %mesh, u, delimiter=',')
		plotstate(u,Cl,Cd,Cp,Cl_con,Res_con,Iter,midpoint,strlineIE,strlineBE,parameter,mesh,start)
	else:
		print('Error! Please see the instruction.')
	
	# plot the case airfoil
	if airfoil_plot == True:
		airfoilplot()

	plt.show()

	return
	
if __name__ == '__main__':
	order, mesh  = 1, 0 # Order: 1 or 2; mesh: 0 ~ 3
	fast_run     = False # Type True or False
	LV_plot      = True # Type True or False
	airfoil_plot = True # Type True or False
	main(order,mesh,fast_run,LV_plot,airfoil_plot)
