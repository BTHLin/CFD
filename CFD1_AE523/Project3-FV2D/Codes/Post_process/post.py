import time
import numpy as np
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from mesh import loadfile
from flux import FluxFunction,FluxBoundary


def postprocess(P,parameter,U0):
	# set up
	Vel = np.linalg.norm([U0[1]/U0[0],U0[2]/U0[0]]) 
	alfa  = parameter[1]*np.pi/180
	chord = parameter[4]
	# sum up the forces for each part of the wing
	F_main = np.zeros(2)
	F_slat = np.zeros(2)
	F_flap = np.zeros(2)
	#P_main = P_slat = P_flap = 0
	Cp = np.zeros(len(P))
	q_inf  = (0.5*U0[0]*Vel**2)
	P_inf  = (parameter[3]-1)*(U0[3]-0.5*U0[0]*Vel**2)
	# P[index no.][pressure value][normal vector][length]

	for i in range(len(P)):
		if P[i][0] == 2: 	# wing main		
			F_main += np.dot(P[i][1],P[i][2])*P[i][3]
			Cp[i] = (P[i][1]-P_inf)/q_inf
		elif P[i][0] == 3:  # wing slat
			F_slat += np.dot(P[i][1],P[i][2])*P[i][3]
			Cp[i] = (P[i][1]-P_inf)/q_inf
		elif P[i][0] == 4:  # wing flap
			F_flap += np.dot(P[i][1],P[i][2])*P[i][3]
			Cp[i] = (P[i][1]-P_inf)/q_inf
	
	# Cl, Cd calculation
	D_main = F_main[1]*np.sin(alfa)+F_main[0]*np.cos(alfa); L_main = F_main[1]*np.cos(alfa)+F_main[0]*np.sin(alfa)
	D_slat = F_slat[1]*np.sin(alfa)+F_slat[0]*np.cos(alfa); L_slat = F_slat[1]*np.cos(alfa)+F_slat[0]*np.sin(alfa)
	D_flap = F_flap[1]*np.sin(alfa)+F_flap[0]*np.cos(alfa); L_flap = F_flap[1]*np.cos(alfa)+F_flap[0]*np.sin(alfa)
	Cd_main = D_main/(q_inf*chord); Cl_main = L_main/(q_inf*chord)     
	Cd_slat = D_slat/(q_inf*chord); Cl_slat = L_slat/(q_inf*chord)
	Cd_flap = D_flap/(q_inf*chord); Cl_flap = L_flap/(q_inf*chord)
	Cd_total = Cd_main+Cd_slat+Cd_flap
	Cl_total = Cl_main+Cl_slat+Cl_flap
	Cl = [Cl_total,Cl_main,Cl_slat,Cl_flap]
	Cd = [Cd_total,Cd_main,Cd_slat,Cd_flap]

	return Cl_total,Cp,Cl,Cd

def plotstate(u,Cl,Cd,Cp,Cl_con,Res_con,Iter,midpoint,strlineIE,strlineBE,parameter,mesh,start):	
	# variabvles set up 
	M_plot = np.zeros(len(u))
	P_plot = np.zeros(len(u))
	gamma = parameter[3]
	midpointBE = midpoint[1]
	v,V,E,IE,BE,area,xcent = loadfile(mesh)
	Vx = v[:,0]; Vy = v[:,1]
	
	n=0
	for i in range(Iter):
		if np.mod(i,10)==1:
			n += 1
	Res_plot = np.zeros(n)
	Cl_plot = np.zeros(n)
	Iter_plot = np.zeros(n)
	n=0
	for i in range(Iter):
		if np.mod(i,10)==0:	
			Res_plot[n] = Res_con[i]
			Cl_plot[n]  = Cl_con[i]
			Iter_plot[n]= n*10
			n += 1

	# Residual vs. Iter
	plt.figure(1)
	plt.semilogy(Iter_plot,Res_plot)
	plt.title('Residual Convergence')
	plt.xlabel('Iteration')
	plt.ylabel('Residual')
	plt.grid()

	# Cl vs. Iter
	plt.figure(2)
	plt.plot(Iter_plot,Cl_plot)
	plt.title('$C_l$ Convergence')
	plt.xlabel('Iteration')
	plt.ylabel('Lift coefficient $C_l$')
	plt.grid()

	# Cp distribution
	f3 = plt.figure(3)
	for i in range(len(Cp)):
		plt.plot(midpointBE[i][0],Cp[i],'bs') 
	plt.xlim(-.1, .7)
	plt.title('Pressure distribution $C_p$')
	plt.xlabel('location $x/c$')
	plt.ylabel('Pressure coefficient $C_p$')
	plt.gca().invert_yaxis()
	plt.grid()

	# Mach contour
	f4 = plt.figure(4)
	elements=np.zeros([len(E),3])
	for i in range(len(E)):
		elements[i] = [E[i][0]-1,E[i][1]-1,E[i][2]-1]
		M_plot[i] = np.sqrt(u[i][1]**2+u[i][2]**2) 
	#triangles = tri.Triangulation(Vx,Vy,elements)
	plt.tripcolor(Vx, Vy, elements, facecolors=M_plot, edgecolors='k')
	plt.xlim(-.3, .9);	plt.ylim(-.35, .25)
	plt.title('Mach contour with mesh elements: %d' %len(elements))
	plt.xlabel('x direction')
	plt.ylabel('y direction')
	plt.gca().set_aspect('equal')
	plt.colorbar().ax.set_ylabel('Mach number')

	# Pressure contour
	f5 = plt.figure(5)
	elements=np.zeros([len(E),3])
	for i in range(len(E)):
		elements[i] = [E[i][0]-1,E[i][1]-1,E[i][2]-1]
		Vel = np.linalg.norm([u[i][1]/u[i][0],u[i][2]/u[i][0]])
		P_plot[i] = (gamma-1)*(u[i][3]-0.5*u[i][0]*Vel**2)#  u[i][1]/np.cos(parameter[1]*np.pi/180) 
	plt.tripcolor(Vx, Vy, elements, facecolors=P_plot, edgecolors='k')
	plt.xlim(-.3, .9); plt.ylim(-.35, .25)
	plt.title('Pressure contour with mesh elements: %d ' %len(elements))
	plt.xlabel('x direction')
	plt.ylabel('y direction')
	plt.gca().set_aspect('equal')
	plt.colorbar().ax.set_ylabel('Pascal')

	# Streamline plotting
	#    psi = [t(1)/f(0), psi] 
	#    strline = [flux*dl]
	f6 = plt.figure(6)	
	psi = np.zeros([len(v),2])
	x = np.zeros(len(BE))
	y = np.zeros(len(BE))	
	farf_num = 0; main_num = 0
	slat_num = 0; flap_num = 0 

	for i in range(len(BE)):
		if BE[i][3] == 1:
			farf_num += 1
		elif BE[i][3] == 2:
			main_num += 1
		elif BE[i][3] == 3:
			slat_num += 1
		elif BE[i][3] == 4:
			flap_num += 1
	main = np.zeros([main_num,2])
	slat = np.zeros([slat_num,2])
	flap = np.zeros([flap_num,2])
	m = 0; n = 0; l = 0

	for i in range(len(BE)):
		x1 = Vx[BE[i][0]-1]
		y1 = Vy[BE[i][0]-1]
		x2 = Vx[BE[i][1]-1]
		y2 = Vy[BE[i][1]-1] 
		if BE[i][3]==2:
			psi[BE[i][0]-1][0] = 1
			psi[BE[i][1]-1][0] = 1			
			main[n] = [x1,y1]
			plt.plot([x1,x2],[y1,y2],'k')
			m += 1
		elif BE[i][3]==3:			
			slat[n] = [x1,y1]
			plt.plot([x1,x2],[y1,y2],'k')
			n += 1
		elif BE[i][3]==4:			
			flap[l] = [x1,y1]
			plt.plot([x1,x2],[y1,y2],'k')
			l += 1
		
	while 1:
		for i in range(len(IE)):
			node1 = IE[i][0]-1
			node2 = IE[i][1]-1
			if psi[node1][0] == 0 and psi[node2][0] == 1:
				psi[node1][1] = psi[node2][1]-strlineIE[i]
				psi[node1][0] = 1 #mark as calculated				
			elif psi[node1][0] == 1 and psi[node2][0] == 0:
				psi[node2][1] = psi[node1][1]+strlineIE[i]
				psi[node2][0] = 1 #mark as calculated				

		for i in range(len(BE)):
			node1 = BE[i][0]-1
			node2 = BE[i][1]-1
			if psi[node1][0] == 0 and psi[node2][0] == 1:
				psi[node1][1] = psi[node2][1]-strlineBE[i]
				psi[node1][0] = 1 #mark as calculated				
			elif psi[node1][0] == 1 and psi[node2][0] == 0:
				psi[node2][1] = psi[node1][1]+strlineBE[i]
				psi[node2][0] = 1 #mark as calculated			

		if psi[:,0].all() != 0:
			break
	
	for i in range(main_num):
		for j in range(slat_num):
			dist = [i,j,np.linalg.norm(main[i]-slat[j])]
			if j>0:
				if dist[2]<mini:
					mini = dist[2]
					dist_main_slat = dist
				else:
					mini = mini
			else:
				mini = dist[2]
				dist_main_slat = dist 
		node1 = BE[farf_num+dist_main_slat[0]][0]         # node on main wing
		node2 = BE[farf_num+main_num+dist_main_slat[1]][0] # node on slat
		mfr_main_slat = abs(psi[node1][1]-psi[node2][1])
		
		for j in range(flap_num):
			dist = [i,j,np.linalg.norm(main[i]-flap[j])]
			if j>0:
				if dist[2]<mini:
					mini = dist[2]
					dist_main_flap = dist
				else:
					mini = mini
			else:
				mini = dist[2]
				dist_main_flap = dist 
		node1 = BE[farf_num+dist_main_flap[0]][0]         # node on main wing
		node2 = BE[farf_num+main_num+slat_num+dist_main_flap[1]][0] # node on slat
		mfr_main_flap = abs(psi[node1][1]-psi[node2][1])

	triangles = tri.Triangulation(Vx,Vy,elements)
	plt.tricontour(triangles,psi[:,1],levels=np.linspace(-.1,.1,50))
	plt.xlim(-.3, .9);	plt.ylim(-.35, .25)
	plt.title('Streamline contour with mesh elements: %d ' %len(elements))
	plt.xlabel('x direction')
	plt.ylabel('y direction')
	plt.gca().set_aspect('equal')
	 
	# Print and Save data to files 
	end = time.clock()
	print((end-start)/3600,'hours')
	print('Cl_total=',Cl[0],'Cl_main=',Cl[1],'Cl_slat=',Cl[2],'Cl_flap=',Cl[3])
	print('Cd_total=',Cd[0],'Cd_main=',Cd[1],'Cd_slat=',Cd[2],'Cd_flap=',Cd[3])
	print('mass flow rate between main wing and slat:',mfr_main_slat)
	print('mass flow rate between main wing and flap:',mfr_main_flap)

	
	file = open('mesh%d.txt' %mesh, 'w')
	file.write('Cl_total=%f; Cl_main=%f; Cl_slat=%f; Cl_flap=%f \n' % (Cl[0],Cl[1],Cl[2],Cl[3]) )
	file.write('Cd_total=%f; Cd_main=%f; Cd_slat=%f; Cd_flap=%f \n' % (Cd[0],Cd[1],Cd[2],Cd[3]) )
	file.write('mass flow rate between main wing and slat:%f \n' % (mfr_main_slat) )
	file.write('mass flow rate between main wing and flap:%f \n' % (mfr_main_flap) )
	file.write('calculation time:%f hours' % ((end-start)/3600) )
	file.write('Iterations:%d ' % Iter )
	file.close()
	
	# show the graphs

	return

def fastrun(order,mesh,LV_plot):
	# load files
	v,V,E,IE,BE,area,xcent = loadfile(mesh)
	Vx = v[:,0]; Vy = v[:,1]; gamma = 1.4
	if order == 1:		u  = np.genfromtxt('u%d.txt' % mesh, delimiter=",")
	elif order == 2:	u  = np.genfromtxt('u%d_2nd.txt' % mesh, delimiter=",")

	if LV_plot == True:
		newu = np.zeros([len(E),4])#u.copy()
		node = np.zeros([len(v),4])
		cont = np.zeros(len(v))
		elements = np.zeros([len(E),3])
		M_plot   = np.zeros(len(v))
		P_plot   = M_plot.copy()
		for i in range(len(E)):
			n1 = E[i][0]-1
			n2 = E[i][1]-1
			n3 = E[i][2]-1
			node[n1] += u[i]
			cont[n1] += 1
			node[n2] += u[i]
			cont[n2] += 1
			node[n3] += u[i]
			cont[n3] += 1
			elements[i] = [E[i][0]-1,E[i][1]-1,E[i][2]-1]
		triangles = tri.Triangulation(Vx,Vy,elements)
		for i in range(len(v)):
			if cont[i] != 0:
				node[i] = node[i]/cont[i]
				Vel = np.linalg.norm([node[i][1]/node[i][0],node[i][2]/node[i][0]])
				M_plot[i] = np.sqrt(node[i][1]**2+node[i][2]**2)
				P_plot[i] = (gamma-1)*(node[i][3]-0.5*node[i][0]*Vel**2)
			else:
				node[i] = node[i]
		# Mach contour
		f1 = plt.figure(1)
		for i in range(len(BE)):
			x1 = Vx[BE[i][0]-1]
			y1 = Vy[BE[i][0]-1]
			x2 = Vx[BE[i][1]-1]
			y2 = Vy[BE[i][1]-1] 
			if BE[i][3]==2:
				plt.plot([x1,x2],[y1,y2],'k')
			elif BE[i][3]==3:			
				plt.plot([x1,x2],[y1,y2],'k')
			elif BE[i][3]==4:			
				plt.plot([x1,x2],[y1,y2],'k')
		plt.tricontourf(triangles,M_plot,levels=np.linspace(0,.56,1000))
		plt.xlim(-.3, .9);	plt.ylim(-.35, .25)
		plt.title('Mach contour with mesh elements: %d ' %len(elements))
		plt.xlabel('x direction')
		plt.ylabel('y direction')
		plt.gca().set_aspect('equal')
		plt.colorbar().ax.set_ylabel('Mach')
		# Pressure contour
		f2 = plt.figure(2)
		for i in range(len(BE)):
			x1 = Vx[BE[i][0]-1]
			y1 = Vy[BE[i][0]-1]
			x2 = Vx[BE[i][1]-1]
			y2 = Vy[BE[i][1]-1] 
			if BE[i][3]==2:
				plt.plot([x1,x2],[y1,y2],'k')
			elif BE[i][3]==3:			
				plt.plot([x1,x2],[y1,y2],'k')
			elif BE[i][3]==4:			
				plt.plot([x1,x2],[y1,y2],'k')
		plt.tricontourf(triangles,P_plot,levels=np.linspace(.5,.75,1000))
		plt.xlim(-.3, .9);	plt.ylim(-.35, .25)
		plt.title('Pressure contour with mesh elements: %d ' %len(elements))
		plt.xlabel('x direction')
		plt.ylabel('y direction')
		plt.gca().set_aspect('equal')
		plt.colorbar().ax.set_ylabel('Pascal')

	elif LV_plot == False:
		M_plot = np.zeros(len(u))
		P_plot = np.zeros(len(u))
		elements=np.zeros([len(E),3])
		# Mach contour
		f1 = plt.figure(1)
		for i in range(len(E)):
			elements[i] = [E[i][0]-1,E[i][1]-1,E[i][2]-1]
			M_plot[i] = np.sqrt(u[i][1]**2+u[i][2]**2)
		triangles = tri.Triangulation(Vx,Vy,elements)
		plt.tripcolor(triangles, facecolors=M_plot, edgecolors='k')
		plt.xlim(-.3, .9);	plt.ylim(-.35, .25)
		plt.title('Mach contour with mesh elements: %d' %len(elements))
		plt.xlabel('x direction')
		plt.ylabel('y direction')
		plt.gca().set_aspect('equal')
		plt.colorbar().ax.set_ylabel('Mach number')
		# Pressure contour
		f2 = plt.figure(2)
		for i in range(len(E)):
			Vel = np.linalg.norm([u[i][1]/u[i][0],u[i][2]/u[i][0]])
			P_plot[i] = (gamma-1)*(u[i][3]-0.5*u[i][0]*Vel**2)
		plt.tripcolor(triangles, facecolors=P_plot, edgecolors='k')
		plt.xlim(-.3, .9); plt.ylim(-.35, .25)
		plt.title('Pressure contour with mesh elements: %d ' %len(elements))
		plt.xlabel('x direction')
		plt.ylabel('y direction')
		plt.gca().set_aspect('equal')
		plt.colorbar().ax.set_ylabel('Pascal')

	else:
		print('Error! Please see the instruction.')
	return

def airfoilplot():
	f7 = plt.figure(7)
	v,V,E,IE,BE,area,xcent = loadfile(3)
	Vx = v[:,0]
	Vy = v[:,1]
	for i in range(len(BE)):
		x1 = Vx[BE[i][0]-1]
		y1 = Vy[BE[i][0]-1]
		x2 = Vx[BE[i][1]-1]
		y2 = Vy[BE[i][1]-1] 
		if BE[i][3]==2:
			plt.plot([x1,x2],[y1,y2],'b')
		elif BE[i][3]==3:			
			plt.plot([x1,x2],[y1,y2],'b')
		elif BE[i][3]==4:			
			plt.plot([x1,x2],[y1,y2],'b')

	plt.xlim(-.3, .9);	plt.ylim(-.35, .25)
	plt.gca().set_aspect('equal')
	plt.title('Case airfoil geometry')
	plt.xlabel('x direction')
	plt.ylabel('y direction')
	plt.grid()
	return