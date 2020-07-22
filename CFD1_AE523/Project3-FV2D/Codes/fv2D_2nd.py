import numpy as np
from mesh import loadfile
from flux import FluxFunction,FluxBoundary
from post import postprocess, plotstate

'''
def initialstate(parameter):
	rho  = parameter[0]; alfa  = parameter[1]*np.pi/180
	Minf = parameter[2]; gamma = parameter[3]
	u = Minf*np.cos(alfa)
	v = Minf*np.sin(alfa)
	E = 1/((gamma-1)*gamma) + Minf**2/2
	U0 = np.array([rho,rho*u,rho*v,rho*E])
	return U0
'''
def Residual(u,U0,norvecIE,norvecBE,dlIE,dlBE,mesh,midpoint,P):
	v,V,E,IE,BE,area,xcent = loadfile(mesh)
	midpointIE = midpoint[0]
	midpointBE = midpoint[1]
	gradu = [[np.zeros(4) for i in range(2)] for j in range(len(E))]
	# Interior Elements
	for i in range(len(IE)): # Build up dradient u
		#define: state in L-element and R-element
		L = IE[i][2]-1; uL = u[L]
		R = IE[i][3]-1; uR = u[R]
		uhat = 0.5*(uL+uR)	
		n = norvecIE[i]
		# + to gradu[L], - to gradu[R]
		gradu[L][0] += (uhat*dlIE[i])*n[0]
		gradu[L][1] += (uhat*dlIE[i])*n[1]
		gradu[R][0] -= (uhat*dlIE[i])*n[0]
		gradu[R][1] -= (uhat*dlIE[i])*n[1]

	# Boundary Elements	
	for i in range(len(BE)): # Build up dradient u
		n = norvecBE[i]
		L = BE[i][2]-1
		uhat = u[L]
		# + to gradu[L], - to gradu[R]
		gradu[L][0] += (uhat*dlBE[i])*n[0]
		gradu[L][1] += (uhat*dlBE[i])*n[1]
	
	Res = [np.zeros(4) for i in range(len(E))]
	s   = np.zeros(len(E))

	# Interior Elements
	for i in range(len(IE)): # calculate the flux
		#define: state in L-element and R-element
		L = IE[i][2]-1;
		R = IE[i][3]-1;
		u[L] += np.dot(np.array(midpointIE[i])-np.array(xcent[L]),np.divide(gradu[L],area[L]))
		u[R] += np.dot((np.array(midpointIE[i])-np.array(xcent[R])),np.divide(gradu[R],area[R]))
		uL = u[L]; uR = u[R]
		# call n 
		n = norvecIE[i]
		#compute the flux
		flux, smag = FluxFunction(uL,uR,n)
		# +Reidual on L, -Reidual on R
		Ri = np.dot(flux,dlIE[i])	
		Res[L] += Ri
		Res[R] -= Ri
		#+wavespeed on tallies on L and R cell
		c  = smag*dlIE[i]
		s[L] += c
		s[R] += c

	# Boundary Elements
	for i in range(len(BE)): # calculate the flux
		if BE[i][3] == 1:  # Farfield
			L = BE[i][2]-1; n = norvecBE[i]
			u[L] += np.dot((np.array(midpointBE[i])-np.array(xcent[L])),np.divide(gradu[L],area[L]))
			ub = U0; uP = u[L]	
			flux, smag = FluxFunction(uP,ub,n)					
		else:             # Inviscid wall
			L = BE[i][2]-1; n = norvecBE[i]
			u[L] += np.dot((np.array(midpointBE[i])-np.array(xcent[L])),np.divide(gradu[L],area[L]))
			uP = u[L]
			flux, smag, Pb = FluxBoundary(uP,n)
			P[i] = [ BE[i][3], Pb, n, dlBE[i]]  # P[index no.][pressure value][normal vector][length]				
		# +Residual on L (for boundaries only process the left side state)					
		Ri =  np.dot(flux,dlBE[i]) 
		Res[L] += Ri
		#+wavespeed on tallies on L and R cell
		c  = smag*dlBE[i]
		s[L] += c
		
	return Res, s

def secondorderfv(maxiter,mesh,U0,CFL,tol,parameter):
	# load in files
	v,V,E,IE,BE,area,xcent = loadfile(mesh)	
	# Set up and initialize the variables(states)
	Iter = 0;P = [np.zeros(4) for i in range(len(IE))] #np.zeros([len(IE),4])
	u    = np.genfromtxt('u%d.txt' % mesh, delimiter=",")
	U 	 = [0 for i in range(len(E))] # U = updated state 
	uFE  = np.zeros([len(E),4])
	dlIE = np.zeros(len(IE))
	dlBE = np.zeros(len(BE))
	midpointIE = np.zeros([len(IE),2]) # change into np.zeros([len(IE),2])
	midpointBE = np.zeros([len(IE),2])
	strlineIE  = np.zeros([len(IE)]) 
	strlineBE  = np.zeros([len(BE)]) 
	norvecIE   = np.zeros([len(IE),2])
	norvecBE   = np.zeros([len(BE),2])	
	Res_con    = np.zeros(maxiter+1)
	Cl_con     = np.zeros(maxiter+1)		
	# normal vector,dl for IE (norvecIE)	
	for i in range (len(IE)):
		n1 = IE[i][0]; n2 = IE[i][1]
		XA = V[n1][0]; YA = V[n1][1]
		XB = V[n2][0]; YB = V[n2][1]
		dlIE[i] = np.sqrt((XA-XB)**2+(YA-YB)**2)
		norvecIE[i] = [(YB-YA)/dlIE[i],(XA-XB)/dlIE[i]]  
		midpointIE[i] = [0.5*(XA+XB),0.5*(YA+YB)]
	# normal vector, dl, midpoint for BE (norvecBE)	
	for i in range (len(BE)):
		n1 = BE[i][0]; n2 = BE[i][1]
		XA = V[n1][0]; YA = V[n1][1]
		XB = V[n2][0]; YB = V[n2][1]
		dlBE[i] = np.sqrt((XA-XB)**2+(YA-YB)**2)
		norvecBE[i] = [(YB-YA)/dlBE[i],(XA-XB)/dlBE[i]] 
		midpointBE[i] = [0.5*(XA+XB),0.5*(YA+YB)]
	midpoint = [midpointIE,midpointBE]

	# START iteration
	for niter in range(maxiter):
		
		# stage 1
		Res,s = Residual(u,U0,norvecIE,norvecBE,dlIE,dlBE,mesh,midpoint,P)	
		for i in range(len(E)):
			dtoverA = 2*CFL/s[i] #compute the dt
			uFE[i]  = u[i] - dtoverA*Res[i]   #updating

		# stage 2
		Res,s = Residual(uFE,U0,norvecIE,norvecBE,dlIE,dlBE,mesh,midpoint,P)	
		for i in range(len(E)):
			dtoverA = 2*CFL/s[i] #compute the dt
			U[i]  = 0.5*(u[i] + uFE[i] - dtoverA*Res[i])   #updating
			
		u = U
		Iter += 1
				
		# Monitor setting
		if np.mod(Iter,10) == 0:
			Cl_total,Cp,Cl,Cd = postprocess(P,parameter,U0)
			Cl_con[Iter]      = Cl_total
			Res_con[Iter]     = np.linalg.norm(Res,np.inf)
			print('Residual=',np.linalg.norm(Res,np.inf),';','Iter=',Iter,';','Cl_total=',Cl[0],'Cd_total= ',Cd[0])
			
		# Check the tolerance
		if np.linalg.norm(Res,np.inf) <= tol: break

	return u,Cl,Cd,Cp,Cl_con,Res_con,Iter,midpoint,strlineIE,strlineBE
'''
def main():
	mesh = 0; maxiter = 10000; CFL = 0.95; tolerence = 10**-7
	parameter = [1, 5, .25, 1.4, 0.5588]   # [rho, alfa(deg), Minf, gamma, chord]
	U0 = initialstate(parameter)
	u,Cl,Cd,Cp,Cl_con,Res_con,Iter,midpoint,strlineIE,strlineBE = secondorderfv(maxiter,mesh,U0,CFL,tolerence,parameter)
	#plotstate(midpoint,Cp,Res_con,Cl_con,Iter,Cl,Cd,u,parameter,mesh)
	return
	
if __name__ == '__main__':
	main()'''





