import numpy as np

#-----------------------------------------------------------
def FluxFunction(UL, UR, n):
# PURPOSE: This function calculates the flux for the Euler equations
# using the Roe flux function
#
# INPUTS:
#    UL: conservative state vector in left cell
#    UR: conservative state vector in right cell
#     n: normal pointing from the left cell to the right cell
#
# OUTPUTS:
#  F   : the flux out of the left cell (into the right cell)
#  smag: the maximum propagation speed of disturbance
    
    gamma = 1.4
    gmi = gamma-1.0

    # process left state
    rL = UL[0]
    uL = UL[1]/rL
    vL = UL[2]/rL
    unL = uL*n[0] + vL*n[1]
    qL = np.sqrt(UL[1]**2 + UL[2]**2)/rL
    pL = (gamma-1)*(UL[3] - 0.5*rL*qL**2)
    if ((pL<0) | (rL<0)): print('Non-physical state!') 
    rHL = UL[3] + pL
    HL = rHL/rL
    cL = np.sqrt(gamma*pL/rL)

    # left flux
    FL = np.zeros(4)
    FL[0] = rL*unL
    FL[1] = UL[1]*unL + pL*n[0]
    FL[2] = UL[2]*unL + pL*n[1]
    FL[3] = rHL*unL

    # process right state
    rR = UR[0]
    uR = UR[1]/rR
    vR = UR[2]/rR
    unR = uR*n[0] + vR*n[1]
    qR = np.sqrt(UR[1]**2 + UR[2]**2)/rR
    pR = (gamma-1)*(UR[3] - 0.5*rR*qR**2)
    if ((pR<0) | (rR<0)): print('Non-physical state!') 
    rHR = UR[3] + pR
    HR = rHR/rR
    cR = np.sqrt(gamma*pR/rR)

    # right flux
    FR = np.zeros(4)
    FR[0] = rR*unR
    FR[1] = UR[1]*unR + pR*n[0]
    FR[2] = UR[2]*unR + pR*n[1]
    FR[3] = rHR*unR

    # difference in states
    du = UR - UL

    # Roe average
    di     = np.sqrt(rR/rL)
    d1     = 1.0/(1.0+di)
    
    ui     = (di*uR + uL)*d1
    vi     = (di*vR + vL)*d1
    Hi     = (di*HR + HL)*d1

    af     = 0.5*(ui*ui+vi*vi )
    ucp    = ui*n[0] + vi*n[1]
    c2     = gmi*(Hi - af)
    if (c2<0): print('Non-physical state!') 
    ci     = np.sqrt(c2)
    ci1    = 1.0/ci
    
    # eigenvalues
    l = np.zeros(3)
    l[0] = ucp+ci; l[1] = ucp-ci; l[2] = ucp
    
    # entropy fix
    epsilon = ci*.1
    for i in range(3):
        if ((l[i]<epsilon) & (l[i]>-epsilon)):
            l[i] = 0.5*(epsilon + l[i]*l[i]/epsilon)

    l = abs(l); l3 = l[2]

    # average and half-difference of 1st and 2nd eigs
    s1    = 0.5*(l[0] + l[1])
    s2    = 0.5*(l[0] - l[1])

    # left eigenvector product generators (see Theory guide)
    G1    = gmi*(af*du[0] - ui*du[1] - vi*du[2] + du[3])
    G2    = -ucp*du[0]+du[1]*n[0]+du[2]*n[1]

    # required functions of G1 and G2 (again, see Theory guide)
    C1    = G1*(s1-l3)*ci1*ci1 + G2*s2*ci1
    C2    = G1*s2*ci1          + G2*(s1-l3)

    # flux assembly
    F = np.zeros(4)
    F[0]    = 0.5*(FL[0]+FR[0]) - 0.5* (l3*du[0] + C1)
    F[1]    = 0.5*(FL[1]+FR[1]) - 0.5* (l3*du[1] + C1*ui + C2*n[0])
    F[2]    = 0.5*(FL[2]+FR[2]) - 0.5* (l3*du[2] + C1*vi + C2*n[1])
    F[3]    = 0.5*(FL[3]+FR[3]) - 0.5* (l3*du[3] + C1*Hi + C2*ucp )
    
    # max wave speed
    smag = max(l)

    return F, smag
    

def FluxBoundary(U,n):
    gamma = 1.4
    gmi   = gamma-1.

    # Process the state
    r  = U[0]
    u  = U[1]
    v  = U[2]
    rE = U[3]
    V = np.array([u,v])
    n = np.array([n[0],n[1]]) 
    # Boundary velocity
    Vb  = V - np.dot(V,n)*n    
    mag_Vb = np.linalg.norm(Vb)

    # Boundary pressure 
    Pb = gmi*(rE-0.5*r*mag_Vb**2)
    
    # Maximum speed 
    c = np.sqrt(gamma*Pb/r)
    smag = np.abs(np.dot(V,n)) + c

    # Combine the flux
    flux = [0,Pb*n[0],Pb*n[1],0]
    return flux, smag, Pb


