In this project, the first and second-order finite volume method will be applied to solve a  
compressible, subsonic flow over a three-element, high-lift device airfoil. In post-process,  
the resulting forces, lift and drag coefficients, Mach and pressure contours and streamlines will  
be presented.  

To run the code, download the code folder and runs main.py  
In main.py, there are several settings can be changed.  
1) In the very button, one can select  
    a) ORDER: first(1) or second(2) order of the finite volume method.  
    b) MESH:  Dfine the mesh resolution, mesh0(1149 elements), mesh1(2116 elements), mesh2(4124 elements), mesh3(8031 elements)  
    c) fast_run: Use the pre-store data to plot contours  
    d) LV_plot: Wether or not to plot the linear-variation Mach and pressure contour.  
    d) airfoil_plot: Wether or not to plot the airfoil in the final result.  
2) In "main" function, one can change  
    a) Maximum iteration: maxiter  
    b) CFL number: CFL  
    c) Error tolerance: tolerance  
    d) Air density: rho  
    e) Angle of attack: alfa(deg)  
    f) Free stream Mach number: Minf
    g) Air Specific heat ratio (CP/CV): gamma
    h) Airfoil chord length: chord
    
