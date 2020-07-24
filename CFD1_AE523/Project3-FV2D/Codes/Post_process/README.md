Define 3 functions:  
1. Postprocess, takes (P, Parameters, U0) returns (Cl_total, Cp, Cl, Cd), where  
    - P is a 2D vector contains [index no., pressure value, normal vector, length] on each edge  
    - Parameter defines air density, angle of attack, free stream Mach number, specific heat ratio, and airfoil chord length  
    - U0 is initial states deines air density, velocity, and kinetic energy in each element  
  
2. plotstate: takes (u,Cl,Cd,Cp,Cl_con,Res_con,Iter,midpoint,strlineIE,strlineBE,parameter,mesh,start) to visulize the solution and save into a file.  
  
3. fastrun: takes (order, mesh, and LV_plot) to plot pressure and mesh contour from the previous solution.  
  
4. airfoilplot: Plots airfoil on the final result.  
