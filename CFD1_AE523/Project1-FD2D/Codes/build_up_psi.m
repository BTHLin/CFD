function [psi] = build_up_psi(N,ds,dtheta,R0,R1,U_inf)

psi = zeros(N+1,N+1);

for i=1:N+1
    for j=1:N+1
        
        if i>=2 && i<=N && j>=2 && j<=N
            s=(j-1)*ds;
            r=R0*(R1/R0)^s;
            y=r*sin((i-1)*dtheta);
            psi(i,j) = U_inf*y;
        end        
    end
    psi(i,N+1) = U_inf*R1*sin((i-1)*dtheta);
end
    
end