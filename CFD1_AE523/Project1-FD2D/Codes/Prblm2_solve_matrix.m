clear;clc;clf;close all;t=cputime;

%========================Intitail Setup===================================%
N=3;
Na=(N+1)^2;
%
A = sparse(1:Na,1:Na,1);
q = zeros(Na,1);
%
R0=1;R1=20*R0;
dr=(R1-R0)/N;
%
theta0=0;theta1=pi;
dtheta=(theta1-theta0)/N;
%
p_inf=0;U_inf=1;rho=1;
Q=0.5*rho*U_inf^2; % dynamic pressure
%

%=========================Build up the matrices===========================%

for j=1:N+1
    
    for i=2:N        
        if j>=2 && j<=N 
        m=(j-1)*(N+1)+i;
        r=R0+(i-1)*dr;
        % A matrix
        A(m,m)=-2*(dr^-2+(r*dtheta)^-2);              % Note:
        A(m,m+1)=dr^-2 + (2*r*dr)^-1; % East side       R0+i*dr = r 
        A(m,m-1)=dr^-2 - (2*r*dr)^-1; % West side
        A(m,m+(N+1))=(r*dtheta)^-2;  % North side
        A(m,m-(N+1))=(r*dtheta)^-2;  % South side
        % Source term matrix
        q(m)=0;
        end        
    end  
    
    %==========================Bondary Conditions=========================%
        %hanger
        m=(j-1)*(N+1)+1;
        q(m,1)=0;
        %farfield
        m=(j-1)*(N+1)+N+1;
        q(m,1)=R1*sin((j-1)*dtheta);
        %grounds
        q(j,1)=0;
        m=N*(N+1)+j;
        q(m,1)=0;
        
end
%==========================Solve the Psi Matrix===========================%

psi=A\q;

for i=1:N+1
    for j=1:N+1
        m=(j-1)*(N+1)+i;
        Psi(N+2-j,i)=psi(m,1);
    end
end

Psi


























