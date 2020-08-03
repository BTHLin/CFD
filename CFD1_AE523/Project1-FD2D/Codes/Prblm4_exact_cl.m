clear;clc;clf;close all;t=cputime;


%========================Intitail Setup===================================%
Nmax=256;Cl=1;
for n=2:log2(Nmax)
    N = 2^n; Na=(N+1)^2;
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


% GRID GENERATE
%============================Global Domain================================%

x=0;x0=0;x1=0;
y=0;y0=0;y1=0;
ang=0:0.01:pi;
d_r=0;figure(1)

for i=1:N+1
    dx = (R0+d_r)*cos(ang);
    dy = (R0+d_r)*sin(ang);
    d_r = d_r+(R1-R0)/N;
    if i==1 || i==N+1
        % plot out the hangar and farfield boundary    
        plot(x+dx,y+dy,'k'); hold on
    else
        % plot out the semi-circles per r+dr
        plot(x+dx,y+dy,'g'); hold on
    end
    
    if i<=N
        % plot out the radial lines per theta+dtheta
        x0 = R0*cos((i/N)*pi);
        y0 = R0*sin((i/N)*pi);
        x1 = R1*cos((i/N)*pi);
        y1 = R1*sin((i/N)*pi);
        plot([x0 x1],[y0 y1],'g'); hold on
    else
        % plot out the ground boundary
        plot([R0 R1],[y y],'k'); hold on
        plot([-R0 -R1],[y y],'k');hold on
    end
end

%
%SOLVING THE PDE EQUATION      
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
        m=(j-1)*(N+1)+1; q(m,1)=0;
        %farfield
        m=(j-1)*(N+1)+N+1; q(m,1)=R1*sin((j-1)*dtheta);
        %grounds
        q(j,1)=0; m=N*(N+1)+j; q(m,1)=0;
        
end

%==========================Solve the Psi Matrix===========================%

psi=A\q;

for j=1:N+1
    for i=1:N+1
        m=(j-1)*(N+1)+i;
        Psi(N+2-j,i)=psi(m,1);
        x(j,i)=(R0+(i-1)*dr)*cos(theta0+(j-1)*dtheta);
        y(j,i)=(R0+(i-1)*dr)*sin(theta0+(j-1)*dtheta);        
    end
end

% POST PROCESSES
% =================== Calculate the Pressure Coefficient ================ %
for i=1:N+1
    V_theta = -(-Psi(i,3)+4*Psi(i,2)-3*Psi(i,1))/(2*dr);  
    p=p_inf + Q - 0.5*rho*V_theta^2; %pressure on the hanger  
    Cp(i)=(p-p_inf)/Q; CP{n-1}=Cp;   
end

% =================== Calculate the Coefficient of lift ================= %
Cl_new=((2*R0)^-1)*trap_int(-Cp.*sin(theta0:dtheta:theta1),theta0,theta1,N);
error_Cl(n-1)=abs(Cl_new-Cl)/Cl_new; Cl=Cl_new;
end

% Plot the streamlines
figure(1)
contour(x,y,Psi); axis equal tight; xlabel('r')
title(['Streamlines for N=' num2str(N)])

% Plot the contour
figure(2)
contourf(x,y,Psi); axis equal tight; xlabel('r')
title(['Contour of Stream Function for N=' num2str(N)]); colorbar

% Plot the Cl converges vs. N number
figure(3)
loglog(2.^(1:n-1),error_Cl);grid on
xlabel('Number of N');ylabel('Error of Cl')
title('Error of Cl vs. Number of N')

% Plot the Pressure Coefficient
figure(4)
x_p=linspace(-R0,R0,N+1);
plot(x_p,Cp,'--');xlabel('r');ylabel('Cp');
title('Coefficient of pressure')

calculating_time=cputime-t














