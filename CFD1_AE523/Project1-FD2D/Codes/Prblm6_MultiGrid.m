clear; clc; close all; t = cputime;

% step 1: set N_max, which must be a number power of 2, but large than 4
% step 2: run the program

% PARAMETERS
%============================Intitail Setup===============================%
Nmax =128;
%
N=Nmax;R0=1;R1=20;ds=1/N; 
theta0=0;theta1=pi;
p_inf=0;U_inf=1;rho=1;
Q=0.5*rho*U_inf^2;Cl=1;
%
% GRID GENERATE
%============================Global Domain================================%

figure(1)
ang = 0:0.01:pi;

for i=1:N+1
    s=(i-1)*ds; r=R0*(R1/R0)^s;
    x=0;dx = r*cos(ang);
    y=0;dy = r*sin(ang);
        
    if i==1 || i==N+1
        %plot out the hangar and farfield boundary  
        plot(x+dx,y+dy,'k');hold on
    else
        %plot out the semi-circles per r+dr  
        plot(x+dx,y+dy,'g');hold on
    end
    
    if i<=N
        %plot out the radial lines per theta+dtheta
        x0 = R0*cos(((i-1)/N)*pi);
        y0 = R0*sin(((i-1)/N)*pi);
        x1 = R1*cos(((i-1)/N)*pi);
        y1 = R1*sin(((i-1)/N)*pi);
        plot([x0 x1],[y0 y1],'g');hold on
    else
        %plot out the ground boundary
        plot([R0 R1],[y y],'k');hold on
        plot([-R0 -R1],[y y],'k');hold on
        axis equal tight
        title(['N=' num2str(N)])
    end
end

%
% SOLVING THE GOVERNING EQUATION    
%=========================MultiGrid Method V Cycle========================%

for n = 3:log2(Nmax)
    N=2^n; V=0; 
    % ======================== Restriction ============================== %
    while 1
        V=V+1; L=0;
        while N~=4
            % Parameters
            if L~=0; N=N/2;end
            L = L+1; %N= 2^n; %L:Level; N:Grid point
            d_s = 1/N; dtheta = (theta1-theta0)/N;
            C1 = (log(20)*d_s)^-2; C2 = dtheta^-2;  
        
            % Initialize the Psi matrix
            if L==1 && V==1
                Psi_h{L} = build_up_psi(N,d_s,dtheta,R0,R1,U_inf);
                fh = zeros(N+1,N+1); f{1} = fh;
            elseif L == 1
                Psi_h{L} =  Psi{V-1}; 
            else
                Psi_h{L} = zeros(N+1,N+1);
            end
   
            % smooth the Psi
            Psi_s{L} = smoothing(Psi_h{L},C1,C2,N,f{L});

            % Calculate the residual 
            r_h{L} = residual(Psi_s{L},C1,C2,N,f{L});

            % Restrcit to the next level
            fH = restriction(r_h{L},N);
            f{L+1} = fH;
        end

    % ======================== Prolongation ============================= %

        while L>=2   
            N = N*2;
            d_s = 1/N; dtheta = (theta1-theta0)/N;
            C1 = (log(20)*d_s)^-2; C2 = dtheta^-2;
            Psi_s{L-1} = Psi_s{L-1} + prolongation(Psi_s{L},N/2);  
            Psi_s{L-1} = smoothing(Psi_s{L-1},C1,C2,N,f{L-1});
            L = L-1;
        end

    % =================== Residual for a completed V cycle ============== %

        Psi{V} = Psi_s{L}; % store the updated Psi_s to the nth V cycle
        tol = 1E-18*ones(N+1,N+1);
        normL1(V) = sum(sum(abs(r_h{L})))/(N+1)^2;
                
    % ======================= Define the Work Unit ====================== %

        WU = 0; L = log2(N);
        for nlevl = 1:L   
            WU = WU + 4*2^(-(nlevl-1));
        end
        WU(n-2)=WU;

    % ================= Decide when to cut off the V-cycle ===============%
        if V == 20; break; end

    end

% POST PROCESSES
% =================== calculate the Pressure Coefficient ================ %
for i=1:N+1

    % Build up x,y matrix for ploting
    ds=1/Nmax; dtheta=pi/Nmax; s=(i-1)*ds;
    for j=1:N+1
        r=R0*(R1/R0)^s;
        x(j,i)=r*cos(theta0+(j-1)*dtheta);
        y(j,i)=r*sin(theta0+(j-1)*dtheta);       
    end
    
    % Calculate Coefficient of Pressure
    C1 = (log(20)*(20)^(2*ds))^-1; C2 = (log(20)*(20)^ds)^-1;
    V_theta = -(-C1*Psi_s{1}(i,3)+C2*4*Psi_s{1}(i,2)-3*Psi_s{1}(i,1))/(2*ds);
    p=p_inf + Q - 0.5*rho*V_theta^2 ;%pressure on the hanger  
    Cp(i)=(p-p_inf)/Q; CP{n-2}=Cp;
end

% =================== calculate the Coefficient of lift ================= %
Cl_new=((2*R0)^-1)*simpson_int(-Cp.*sin(theta0:pi/N:theta1),theta0,theta1,N);
error_Cl(n-2)=abs(Cl_new-Cl)/Cl; Cl=Cl_new;
norm{n-2} = normL1;
end

% Plot the streamlines
figure(1)
contour(x,y,Psi_s{1}); axis equal tight; xlabel('r')
title(['Streamlines for N=' num2str(N)]); 

% Plot the contour
figure(2)
contourf(x,y,Psi_s{1}); axis equal tight; xlabel('r')
title(['Contour of Stream Function for N=' num2str(N)]); colorbar
 
% Plot the error_Cl converges vs. N number
figure(3)
semilogy(2.^(3:n),error_Cl,'-.');grid on;
xlabel('Number of N');ylabel('Error of Cl');
title('Error of Cl vs. Number of N')

% Plot the Pressure Coefficient
figure(4)
x_p=linspace(-R0,R0,N+1); %x_p = x axis for pressure diagram
plot(x_p,Cp,'*');xlabel('location on the hanger');ylabel('Cp');
title('Coefficient of pressure')

% Plot the Norm vs. # of V Cycle
figure(5)
for n=1:n-2
     N = 2^(n+2)+1;
     x_p=linspace(-R0,R0,N); %x_p = x axis for pressure diagram
     Legend{n}=strcat('N=', num2str(2^(n+2)));
     semilogy(WU(n).*(1:V),norm{n},'-');hold on
end
title(['N=' num2str(Nmax)]);xlabel('Work Units');ylabel('||r_h|| L1 norm');
legend(Legend,'Location','NorthEast');hold on

% Plot the Norm vs. # of V Cycle
figure(6)
n=n+2;
for n=1:n-2
     N = 2^(n+2)+1;
     x_p=linspace(-R0,R0,N); %x_p = x axis for pressure diagram
     Legend{n}=strcat('N=', num2str(2^(n+2)));
     semilogy(1:V,norm{n},'-');hold on
end
title(['N=' num2str(Nmax)]);xlabel('Number of V cycle');ylabel('||r_h|| L1 norm');
legend(Legend,'Location','NorthEast');hold on

Cl
runtime = cputime-t








