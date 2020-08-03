clear;clc;clf

%=========================== Initial Setup ===============================%
N=16;
%
R0=1;R1=20*R0;
x=0;x0=0;x1=0;
y=0;y0=0;y1=0;
ang = 0:0.01:pi;

%============================Global Domain================================%

dr=0;
for i=1:N+1
    dx = (R0+dr)*cos(ang);
    dy = (R0+dr)*sin(ang);
    dr = dr+(R1-R0)/N;
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
        plot([-R0 -R1],[y y],'k')
        axis equal tight
    end
end

title(['N=' num2str(N)])
xlabel('r'); 


       





