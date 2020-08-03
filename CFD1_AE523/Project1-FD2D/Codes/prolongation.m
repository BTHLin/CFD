function [Psi_h] = prolongation(Psi_H,N)

A = zeros(2*N+1,2*N+1);
Psi_h = zeros(2*N+1,2*N+1);
for i=1:N+1
    for j=1:N+1
        
        A(1+(i-1)*2,1+(j-1)*2) = Psi_H(i,j);     

    end
end

for i=2:2*N
    for j=2:2*N
        if mod(i,2)~=0 && mod(j,2)~=0
            Psi_h(i,j) = A(i,j);
        elseif mod(i,2)==0 && mod(j,2)==0
            Psi_h(i,j) = 0.25*...
                (A(i-1,j-1)+A(i+1,j-1)+A(i-1,j+1)+A(i+1,j+1));
        elseif mod(j,2)~=0
            Psi_h(i,j) = 0.5*(A(i+1,j)+A(i-1,j));
        else
            Psi_h(i,j) = 0.5*(A(i,j+1)+A(i,j-1));
        end
    end
end


end