function [r] = residual(Psi_s,C1,C2,N,f)

r = zeros(N+1,N+1);
for i=2:N
    for j=2:N
        r(i,j) = f(i,j) - (-2*(C1+C2)*Psi_s(i,j) + ...
            (C2*(Psi_s(i+1,j)+Psi_s(i-1,j)) + ...
            C1*(Psi_s(i,j+1)+Psi_s(i,j-1))));
    end
end

% m=i*(N+1)+j;
% r(m,1)=r(i,j);

end