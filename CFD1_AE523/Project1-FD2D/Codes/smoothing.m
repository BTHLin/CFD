function [Psi] = smoothing(Psi,C1,C2,N,f)


if N==4
    for iter=1:10
        for i=2:N
            for j=2:N
                % Red points
                if mod(i+j,2)==0
                    Psi(i,j) = (2*(C1+C2))^-1*...
                        (C2*(Psi(i+1,j)+Psi(i-1,j)) + ...
                        C1*(Psi(i,j+1)+Psi(i,j-1)) - f(i,j));
                end
            end
        end
    
        for i=2:N
            for j=2:N
                % Black points
                if mod(i+j,2) ~= 0
                    Psi(i,j) = (2*(C1+C2))^-1*...
                        (C2*(Psi(i+1,j)+Psi(i-1,j)) + ...
                        C1*(Psi(i,j+1)+Psi(i,j-1)) - f(i,j));
                end
            end
        end
    end 
    
else
    
    for iter=1:2     
        for i=2:N
            for j=2:N        
                % Red points
                if mod(i+j,2)==0
                    Psi(i,j) = (2*(C1+C2))^-1*...
                        (C2*(Psi(i+1,j)+Psi(i-1,j)) +...
                        C1*(Psi(i,j+1)+Psi(i,j-1)) - f(i,j));
                end
            end
        end
        
        for i=2:N
            for j=2:N
                % Black points
                if mod(i+j,2) ~= 0
                    Psi(i,j) = (2*(C1+C2))^-1*...
                        (C2*(Psi(i+1,j)+Psi(i-1,j)) + ...
                        C1*(Psi(i,j+1)+Psi(i,j-1)) - f(i,j));
                end
            end
        end
    end
end

end