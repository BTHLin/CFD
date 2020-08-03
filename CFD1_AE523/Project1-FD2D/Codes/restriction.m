function [r_H] = restriction(r_h,N)


% step 1: calculate the update
% step 2: restrict the metrix

A = r_h;
% step 1
for i=2:N
    for j=2:N
        A(i,j)=...
            0.0625*(r_h(i-1,j-1)+r_h(i-1,j+1)+r_h(i+1,j-1)+r_h(i+1,j+1))+...
            0.125*(r_h(i-1,j)+r_h(i+1,j)+r_h(i,j-1)+r_h(i,j+1)) + ...
            0.25*r_h(i,j);
    end
end

% step 2
N = N/2;
r_H = zeros(N+1,N+1);
for i=1:N+1
    for j=1:N+1
%         if i>=2 && i<=N && j>=2 && j<=N
%             r_H(i,j) = A(1+(i-1)*2,1+(j-1)*2);
%         else
            r_H(i,j) = A(1+(i-1)*2,1+(j-1)*2);
%         end
    end
end



end


