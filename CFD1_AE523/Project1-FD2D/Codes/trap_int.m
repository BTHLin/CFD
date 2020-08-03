function [fnc_ans]=trap_int(int_fnc,a,b,n)

% int_fnc = integral function
% a = start value
% b = end value
% n = number of segments

dx=(b-a)/n;

sum_f=0;

for i=2:n-1 
    sum_f=sum_f+int_fnc(i);
end

fnc_ans=(dx/2)*(int_fnc(1)+int_fnc(n)+2*sum_f);

end