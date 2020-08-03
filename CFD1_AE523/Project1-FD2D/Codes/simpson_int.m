function [fnc_ans]=simpson_int(int_fnc,a,b,n)

% int_fnc = integral function
% a = start value
% b = end value
% n = number of segments

dx=(b-a)/n;

sum_f=0;
for i=2:n-1
    if mod(i,2)==0
        sum_f=sum_f+4*int_fnc(i); % even part
    else
        sum_f=sum_f+2*int_fnc(i); % odd part
    end       
end

fnc_ans=(dx/3)*(int_fnc(1)+int_fnc(n)+sum_f);
end