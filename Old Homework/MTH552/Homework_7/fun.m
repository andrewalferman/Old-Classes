function f=fun(x,y,A)
%derivatives for IVP 
global count;% counts how often fun is called
f = zeros(size(y));
f = A * y;
count = count+1;