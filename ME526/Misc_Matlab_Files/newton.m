
% NEWTON.M 
% Use Newton?s method to find the roots of a function f. 
% As new iterates are computed, they are added to the 
% bottom of the array named iters. 

f = 'x.^3 - 4*x.^2 + 3*x + 1'; 
fprime = '3*x.^2 - 8*x + 3'; 
x = [0 1.5 3]; 

iters = x; 

for k = 1:20, 
    
    xnew = x - eval(f) ./ eval(fprime); 
    iters = [iters; xnew]; 
    
    if max(abs( (xnew - x)./xnew )) < 10*eps,
        break, 
    end; 
    
    x = xnew; 

end

% Switch to long format, display the final answer, and then
% turn off the long format. 

format long, iters, format 
