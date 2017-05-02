function u=exact_sol_linear_ivp(t,A,t0,eta) 
%Find the solution u(t) of the linear 
%IVP u?=Au, u(t0) = eta % A must be diagonalizable
[V,D]= eig(A); 
lam = diag(D); %eigenvalues 
c = V\eta; 
lt = length(t); 
trow = zeros(1,lt); 
trow(1:end) = t(1:end); 
m = max(size(A)); 
u = zeros(m,lt); 
for j = 1:m 
    erow = exp(lam(j)*(trow-t0)); 
    u = u+ c(j)*V(:,j)*erow; 
end
end