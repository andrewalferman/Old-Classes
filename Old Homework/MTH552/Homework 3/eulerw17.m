function [t,U] = eulerw17(odefun, TSPAN ,U0,NSTEP, RKmethod);
T0 = TSPAN(1); TFINAL = TSPAN(2);
dt = (TFINAL-T0)/NSTEP;
m = length(U0);
U = zeros(m,NSTEP+1);
U(:,1) = U0;
t = T0:dt:TFINAL;
if RKmethod == 'RK4'
    A = [0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0];
    b = [1/6;1/3;1/3;1/6];
    c = [0;0.5;0.5;1];
end
for k = 1:NSTEP
    U(:,k+1) = RKexplicitstep(odefun,t(k),U(:,k),dt, A, b, c);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function U = eulerstep(odefun,t,U0,dt)
U = U0 + dt*feval(odefun,t,U0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yprime = myode1(t,y)
% Implements ODE y''= -y as a first order system
yprime = zeros(size(y));
yprime(1) = y(2);
yprime(2) = -y(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = RKexplicitstep(odefun,t,U0,dt, A, b, c)
% One step of a general explicit Runge Kutta method
% A is assumed to be strictly lower triangular
% b and c must be column vectors
%
r = length(b);
m = length(U0);
K = zeros(m,r);%matrix whose j-th column is K_j
K(:,1) = feval(odefun,t,U0);
for j = 2:r
    Y = U0 + dt*K(:,1:j-1)*(A(j,1:j-1).');
     % Y = U0 + dt*sum_{l=1}^{j-1} a_{jl}K_l
    K(:,j) = feval(odefun,t + c(j)*dt, Y);
end

U = U0 + dt*K*b;
end