clear;
close all;

scatterpoints = 1000000;

X = zeros(1,scatterpoints);
Y = zeros(1,scatterpoints);

xmin = 999;
mindeviance = 9999;

% Be sure not to use the variable i for indexing here, otherwise you will
% surely spend way too much time figuring out why values of z are so huge.
for k = 1:scatterpoints
   angle = (k / scatterpoints) * 2 * pi;
   z = ((-3/4)-(exp(i*angle)/2)+(exp(2i*angle)/4)+exp(3i*angle))/...
       ((5/8)+(19*exp(2i*angle)/8));
   X(k) = real(z);
   Y(k) = imag(z);
   
   check = abs(real(z) - imag(z));

   if real(z) < xmin
       xmin = real(z);
       yleft = imag(z);
   end
   
   if ((check < mindeviance) && (X(k) < -0.05))
       xval = X(k);
       yval = Y(k);
       mindeviance = check;
   end
end

hold on;
scatter(X,Y,0.3);
grid on;
xlabel('Real Part', 'FontSize', 14);
ylabel('Imaginary Part', 'FontSize', 14);
set(gca,'FontSize',12)

A = [200 398 198; -500 -696 -296; 500 694 294];
U0 = [2.6726e-1; -5.3452e-1; 8.0178e-1];

eigenvals = eig(A);