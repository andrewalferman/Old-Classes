
% qroots.m

function [x1,x2] = qroots(a,b,c)


delta = b^2 - 4*a*c;

if ( delta > 0 ),
    x1 = ( - b + sqrt(delta) ) / ( 2 * a );
    x2 = ( - b - sqrt(delta) ) / ( 2 * a );
elseif ( delta < 0 ),
    x1 = ( - b ) / ( 2 * a ) + sqrt(delta) / ( 2 * a ) * i;
    x2 = ( - b ) / ( 2 * a ) - sqrt(delta) / ( 2 * a ) * i;
else,
    x1 = ( - b ) / ( 2 * a );
    x2 = ( - b ) / ( 2 * a );
end
