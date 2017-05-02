
% mult_matrix.m

function [C] = add_matrix(A,B)

% check if they are of the same size...
[m1,n1] = size(A);
[m2,n2] = size(B);

if ( m1 ~= m2 | n1 ~= n2 ),
    disp('ERROR: Rank of matrices do not match for addition.');
    return;
else
%    C = zeros(m1,n1);
end;


% now loop over elements and add them...
for ( i = 1:m1 ),
    for ( j = 1:n2 ),
        C(i,j) = A(i,j) + B(i,j);
    end;
end;


