
% mult_matrix.m

function [C] = add_matrix(A,B)

% check if they are of the same size...
[m1,n1] = size(A);
[m2,n2] = size(B);

if ( n1 ~= m2 ),
    disp('ERROR: Rank of matrices do not match to multiply.');
    return;
else
    C = zeros(m1,n2);
end;


% now loop over elements and add them...
for ( i = 1:m1 ),
    for ( j = 1:n2 ),
        C(i,j) = 0;
        for ( k = 1:n1 ),
            C(i,j) = C(i,j) + A(i,k)*B(k,j);
        end;
    end;
end;


