function [ is_int ] = is_integer_matrix( A )
% written by bz
% tests if each element in the matrix A is an integer

C = abs( round(A) - A ) > eps(A);

if sum(C(:)) > 0
    is_int = false;
else
    is_int = true;
    
end

