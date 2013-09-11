function v = isMatrixValid(M)
% isMatrixValid(M) check whether M is integral, regular and quadratic.
%
%    INPUT
%        M: The matrix to check
%    OUTPUT
%        v: true is the matrix is valid, else 0,
%           but the function will then also abort with an error.
%
%--------------------------------------------------------------------------
% MPAWL 1.0, written on 2013-09-11 by Ronny Bergmann
v = 0;
if size(size(M)) ~= [1,2] %2D array?
    error('The input M is not a matrix');  
end
if size(M,1) ~= size(M,2) %both dimensions identicyl?
    error('The matrix is not quadratic');
end
if any(any(round(M) ~= M)) %Integer?
    error('The matrix is not an integer matrix');
end
if det(M) == 0 %regular matrix?
    error('The matrix is not regular');
end
v = 1;
end

