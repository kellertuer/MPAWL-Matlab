function isMatrixValid(M)
% isMatrixValid(M) check whether M is integral, regular and quadratic.
%
% INPUT
%   M: The matrix to check
%
% OUTPUT
%   none, because if anithing is invalid, the method will provide an error
%
% ---
% MPAWL 1.0, R. Bergmann ~ 2013-09-11 ~last edited 2014-08-20

%2D array?
assert(all(size(size(M)) == [1,2]),'MPAWL:isMatrixValid', ...
    'The input matrix M is not a matrix');
%both dimensions identicyl?
assert(size(M,1) == size(M,2),'MPAWL:isMatrixValid', ...
    'The matrix is not quadratic');
%Integer?
assert(all(all(round(M) == M)),'MPAWL:isMatrixValid',...
    'The matrix is not an integer matrix');
%regular matrix?
assert(det(M) ~= 0, 'MPAWL:isMatrixValid', 'The matrix is not regular');
end

