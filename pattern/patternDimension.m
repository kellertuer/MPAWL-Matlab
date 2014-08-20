function dM = patternDimension(M)
% dM = patternDimension(M)
% Returns the number of elementary divisors of mM, that are greater than 1.
% This corresponds to the number of basis vectors of the pattern and hence
% the dimension of the corresponding lattice.
%
% INPUT
%   M : a quadratic integer matrix 
%
% OUTPUT
%  dM : dimension of the lattice corresponding to M
%
% ---
% MPAWL 1.0, R. Bergmann ~ 2014-08-19

dM = sum(diag(snf(M)) > 1);

end

