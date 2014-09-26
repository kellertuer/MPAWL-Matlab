function [hata,hatb] = dirichletKernelSubspaces(M,varargin)
% [hata, hatb] = dirichletKernelSubspaces(M,J)
% Generate the Dirichlet Kernel complement with respect to % an M-invariant
% space, where the coefficients are given with respect to a dirichlet
% Kernel of the translates of M.
%
% INPUT
%    M : matrix for the translates of the de la Vallée Poussin means
%    J ; integer matrix indicating the subpattern, i.e. such that inv(J).M
%        is integer and |det(J)|=2.
%
% OUTPUT
%   hata : coeffs of the Dirichlet Kernel w.r.t. N=inv(J).M
%   hatb : coeffs of the unique orhogonal complement of hata
%
% OPTIONAL PARAMETERS
%   'Validate'     : (true) whether or not to validate the input
%   'File'         : (string) or ({string,string}) whether or not to load
%                    (or if not possible try to save) the coefficients.
%  'Orthonormalize': (true) whether or not to orthonormalize the both
%                    functions w.r.t. their translatestranslates. If both
%                    spaces are succesfully loaded, this option is ignored
% ---
% MPAWL, R. Bergmann ~ 2014-09-26
[hata,hatb] = delaValleePoussinSubspaces(0,M,J,varargin);
end

