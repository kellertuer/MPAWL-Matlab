function [ckphi, BSums] = dirichletKernel(M,varargin)
% DirichletKernel(M)
% Generate the Dirichlet Kernel with respect to the translates of 
% pattern(M).
%
% INPUT
%    M : matrix for the translates of the de la Vallée Poussin means
%
% OUTPUT
%   ckphi : Fourier coefficients of the de la Vallée Poussin means, where
%   all dimensions have odd number of entries, and c_0 is at the center
%   BSums : Corresponding Bracket sums
%
% OPTIONAL PARAMETERS
%   'Validate' : (true) whether or not to validate the input
%   'File'     : (string) or ({string,string}) whether or not to load (or
%                if not possible ty try to save) the coefficients. If a
%                second file is given, the same holds for the bracket sums.
%  'Orthonormalize': (true) whether or not to orthonormalize the translates
%
% ---
% MPAWL, R. Bergmann ~ 2014-09-19
if (nargout==1)
    ckphi = delaValleePoussinMean(0,M,varargin);
elseif (nargout==2)
    [ckphi,BSums] = delaValleePoussinMean(0,M,varargin);
else
    error('Too many return values');
end