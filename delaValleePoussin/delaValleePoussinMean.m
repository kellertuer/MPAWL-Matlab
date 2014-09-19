function [ckphi, BSums] = delaValleePoussinMean(g,M,varargin)
% delaValleePoussin(g,M)
% Generate the de la Vallée Poussin mean based on the function g with
% respect to the translates of pattern(M), where g has to be a partition of
% unity in the d-dimensional space, nonnegative and positive on the unit
% cube. For simplicity g may also be a number or vector, which reduces to
% using the pyramidFunction
%
% INPUT
%    g : a function or vector characterizing the de la Vallée Poussin mean
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
%  'Support'       : cube indicating the support, if g is a function. for
%                    the vector case, the support is determined
%                    automatically. If it is not given for g, the area
%                    M*unit cube is taken.
%
% ---
% MPAWL, R. Bergmann ~ 2014-09-18
p = inputParser;
addParamValue(p, 'Validate',true,@(x) islogical(x));
addParamValue(p, 'File','');
addParamValue(p, 'Orthonormalize',true,@(x) islogical(x));
addParamValue(p, 'Support',[]);
parse(p, varargin{:});
ppV = p.Results.Validate;
if (ppV)
    isMatrixValid(M);
end

end

