function b = changeBasis(M, a, bracketSums,varargin)
% hatb = changeBasis(M,hata, bracketSums)
% perform a change of basis on the coeffs a or hata, which are a DFT of a,
% coeffs w.r.t. a basis of translates and the bracket sums are in order to
% perform the change of basis and return the coefficients b or hatb of the
% result.
%
% INPUT
%   M           : a matrix indicating the pattern
%   a           : coeffs or their Fourier transform representing a function
%                 in the translates of phi
%   bracketSums : corresponding bracket sums for the change of basis
%
% OUTPUT
%   b           : coeffs or their Fourier transform representing a function
%                 in the translates of the new basis
%
% OPTIONAL PARAMETERS
%   'Validate' : (true) whether or not to validate the input
%   'Input'    : ('time') or 'Fourier' Domain of the input coefficients a
%   'Output'   : ('time') or 'Fourier' Domain of the output coefficients b
%
% ---
% MPAWL, R. Bergmann ~ 2014-09-17

p = inputParser;
addParameter(p, 'Validate',true,@(x) islogical(x));
addParameter(p, 'Input','time');
addParameter(p, 'Output','time');
parse(p, varargin{:});
ppV = p.Results.Validate;
if (ppV)
    isMatrixValid(M);
end
d = size(M,1);
dM = patternDimension(M);
epsilon = diag(snf(M)); epsilon = epsilon(d-dM+1:d);
if (ppV)
    if length(epsilon)>1
        assert(all(size(a)==epsilon'),...
            ['The required size for the coefficients is ''',num2str(epsilon'),''' but a is of size ''',num2str(size(a)),'''.']);
        assert(all(size(bracketSums)==epsilon'),...
            ['The required size for the Bracket sums is ''',num2str(epsilon'),''' but bracketSums is of size ''',num2str(size(bracketSums)),'''.']);
    else %1D case
        assert(all(numel(a)==epsilon'),...
            ['The required size for the coefficients is ''',num2str(epsilon'),''' but a is of size ''',num2str(size(a)),'''.']);
        assert(all(size(bracketSums)==epsilon'),...
            ['The required size for the Bracket sums is ''',num2str(epsilon'),''' but bracketSums is of size ''',num2str(size(bracketSums)),'''.']);
    end
end
if strcmp(p.Results.Input,'time')
    debug('text',3,'Text','Performing Fourier Transform on the input from time to Fourier domain');
    hata = patternFFT(M,a,'Validate',false);
elseif strcmp(p.Results.Input,'Fourier')
    hata = a;
else
    error('Unknoen domain for the input coefficients a');
end
debug('time',3,'StartTimer','changeBasis');
hatb = 1/abs(det(M))*bracketSums.*hata;
debug('time',3,'StopTimer','changeBasis');
if strcmp(p.Results.Output,'time')
    debug('text',3,'Performing inverse Fourier Transform on the output from time to Fourier domain');
    b = patternIFFT(M,hatb,'Validate',false);
elseif strcmp(p.Results.Output,'Fourier')
    b = hatb;
else
    error('Unknown domain for the output coefficients b');
end
end
