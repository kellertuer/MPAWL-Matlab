function hatb = patternFFT(M,b,varargin)
% patternFFt
% Compute the Fourier transform with respect to the pattern P(M), where
% b is the data given as a matrix that is running through the basis vectors
% w.r.t patternBasis(M) or the same corresponding reshaped vector b.
%
% INPUT
%   M : regular integer matrix
%   b : data vector or matrix (or b(:))
%
% OUTPUT
%   hatb : resulting Fourier transform
%
% ---
% MPAWL 1.0, R. Bergmann ~ 2014-08-28
    p = inputParser;
    addParamValue(p, 'Validate',true,@(x) islogical(x));
    parse(p, varargin{:});
    opt = p.Results;
    if (opt.Validate)
        isMatrixValid(M);
    end
    dM = patternDimension(M);
    d = size(M,1);
    S = diag(snf(M));
    e = S(d-dM+1:d);
    detM = prod(e);
    assert(numel(b)==detM,...
        'wrong number of elements (%d) in the data vector b. %d elements required',...
        numel(b),detM);
   if isvector(b)
       intb = reshape(b,e');
   else
       intb = b;
   end
   debug('text',3,'Text',['Starting the ',num2str(dM),'-dim. Fourier Transform for b (size: ',num2str(size(b)),').']);
   debug('time',3,'StartTimer','pfft');
   inthatb = 1/sqrt(det(M))*fftn(intb); %dM-dimensional Fourier transform
   debug('time',3,'StopTimer','pfft');
   if isvector(b)
       hatb = inthatb(:);
   else
       hatb = inthatb;
   end
end