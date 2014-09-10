function ckf = coeffsSpace2Fourier(hata,ckphi,origin,M,varargin)
% coeffsSpace2Fourier(hata,ckphiorigin,M)
% compute the Fourier coefficients of f, ckf, based on hata, the Fourier
% transform of its coefficients a wrt the translates T(y)phi, y from the
% pattern(M). 
%
%    INPUT
%        hata    : Fourier transform of the coefficients of the sum of
%                  T(y)phi yielding f.
%         ckphi  : Fourier coefficitents of phi
%         origin : origin, i.e. the index corr. to c_0 iin both above
%                  parameters ckf and ckphi
%         M      : matrix indicating the pattern for the translates of phi
%
%    OUTPUT
%         ckf    : Fourier coefficients of f
%
% ---
% MPAWL, R. Bergmann, 2014-09-10

p = inputParser;
addParamValue(p, 'Validate',true,@(x) islogical(x));
parse(p, varargin{:});
pp = p.Results;
if (pp.Validate)
    isMatrixValid(M);
end
hM = generatingSetBasis(transpose(M),'Target','symmetric','Validate',false);
d = size(M,1);
dM = patternDimension(M);
epsilon = diag(snf(M));
epsilon = epsilon(d-dM+1:d);

tmax = getMaxIndex(transpose(M));
torigin = tmax+1;

coeffsOI = Inf(2*tmax+1);
summation = nestedFor(zeros(1,dM),epsilon-ones(1,dM));
% reorder
while (summation.hasNext())
    ind = num2cell(summation.next()+1);
    sumInd = num2cell(modM(epsilon*hM,transpose(M),'Target','symmetric','Validate',false)'+torigin);
    coeffsOI(sumInd{:}) = hata(ind{:});
end
ckf = zeros(size(ckphi));
summation = nestedFor(ones(1,dM),size(ckphi));
while (summation.hasNext())
    ind = summation.next();
    indc = num2cell(ind);
    sumIndc = num2cell(modM((ind-origin)',transpose(M),'Target','symmetric','Validate',false)'+torigin);
    ckf(indc{:}) = coeffsOI(sumIndc{:})*ckphi(indc{:});
end
end

