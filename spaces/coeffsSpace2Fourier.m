function ckf = coeffsSpace2Fourier(M,hata,ckphi,origin,varargin)
% coeffsSpace2Fourier(M,hata,ckphi,origin)
% compute the Fourier coefficients of f, ckf, based on hata, the Fourier
% transform of its coefficients a wrt the translates T(y)phi, y from the
% pattern(M). 
%
%    INPUT
%         M      : matrix indicating the pattern for the translates of phi
%         hata   : Fourier transform of the coefficients of the sum of
%                  T(y)phi yielding f.
%         ckphi  : Fourier coefficitents of phi
%         origin : origin, i.e. the index corr. to c_0 iin both above
%                  parameters ckf and ckphi
%
%    OUTPUT
%         ckf    : Fourier coefficients of f
%
%    OPTIONAL ARGUMENTS
%        'Validate' (true) : whether to validate input or not
%
%     NOTE 
%      The corresponding Mathematica function is called
%      'getFourierfromSpace' and was renamed to fit Matlab conventions
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
% reorder
debug('time',3,'StartTimer','Generating Fourier coefficients from space coefficients');
ckf = zeros(size(ckphi));
if (d==1)
    ind = 1:length(ckf);
    gSetInds = generatingSetBasisDecomp(ind-torigin,transpose(M))+1;
    indc = sub2ind(size(ckf),1:length(ckphi));
    ckf(indc) = hata(gSetInds).*(ckphi(:)');
elseif (d==2)
    [ind1,ind2] = meshgrid(1:size(ckf,1),1:size(ckf,2));
    gSetInds = generatingSetBasisDecomp([ind1(:)';ind2(:)']-repmat(origin',[1,numel(ckf)]),transpose(M),'Validate',false)+1;
    [ind1c,ind2c] = meshgrid(1:size(ckphi,1),1:size(ckphi,2));
    ckf(sub2ind(size(ckf),ind1c(:)'-origin(1)+torigin(1),ind2c(:)'-origin(1)+torigin(1))) = hata(sub2ind(size(hata),gSetInds(1,:),gSetInds(2,:))).*(ckphi(:)');
else %for higher dimensions i only now the slow (mathematica style) way of nested fors
    coeffsOI = Inf(2*tmax+1);
    summation = nestedFor(zeros(1,dM),epsilon'-ones(1,dM));
    while (summation.hasNext())
        ind = summation.next();
        indc = num2cell(ind'+1);
        sumIndc = num2cell(modM(hM*ind',transpose(M),'Target','symmetric','Validate',false,'Index',true)'+torigin);
        coeffsOI(sumIndc{:}) = hata(indc{:});
    end
    summation = nestedFor(ones(size(size(ckphi))),size(ckphi));
    while (summation.hasNext())
        ind = summation.next();
        indc = num2cell(ind');
        sumIndc = num2cell(modM((ind-origin)',transpose(M),'Target','symmetric','Validate',false,'Index',true)'+torigin);
        ckf(indc{:}) = coeffsOI(sumIndc{:})*ckphi(indc{:});
    end
end
debug('time',3,'StopTimer','Generating Fourier coefficients from space coefficients');
end