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
% MPAWL, R. Bergmann, 2014-10-05

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
if (pp.Validate)
    assert(all(size(hata)==epsilon'),['The coefficient array hata (',...
        num2str(size(hata)),') is not of the correct size with respect to the cycles of M (',...
        num2str(epsilon')]);
end
% reorder
debug('time',3,'StartTimer','Generating Fourier coefficients from space coefficients');
ckf = zeros(size(ckphi));
griddims = cell(1,d);
for i=1:d
    griddims{i} = 1:size(ckf,i);
end
gridmeshes = cell(1,d);
[gridmeshes{:}] = ndgrid(griddims{:});
inds = zeros(d,numel(ckf));
for i=1:d
    inds(i,:) = reshape(gridmeshes{i},1,[]);
end
    gSetInds = generatingSetBasisDecomp(inds-repmat(origin',[1,numel(ckf)]),transpose(M),'Validate',false)+1;
gSetIndsc = cell(1,d);
indsc = cell(1,d);
for i=1:d
    gSetIndsc{i} = gSetInds(i,:);
    indsc{i} = inds(i,:);
end
ckf(sub2ind(size(ckf),indsc{:})) = hata(sub2ind(size(hata),gSetIndsc{:})).*ckphi(sub2ind(size(ckf),indsc{:}));
debug('time',3,'StopTimer','Generating Fourier coefficients from space coefficients');
end