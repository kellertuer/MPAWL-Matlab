function hata = coeffsFourier2Space(M,ckf,ckphi,origin, varargin)
% coeffsFourier2Space(M,ckf,ckphi,origin)
%  compute the coefficients of the space of translates (w.r.t. pattern(M))
%  of phi for f based on the given Fourier coefficients ckf and ckphi,
%  if possible (else a NaN-matrix is returned). 
%
%     INPUT
%         M      : matrix indicating the pattern for the translates of phi
%         ckf    : Fourier coefficients of f
%         ckphi  : Fourier coefficitents of phi
%         origin : origin, i.e. the index corr. to c_0 iin both above
%                  parameters ckf and ckphi
%
%    OUTPUT
%        hata    : Fourier transform of the coefficients of the sum of
%                  T(y)phi yielding f, if existent, else an entry of the
%                  result indicates the failure by NaN at the corresponding index.
%
%    OPTIONAL ARGUMENTS
%        'Validate' (true) : whether to validate input or not
%
%     NOTE 
%      The corresponding Mathematica function is called
%      'getFourierfromSpace' and was renamed to fit Matlab conventions
% ---
% MPAWL 1.0, R. Bergmann, 2014-09-10
p = inputParser;
addParameter(p, 'Validate',true,@(x) islogical(x));
parse(p, varargin{:});
pp = p.Results;
if (pp.Validate)
    isMatrixValid(M);
end
d = size(M,1);
dM = patternDimension(M);
epsilon = diag(snf(M));
epsilon = epsilon(d-dM+1:d);
debug('time',3,'StartTimer','space coefficients from Fourier.');
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
gSetInds = round(generatingSetBasisDecomp(inds-repmat(origin',[1,numel(ckf)]),transpose(M),'Validate',false)+1);
data = ckf./ckphi;
hata = conj(accumarray(gSetInds',data(:),epsilon',@checkGroup,NaN)');
%if(dM>1)
%    hata = reshape(hata,epsilon');
%end
if any(isnan(hata))
    warning('the given data does not seem to be in the space of translates, the coefficients contain NaNs!');
end
debug('time',3,'StopTimer','space coefficients from Fourier coefficients.');
end

function v = checkGroup(x)
% checkGroup(x)
% check whether a group of Fourier coefficients (belonging to the same
% congruence class) is a valid group, i.e. there is only a unique nonzero
% value or all are zero, then this value is the result, else NaN is
% returned
if (numel(x)==0)
    v = 0;
elseif (numel(unique(x))==1)
    v=unique(x);
else
    if any(isinf(x)) %ckf nonzero ckphi zero
        v = NaN;
    elseif any(isnan(x)) %only okay is all are NaN
        if all(isnan(x))
            v = 0;
        else
            v = NaN;
        end
    else
        %the only things that may happen is nans and we can ignore these
        v = unique(x(~isnan(x)));
        if (numel(v)>1) %if v is not unique, we return nan.
          if max(abs(mean(v)-v))>eps
              v = NaN;
          else
            v = mean(v);
          end
        end
    end
end
end