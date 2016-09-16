function hatb = bracketSums(data,origin,M,varargin)
% bracketSums(data,originindex,M)
% Compute bracket sums for a given array of data, where the index 0 is
% placed at origin values. These sums are computed for every index h from
% the generating Set (index by the generatingSetBasis).
%
% INPUT
%   data   : an array of values, which is d-dimensional
%   origin : entry in data that is adressing the index 0 in data
%   M      : an integer valued matrix indicating the pattern
%
% OUTPUT
%   hatb   : array of bracket sums adressed with respect to the basis of
%            M^T (given by generatingSetBasis(transpose(M))
%
% OPTIONAL PARAMETERS
%   'Validate' : (true) Whether to check the validtidy of input M and origin
%   'Compute'  : ('Bracket') compute different sums, other values:
%                 'absolute Squares', 'absolute Values'
%
% ---
% MPAWL 1.0, R. Bergmann ~ 2014-08-30 ~ 2014-09-15

p = inputParser;
addParameter(p, 'Validate',true,@(x) islogical(x));
addParameter(p, 'Compute','Bracket');
parse(p, varargin{:});
pp = p.Results;
if (pp.Validate)
    isMatrixValid(M);
end
d = size(M,1);
if (pp.Validate)
    assert(isvector(origin),'origin is not a vector.');
    assert(length(origin)==d,'vector is of invalid size (%d), %d required',length(origin),d);
    assert(all(origin>0),'origin has to be a positive vector');
    assert(all(origin<size(data)),'origin has to be in range of data');
end
dM = patternDimension(M);
epsilon = diag(snf(M));
epsilon = epsilon(d-dM+1:d);
%Compute maximal values
%run over all indices of data
debug('time',3,'StartTimer','the Bracket Sum');
griddims = cell(1,d);
for i=1:d
    griddims{i} = 1:size(data,i);
end
gridmeshes = cell(1,d);
[gridmeshes{:}] = ndgrid(griddims{:});
inds = zeros(d,numel(data));
for i=1:d
    inds(i,:) = reshape(gridmeshes{i},1,[]);
end
gSetInds = round(generatingSetBasisDecomp(inds-repmat(origin',[1,numel(data)]),transpose(M),'Validate',false)+1);
if strcmp(pp.Compute,'absolute Squares')
    hatb = accumarray(gSetInds',data(:).^2)';
else % else default: brackets
    hatb = accumarray(gSetInds',data(:))';
end
if(dM>1)
    hatb = reshape(hatb,epsilon');
end
debug('time',3,'StopTimer','the Bracket Sum');
end
