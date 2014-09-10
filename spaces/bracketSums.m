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
%   'Compute'  : ('bracket') compute different sums, other values:
%                 'absolute suqares'
%
% ---
% MPAWL 1.0, R. Bergmann ~ 2014-09-30

p = inputParser;
addParamValue(p, 'Validate',true,@(x) islogical(x));
addParamValue(p, 'Compute','Bracket');
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
m = abs(det(M))
d = size(M,1);
dM = patternDimension(M);
epsilon = diag(snf(M));
epsilon = epsilon(d-dM+1:d);
%Compute maximal values
tmax = getMaxIndex(transpose(M));
torigin = tmax+1;
sums = zeros(2*tmax+1);
%run over all indices of data
summation = nestedFor(ones(1,size(size(data),2)),size(data));
while summation.hasNext()
    index = summation.next();
    epsMod = num2cell(modM((index-origin)',transpose(M),'target','symmetric')'+torigin);
    index = num2cell(index);
    if strcmp(pp.Compute,'absolute squares')
        sums(epsMod{:}) = sums(epsMod{:}) + abs(data(index{:}))^2;
    else % else default: brackets
        sums(epsMod{:}) = sums(epsMod{:}) + data(index{:});
    end
end
% put result in right circle order
if (sum(size(epsilon))==2) % one cycle
    hatb = zeros(1,epsilon);
else
    hatb = zeros(epsilon);
end
hM = generatingSetBasis(transpose(M));
summation = nestedFor(zeros(1,dM),epsilon-ones(1,dM));
while (summation.hasNext())
    ind = num2cell(summation.next()+1);
    sumInd = num2cell(modM(epsilon*hM,transpose(M),'Target','symmetric')'+torigin);
    hatb(ind{:}) = sums(sumInd{:});
end
end