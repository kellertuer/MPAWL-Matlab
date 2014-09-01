function hatb = bracketSums(data,origin,M,varargin)
% bracketSums(data,originindex,M)
% Compute bracket sums for a given array of data, where the index 0 is
% placed at origin values. These sums are computed for every index h from
% the generating Set (index by the generatingSetBasis).
%
% INPUT
%   data   : an array of values, which is d-dimensional
%   origin : entry in data that is adressing the index 0 in data
%   M      : matrix 
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
    assert(length(vector)==d,'vector is of invalid size (%d), %d required',length(origin),d))
    assert(all(origin>0),'origin has to be a positive vector');
    assert(all(origin<size(data)),'origin has to be in range of data');
end
m = abs(det(M))
d = size(M,1);
dM = patternDimension(M);
epsilon = diag(snf(M));
epsilon = epsilon(d-dM+1:d);
hM = generatingSetBasis(transpose(M))
%Compute maximal values
tmax = zeros(1,d);
tmaxsum = NestedFor(zeros(d,1),ones(d,1));
while tmaxsum.hasNext())
    omega = abs(transpose(M)*(tmaxsum.next - 0.5*ones(d,1)));
    tmax = max(tmax,omega);
end
torigin = tmax+1;
sums = zeros(2*tmax+1);
%run over all indices of data
summation = nestedFor(ones(1,size(size(data))),size(data));
while (summation.hasNext())
    epsilon = summation.next();
    epsMod = modM(epsilon-origin,transpose(M),'target','symmetric')+torigin;
    if strcmp(pp.Compute,'absolute squares')
        sums(subd2ind(epsMod,2*tmax+1)) = sums(subd2ind(epsMod)) + abs(data(epsilon))^2;    
    else % else default: brackets
        sums(subd2ind(epsMod,2*tmax+1)) = sums(subd2ind(epsMod)) + data(epsilon);
    end
    % put result in right circle order
    hatb = zeros(epsilon);
    hM = generatingSetBasis(transpose(M));
    summation = NestedFor(ones(1,dM),epsilon);
    while (summation.hasNext())
        ind = summation.next();
        hatb(sub2ind(ind,epsilon)) = sums(sub2ind(symMod(epsilon*hM,transpose(M),'Target','symmetric')+torigin,2*tmax+1));
    end
end

%collect
summation = nestedFor(zeros(size(epsilon)),epsilon);

end