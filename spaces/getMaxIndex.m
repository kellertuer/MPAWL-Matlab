function maxInd = getMaxIndex(M,varargin)
% getMaxIndex(M)
%    returns the maximal index (vector) an entry of the (symmetric)
%    generating Set uses when looking at G(M).
%
%    INPUT
%       M     : a regular integer matrix with nonzero determinant
%
%    OUTPUT
%      maxInd : the maximum integer (elementwise) vector used in G(M)
%               using the symmetric definition
%
%    OPTIONAL ARGUMENTS
%        'Target' : ('symmetric') specifies whether the unit cube ('unit')
%                   or the nearly symmetric cube ('symmetric') is used for
%                   the result.
%        'Validate' (true) : whether to validate input or not

p = inputParser;
addParamValue(p, 'Validate',true,@(x) islogical(x));
addParamValue(p, 'Target','Bracket');
parse(p, varargin{:});
pp = p.Results;
if (pp.Validate)
    isMatrixValid(M);
end
d = size(M,1);
shift = 0.5*ones(1,d);
if (strcmp(pp.Target,'unit')) %lazy
    shift = zeros(1,d);
end
maxInd = zeros(1,d);
mysum = nestedFor(zeros(1,d),ones(1,d));
while mysum.hasNext()
    omega = abs(ceil(M*(mysum.next()' - shift')))'+1;
    maxInd = max(maxInd,omega);
end
end

