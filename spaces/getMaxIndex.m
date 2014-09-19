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
%        'CubeSize' : (ones(1,d)) cube size. 

p = inputParser;
addParamValue(p, 'Validate',true,@(x) islogical(x));
addParamValue(p, 'Target','symmetric');
d = size(M,1);
addParamValue(p, 'CubeSize',ones(1,d));
parse(p, varargin{:});
pp = p.Results;
if (pp.Validate)
    isMatrixValid(M);
end
shift = 0.5*pp.CubeSize;
if (strcmp(pp.Target,'unit')) %lazy
    shift = pp.CubeSize;
end
maxInd = zeros(1,d);
mysum = nestedFor(zeros(1,d),ones(1,d));
while mysum.hasNext()
    omega = abs(ceil(M*(mysum.next()' - shift')))'+1;
    maxInd = max(maxInd,omega);
end    
end

