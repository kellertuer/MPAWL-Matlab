function s = pattern(M,varargin)
% pattern(M)
%   generate the pattern of a matrix M.
%
%    INPUT
%        M: Input matrix, which has to be in pattern normal form
%
%    OUTPUT
%        s: set of vectors that belong to the pattern of M
%
%    OPTIONAL ARGUMENTS
%        'Target' : ('symmetric') specifies whether the unit cube ('unit')
%                   or the nearly symmetric cube ('symmetric') is used for
%                   the resulting set of pattern points.
%        'Validate' : (true) whether to validate input or not (for
%        performance reasons)
%
% ---
% MPAWL 1.0 ~ R. Bergmann ~ 2013-11-15 ~ last edit: 2014-08-19
p = inputParser;
addParameter(p, 'Target','symmetric');
addParameter(p, 'Validate',true,@(x) islogical(x));
parse(p, varargin{:});
pp = p.Results;
if (pp.Validate)
    isMatrixValid(M);
    if ~all(all(M == patternNormalForm(M)))
        error('MPAWL:pattern',...
            'The matrix M is not provided in pattern normal form');
    end
end

d = size(M,1);
step = 1/abs(M(d,d));
min = ceil( - (1/step)/2);
max = ceil( (1/step)/2 )-1;
s = (min:max)*step;
if (d>1)
    s = recPattern(M,d-1,s);
end
if (strcmp(pp.Target,'unit')) %lazy, sry
    s = mod(s,1);
end
end
function s=recPattern(M,dim,s_pre)
% recursively generate patterns in each dimension starting from the last
%
d = size(M,1);
l = d-dim+1;
step = 1/abs(M(dim,dim));
s = zeros(l,size(s_pre,2)/step);
% shifts by elements of the upper right of the matrix
for i=1:size(s_pre,2)
    tsum = M(dim,dim+1:d)*s_pre(:,i);
    min = ceil( - (1/(2*step)) + tsum);
    max = ceil( (1/(2*step)+tsum) ) - 1;
    s(1,(1/step*(i-1)+1):(1/step*i)) = ((min:max) - tsum)*step;
    s(2:l,(abs(M(dim,dim))*(i-1)+1):(abs(M(dim,dim))*(i))) = repmat(s_pre(:,i),1,M(dim,dim));
end
if (dim>1)
    s = recPattern(M,dim-1,s,target);
end
end