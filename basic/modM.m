function h = modM(k,M,varargin)
%modM(k,M) Compute k mod M, there the target specifies the set of 
%      congruence class representants
%
%   INPUT
%       k: a d-dimensional vector
%       M: a dxd regular integral matrix
%       target (optional): specify the congrunce class representants
%           'unit' (standard): [0,1)^d or 'symmetric': [-0.5,0.5)^d
%
%   OUTPUT
%       h: the vector congruence class representant congruent to k
%
% OPTIONAL PARAMETERS
%   'Validate' : (true) whether or not to validate the input
%   'Target'   : ('unit') whether to produce a basis for the 'unit' cube or
%               the nearly 'symmetric' block. 
%   'Index'    : (false) whether or not the result is used as an index an
%                hence required to be integer. For generatingSetElements,
%                there might be rounding errors.
%--------------------------------------------------------------------------
% MPAWL 1.0, written on 2013-09-11 by Ronny Bergmann

p = inputParser;
addParameter(p, 'Validate',true,@(x) islogical(x));
addParameter(p, 'Target','unit');
addParameter(p, 'Index',false);
parse(p, varargin{:});
ppV = p.Results.Validate;
target = p.Results.Target;
index = p.Results.Index;
if (ppV)
    isMatrixValid(M);
end
if strcmp(target,'unit')
    h = M*mod(M\k,1);
elseif strcmp(target,'symmetric')
    h = M*(mod(M\k+0.5,1)-0.5);
else
    error('Congruence class type (Target=''%s'') is unknown.',target);
end
if index
    h = round(h); %just ensure integers for adressing in array
end
end