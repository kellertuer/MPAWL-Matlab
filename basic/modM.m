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
%--------------------------------------------------------------------------
% MPAWL 1.0, written on 2013-09-11 by Ronny Bergmann

p = inputParser;
addParamValue(p, 'Validate',true,@(x) islogical(x));
addParamValue(p, 'Target','unit');
parse(p, varargin{:});
ppV = p.Results.Validate;
target = p.Results.Target;
if (ppV)
    isMatrixValid(M);
end
if strcmp(target,'unit')
    h = M*mod(M\k,1);
elseif strcmp(target,'symmetric')
    h = M*(mod(M\k+0.5,1)-0.5);
else
    error('Congruence class type (target=''%s'') is unknown.',target);
end
end