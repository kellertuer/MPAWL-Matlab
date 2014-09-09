function s = generatingSet(M,varargin)
%generatingSet(M) 
%  returns the set of points belonging to a generating set of M.
% 
%   INPUT
%       M : regular integral matrix
%
%   OUTPUT
%      s  : set of (integer) vectors of the generating set
%
%    OPTIONAL ARGUMENTS
%        'Target' : ('symmetric') specifies whether the unit cube ('unit')
%                   or the nearly symmetric cube ('symmetric') is used for
%                   the result.
%        'Validate' (true) : whether to validate input or not
%
% ---
% MPAWL 1.0, R. Bergmann ~ 2013-11-15 ~ last edit: 2014-08-19

s = M*pattern(patternNormalForm(M),varargin{:});

end