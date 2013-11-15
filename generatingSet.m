function s = generatingSet(M,varargin)
%generatingSet(M) 
%  returns the set of points belonging to a generating set of M.
% 
%    INPUT
%       M:: regular integral matrix
%
%    OUTPUT
%      s: set of (integer) vectors of the generating set
%
%    OPTIONAL ARGUMENTS
%        'Target' [='symmetric']
%            specifies whether the unit cube ('unit') or the nearly 
%            symmetric cube ('symmetric') is used for the result.
%        'Validate' [False]
%            whether to validate input or not
%        'Debug' ['None']
%            Different styles (&-concatenated strings) of Debug, e.g.
%            'Time', 'Text', 'Verbose', 'Figures',...
%   
%
% ---
% MPAWL 1.0 written by Ronny Bergmann
% created 2013-11-15

s = M*pattern(patternNormalForm(M),varargin{:});

end

