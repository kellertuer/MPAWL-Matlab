function V = patternBasis(M,varargin)
% patternBasis(M)
% creates a set of vectors (columns of V), whose integer multiples (up to
% the corresponding elementary divisor - 1) create the pattern. The vectors
% are ordered w. r. t. the ordering of the elementary divisors, i.e. w.r.t.
% their nondecreasing ordering from the diagonal of snf(M)
%
% INPUT
%   M : an integer matrix n times n
%   
% OUTPUT
%   V : a matrix containing the pattern basis vectors, hence its format is
%   n times patternDimension(M)
%
% OPTIONAL PARAMETERS
%   'Validate' : (true) whether or not to validate the input
%   'Target'   : ('unit') whether to produce a basis for the 'unit' cube or
%               the nearly 'symmetric' block. 
%
% ---
% MPAWL 1.0, R. Bergmann ~ 2014-08-20

% Transcribed from Mathematica, original code
% localPatternBasis[mM_,False] := Module[{mE,mP,mS,d,dM,j},
% 	{mE,{mP,mS}} = IntegerSmithForm[mM, ExtendedForm-> True];
% 	d = Dimensions[mM][[1]];
% 	dM = patternDimension[mM, validateMatrix -> False];
% 	Return[Table[mS.(1/(Diagonal[mE][[d-dM+j]])*UnitVector[d,d-dM+j]), {j,1,dM}] ]; 
% ]

p = inputParser;
addParamValue(p, 'Validate',true,@(x) islogical(x));
addParamValue(p, 'Target','unit');
parse(p, varargin{:});
ppV = p.Results.Validate;
target = p.Results.Target;
if (ppV)
    isMatrixValid(M);
end
d = size(M,1);
dM = patternDimension(M);
[~,S,V] = snf(M);
V = V*(diag(1./diag(S)));
V = V(:,d-dM+1:d);
for i = 1:dM
    V(:,i) = modM(V(:,i),eye(d),'Target',target,'Validate',false);
end

end

