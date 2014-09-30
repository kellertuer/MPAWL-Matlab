function V = generatingSetBasis(M,varargin)
% generatingSetBasis(mM)
% Compute the vectors, whose integral multiples (up to each corresponding
% elementary divisor - 1) reproduce the generating set. The basis vectors
% are returned as columns of a n times patternDimension matrix and are
% orderes w.r.t. the ordering of the elementary divisors on the diagonal.
%
% INPUT
%   M : an integer valued n times n matrix
%
% OUTPUT
%   V : matrix containing the patternDimension(M) basis vectors as columns
%
% OPTIONAL PARAMETERS
%   'Validate' : (true) whether or not to validate the input
%   'Target'   : ('unit') whether to produce a basis for the 'unit' cube or
%               the nearly 'symmetric' block. 
%
% ---
% MPAWL 1.0, R. Bergmann ~ 2014-08-20

% Transcribed from Mathematica
% localGeneratingSetBasis[mM_,t_,False] := Module[{mE,mP,mS,dM,d},
% 	{mE,{mP,mS}} = IntegerSmithForm[Transpose[mM], ExtendedForm -> True];
% 	d = Dimensions[mM][[1]];
% 	dM = patternDimension[mM];
% 	Return[Table[
% 		modM[Transpose[Inverse[mS]].UnitVector[d,d-dM+j],mM,
% 			Target -> t, validateMatrix -> False]
% 		, {j,1,dM}] ];
% ];

p = inputParser;
addParamValue(p, 'Validate',true,@(x) islogical(x));
addParamValue(p, 'Target','symmetric');
parse(p, varargin{:});
ppV = p.Results.Validate;
target = p.Results.Target;
if (ppV)
    isMatrixValid(M);
end
d = size(M,1);
dM = patternDimension(M);
[~,~,V] = snf(transpose(M));
V = transpose(inv(V));
V = V(:,d-dM+1:d);
for i = 1:dM
    V(:,i) = modM(V(:,i),M,'Target',target,'Validate',false,'Index',true);
end
end

