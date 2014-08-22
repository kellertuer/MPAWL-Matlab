function v = generatingSetBasisDecomp(k,M,varargin)
% generatingSetBasisDecomp(k,M)
% Decompose the vector k with respect to the geerating set basis of M, i.e.
% k = modM(V.v,M), where V = geteratingSetBasis(M)
%
% INPUT
%   k : the integer vector (dimension d) that is to be decomposed
%   M : the integer matrix (d times d) for the generating set
%
% OUTPUT
%  v : the basis vector factors, i.e. k = V.v mod M
%
% OPTIONAL PARAMETERS
%   'Validate' : (true) whether or not to validate the input
%   'Target'   : ('unit') whether to produce a basis for the 'unit' cube or
%               the nearly 'symmetric' block. 
%
% ---
% MPAWL 1.0, R. Bergmann ~ 2014-08-22

% Transcribed from the following Mathematica Code
% localgenSetBDecomp[k_,mM_,t_,False] := Module[{mE,mS,mP,d,dM,aBV},
% 	{mE,{mP,mS}} = IntegerSmithForm[Transpose[mM], ExtendedForm-> True];
% 	d = Dimensions[mM][[1]];
% 	dM = patternDimension[mM];
% 	aBV = Transpose[Inverse[mS]];
% 	Return[(modM[Inverse[aBV].k,mE, Target -> t])[[d-dM+1;;d]]];
% ];
p = inputParser;
addParamValue(p, 'Validate',true,@(x) islogical(x));
addParamValue(p, 'Target','unit');
parse(p, varargin{:});
if (p.Results.Validate)
    isMatrixValid(M);
end
assert(isvector(k),'The input k has to be a vector');
if isrow(k)
    vk = k';
else
    vk=k;
end 
d = size(M,1);
dM = patternDimension(M);
[~,S,V] = snf(transpose(M));
v = modM(transpose(V)\vk,S,'Target',p.Results.Target);
v = v(d-dM+1:d);
end