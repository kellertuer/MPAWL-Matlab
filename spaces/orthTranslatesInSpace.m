function hatb = orthTranslatesInSpace(hata, M,J,varargin)
% orthTranslatesInSpace(hata,M,J)
% Let V be a TI space w.r.t orthonormal translates of the function ci wrt
% the shifts y of P(M). Then hata represents the coeffs of a function phi
% a w.r.t these translates. Suppose phi has linear independent translates
% w.r.t x from P(N), N=inv(J)M. Then this function orthonormalizes the
% translates of phi w.r.t P(N).
%
% INPUT
%     hata : FT of coefficients of phi w.r.t translates T(y)xi, y in P(M)
%     M    : integer matrix; basis for the pattern w.r.t of xi
%     J    ; integer matrix indicating the subpattern, i.e. such that
%            inv(J).M is integer
%
% OUTPUT
%     hatb : coeffs of phi_2, which spans the same space of phi, but whose
%            translates are orthonormal
%
% OPTIONAL PARAMETERS
%   'Validate' : (true) Whether to check the validtidy of input M and origin
%
% SEE ALSO
%   Lemma 1.30 (b) in
%   R. Bergmann, Translationsinvariante Räume multivariater anisotroper
%       Funktionen auf dem Torus, Dissertation, University of Lübeck, 2013
%       (in german)
%
% ---
% MPAWL, R. Bergmann 
p = inputParser;
addParamValue(p, 'Validate',true,@(x) islogical(x));
parse(p, varargin{:});
pp = p.Results;
if (pp.Validate)
    isMatrixValid(M);
    isMatrixValid(J);
end
d = size(M,1);
N = round(inv(J)*M); %round for security reasons?
assert(det(M) == det(N)*det(J),'N = inv(J)M hast to be integer valued');
assert(patternDimension(J)==1,['The pattern dimension of J has to be 1 not',num2str(patternDimension(J)),'.']);
dM = patternDimension(M);
dN = patternDimension(N);
debug(text,2,'Orthonormalizing Coefficients...');
epsilon = diag(snf(M)); epsilon = epsilon(d-dM+1:d);
mu = diag(snf(N)); mu = mu(d-dN+1:d);
NTg = transpose(N)*generatingSetBasis(transpose(J));
InvNy = N\patternBasis(patternNormalForm(J));
hN = generatingSetbasis(transpose(N));
lambdag = generatingSetBasisDecomp(NTg,transpose(M),'Target','symmetric','Validate',false);
P = zeros(dM,dN);
for i=1:dN
    P(:,i) = generatingSetBasisDecomp(hN(:,i),transpose(M),'Target','symmetric','Validate',false);
end
hatb = inf(epsilon);
summation = NestedFor(zeros(1,dN),mu-1);
% Timer Debug.
debug('time',3,'StartTimer','orthogonalizeInTranslatestiming');
while summation.hasNext()
    ind = summation.next();
    indM = modM(P*(ind'),diag(epsilon));
    indM2 = modM(P*(ind')+lambdag,diag(epsilon));
    indMc = num2cell(indM);
    indM2c = num2cell(indM2);
    actBSq = abs(hata(indMc))^2 + abs(hata(indM2c))^2;
    assert(actBSq ~= 0,'The translates of phi are not linearly independent w.r.t pattern(N)');
    hatb(indMc) = hata(indMc)*sqrt(abs(det(J))/actBSq);
    hatb(indM2c) = hata(indM2c)*sqrt(abs(det(J))/actBSq);
end
debug('time',3,'StopTimer','orthogonalizeInTranslatestiming');