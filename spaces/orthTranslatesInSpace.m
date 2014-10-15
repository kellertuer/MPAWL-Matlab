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
addParameter(p, 'Validate',true,@(x) islogical(x));
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
debug('text',2,'Text','Orthonormalizing Coefficients...');
epsilon = diag(snf(M)); epsilon = epsilon(d-dM+1:d);
NTg = transpose(N)*generatingSetBasis(transpose(J),'Target','unit');
hN = generatingSetBasis(transpose(N));
lambdag = round(generatingSetBasisDecomp(NTg,transpose(M),'Target','symmetric','Validate',false));
P = zeros(dM,dN);
for i=1:dN
    P(:,i) = generatingSetBasisDecomp(hN(:,i),transpose(M),'Target','symmetric','Validate',false);
end
P = round(P);
debug('time',3,'StartTimer','orthogonalization of a functon in space coeffiecients');
griddims = cell(1,dM);
for i=1:dM
    griddims{i} = 1:epsilon(i);
end
gridmeshes = cell(1,dM);
[gridmeshes{:}] = ndgrid(griddims{:});
inds = zeros(dM,prod(epsilon));
for i=1:dM
    inds(i,:) = reshape(gridmeshes{i},1,[]);
end
indMs = modM(P*(inds-1),diag(epsilon),'Index',true)+1;
indMs2 = modM(P*(inds-1) + repmat(lambdag,[1,size(inds,2)]),diag(epsilon),'Index',true)+1;
indMsc = cell(1,d);
indMs2c = cell(1,d);
for i=1:d
    indMsc{i} = indMs(i,:);
    indMs2c{i} = indMs2(i,:);
end
allSq = abs(hata(sub2ind(size(hata),indMsc{:}))).^2 + abs(hata(sub2ind(size(hata),indMs2c{:}))).^2;
assert(all(allSq(:)~=0),'The translates of phi are not linearly independent w.r.t pattern(N)');
hatb = zeros(epsilon');
hatb(sub2ind(size(hatb),indMsc{:})) = 1/abs(det(J)) * hata(sub2ind(size(hata),indMsc{:}))./allSq;
hatb(sub2ind(size(hatb),indMs2c{:})) = 1/abs(det(J)) * hata(sub2ind(size(hata),indMs2c{:}))./allSq;
debug('time',3,'StopTimer','orthogonalization of a functon in space coeffiecients');