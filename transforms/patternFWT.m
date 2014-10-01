function [hatf_s,hatf_w] = patternFWT(M,J,phata,phatbS,phatbW,varargin)
% [hatf_s,hatf_w] = patternFWT(M,J,hata, hatbS, hatbW)
% [hatf_s,hatf_w] = patternFWT(M,J,ckf, ckphi, ckpsi, ckxi, origin)
% compute the wavelet transform of a function f with respect to a scaling
% space Xi (w.r.t. pattern(M)) into phi/psi w.r.t. N=inv(J)M. Though there
% are two variants for the input
%
% INPUT
%   M : matrix indicating the translates pattern(M) of the Xi-translates
%   J : matrix indication the pattern decimation, i.e. pattern(inv(J)M) is
%   the pattern for phi/psi
%
% Variant (a) all coefficients are given as Fourier transform of a
% coefficient set w.r.t. to the translates of Xi:
%   hata   : Fourier transform of the coeffs of f
%   hatbS  : Fourier transform of the coeffs of phi
%   hatbW  : Fourier transform of the coeffs of psi
% 
% Variant (b) all coefficients are Fourier coefficients arrays of same size
%   ckf     : Fourier coefficients of f
%   ckphi   : Fourier coefficients of the smaller scaling function phi
%   ckpsi   : Fourier coefficients of the wavelet function psi
%   ckXi    : Fourier coefficients of the mother scaling function Xi
%   origin  : indicates the c_0 index in all four arrays.
%
% OUTPUT
%     hata_s : coefficients of the low pass part of f, i.e. w.r.t. phi 
%     hata_w : coefficients of the high pass part of f, i.e. w.r.t. psi
% These can easily be transformed into Fourier coefficients using
% coeffsSpace2Fourier.
%
% OPTIONAL PARAMETERS
%   'Validate' : (true) Whether to check the validtidy of input M and origin
%
% ---
% MPAWL, R. Bergmann ~2014-09-20

p = inputParser;
addOptional(p,'ckxi',[]);
addOptional(p,'origin',[]);
addParamValue(p, 'Validate',true,@(x) islogical(x));
parse(p, varargin{:});
pp = p.Results;
if (pp.Validate)
    isMatrixValid(M);
    isMatrixValid(J);
end
    N = inv(J)*M; %#ok<MINV>
if (pp.Validate)
    isMatrixValid(N);
end
d = size(M,1);
dM = patternDimension(M);
epsilon = diag(snf(M));
epsilon = epsilon(d-dM+1:d);
if numel(pp.origin)*numel(pp.ckxi>0)
    if (pp.Validate)
        assert(all(size(pp.ckxi)==size(phata)),['Either ck(xi) (',num2str(size(ckxi)),') or ck(f) (',num2str(size(phata)),') is of wrong size']);
        assert(all(size(pp.ckxi)==size(phatbS)),['Either ck(xi) (',num2str(size(ckxi)),')or ck(phi) (',num2str(size(phatbS)),')is of wrong size']);
        assert(all(size(pp.ckxi)==size(phatbW)),['Either ck(xi) (',num2str(size(ckxi)),')or ck(psi) (',num2str(size(phatbW)),')is of wrong size']);
        assert(all(origin<=size(pp.ckxi)),'Origin of wring dimension or out of range for Fourier coeffcients of xi.');
    end
    hata = coeffsFourier2Space(M,phata,ckxi,pp.origin);
    hatbS = coeffsFourier2Space(M,phatbS,ckxi,pp.origin);
    hatbW = coeffsFourier2Space(M,phatbW,ckxi,pp.origin);
else
    hata = phata;
    hatbS = phatbS;
    hatbW = phatbW;
    if (pp.Validate)
        assert(all(size(hata)==(epsilon')),['Coefficients of f are of dimension ',num2str(size(hata)),' but M requires ',num2str(epsilon'),'.']);
        assert(all(size(hatbS)==(epsilon')),['Coefficients of f are of dimension ',num2str(size(hatbS)),' but M requires ',num2str(epsilon'),'.']);
        assert(all(size(hatbW)==(epsilon')),['Coefficients of f are of dimension ',num2str(size(hatbW)),' but M requires ',num2str(epsilon'),'.']);
    end
end
debug('text',2,'Text','Performing the Wavelet Transform');
% Further parameters
dN = patternDimension(N);
mu = diag(snf(N)); mu = mu(d-dN+1:d);
NTg = transpose(N)*generatingSetBasis(transpose(J));
hN = generatingSetBasis(transpose(N));
lambdag = generatingSetBasisDecomp(NTg,transpose(M),'Target','symmetric','Validate',false);
P = zeros(dM,dN);
for i=1:dN
    P(:,i) = generatingSetBasisDecomp(hN(:,i),transpose(M),'Target','symmetric','Validate',false);
end
hatf_s = zeros(mu');
hatf_w = zeros(mu');
summation = nestedFor(zeros(1,dN),mu'-1);
% Timer Debug.
debug('time',3,'StartTimer','patternWavelet-Transform');
while summation.hasNext()
    ind = summation.next();
    indcp1 = num2cell(ind'+1);
    indM = modM(P*(ind'),diag(epsilon));
    indM2 = modM(P*(ind')+lambdag,diag(epsilon));
    indMcp1 = num2cell(indM+1);
    indM2cp1 = num2cell(indM2+1);
    hatf_s(indcp1{:}) = 1/abs(det(J)) * ( ...
        conj(hatbS(indMcp1{:}))*hata(indMcp1{:}) + conj(hatbS(indM2cp1{:}))*hata(indM2cp1{:})...
        );
    hatf_w(indcp1{:}) = 1/abs(det(J)) * ( ...
        conj(hatbW(indMcp1{:}))*hata(indMcp1{:}) + conj(hatbW(indM2cp1{:}))*hata(indM2cp1{:})...
        );
end
debug('time',3,'StopTimer','patternWavelet-Transform');
end