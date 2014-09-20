function hata = coeffsFourier2Space(M,ckf,ckphi,origin, varargin)
% coeffsFourier2Space(ckf,ckphi,origin,M)
%  compute the coefficients of the space of translates (w.r.t. pattern(M))
%  of phi for f based on the given Fourier coefficients ckf and ckphi,
%  if possible (else a NaN-matrix is returned). 
%
%     INPUT
%         M      : matrix indicating the pattern for the translates of phi
%         ckf    : Fourier coefficients of f
%         ckphi  : Fourier coefficitents of phi
%         origin : origin, i.e. the index corr. to c_0 iin both above
%                  parameters ckf and ckphi
%
%    OUTPUT
%        hata    : Fourier transform of the coefficients of the sum of
%                  T(y)phi yielding f, if existent, else an entry of the
%                  result indicates the failure by NaN at the corresponding index.
%
%    OPTIONAL ARGUMENTS
%        'Validate' (true) : whether to validate input or not
%
%     NOTE 
%      The corresponding Mathematica function is called
%      'getFourierfromSpace' and was renamed to fit Matlab conventions
% ---
% MPAWL 1.0, R. Bergmann, 2014-09-10

p = inputParser;
addParamValue(p, 'Validate',true,@(x) islogical(x));
parse(p, varargin{:});
pp = p.Results;
if (pp.Validate)
    isMatrixValid(M);
end
hM = generatingSetBasis(transpose(M),'Target','symmetric','Validate',false);
m = det(M);
d = size(M,1);
dM = patternDimension(M);
epsilon = diag(snf(M));
epsilon = epsilon(d-dM+1:d);
%Compute maximal values
tmax = getMaxIndex(transpose(M));
torigin = tmax+1;
checks = Inf(2*tmax+1);

sumObj = nestedFor(ones(size(size(ckphi))),size(ckphi)); %run through ckspace
while sumObj,hasNext()
    Ind = sumObj.next();
    Indc = num2cell(Ind);
    %INdex in checks array
    checkIndc = num2cell(modM(Ind-origin,transpose(M),'Target','symmetric','Validate',false,'Index',true)+torigin);
    if ckphi(Indc{:}) == 0
        if all(Indc<=size(ckf)) && (ckf(Indc{:}) ~= 0) %lazy inRagne & nonzero
            checks(checkIndc{:}) = NaN; %error
        end %if it is zero, we have no contradiction and everything stays as it is
    else
        if all(Indc<=size(ckf)) %valid index
            actfactor = ckf(Indc{:})/ckphi(Indc{:});
        else
            actfactor=0;
        end
        % Was a checks factor set before?
        if isinf(checks(checkIndc{:})) %No
            checks(checkIndc{:}) = actfactor;
        elseif checks(checkIndc{:}) ~= actfactor %Yes but different than ours
            checks(checkIndc{:}) = NaN; %error, contradiction
        end
    end 
end %end while running over ckphi

%collect result in right order
if (sum(size(epsilon))==2) % one cycle
    hata = zeros(1,epsilon);
else
    hata = zeros(epsilon);
end
hM = generatingSetBasis(transpose(M));
summation = nestedFor(zeros(1,dM),epsilon-ones(1,dM));
while (summation.hasNext())
    ind = summation.next();
    indc = num2cell(ind+1);
    sumInd = num2cell(modM(ind*hM,transpose(M),'Target','symmetric','Index',true)'+torigin);
    hata(indc{:}) = checks(sumInd{:});
end
end