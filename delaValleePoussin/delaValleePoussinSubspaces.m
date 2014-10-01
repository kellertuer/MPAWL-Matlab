function [ hata, hatb ] = delaValleePoussinSubspaces(g,M,J,varargin)
% [hata, hatb] = delaValleePoussinSubspaces(g,M)
% Generate the de la Vallée Poussin mean and its complement with respect to
% an M-invariant space, where the coefficientsbased of the mean are based
% on the function g with respect to the translates of pattern(inv(J).M),
% where g has to be a partition of unity in the d-dimensional space, non-
% negative everywhere and positive on the unit cube. For simplicity g may
% also be a number or vector, which reduces to using the pyramidFunction.
%
% INPUT
%    g : a function or vector characterizing the de la Vallée Poussin mean
%    M : matrix for the translates of the de la Vallée Poussin means
%    J ; integer matrix indicating the subpattern, i.e. such that inv(J).M
%        is integer and |det(J)|=2.
%
% OUTPUT
%   hata : coeffs of the de la Vallée Poussin mean w.r.t. N=inv(J).M
%   hatb : coeffs of the unique orhogonal complement
%
% OPTIONAL PARAMETERS
%   'Validate'     : (true) whether or not to validate the input
%   'File'         : (string) or ({string,string}) whether or not to load
%                    (or if not possible try to save) the coefficients.
%  'Orthonormalize': (true) whether or not to orthonormalize the both
%                    functions w.r.t. their translatestranslates. If both
%                    spaces are succesfully loaded, this option is ignored
% ---
% MPAWL, R. Bergmann ~ 2014-09-26
p = inputParser;
addParamValue(p, 'Validate',true,@(x) islogical(x));
addParamValue(p, 'File',{});
addParamValue(p, 'Orthonormalize',true,@(x) islogical(x));
addParamValue(p, 'Support',1);
parse(p, varargin{:});
pp = p.Results;
if (pp.Validate)
    isMatrixValid(M);
    isMatrixValid(J);
    assert(patternDimension(J)==1,['The pattern dimension of J has to be 1 not',num2str(patternDimension(J)),'.']);
    N = round(inv(J)*M); %#ok<MINV> %round for security reasons?
    isMatrixValid(N);
else
    N = round(inv(J)*M); %#ok<MINV> %round for security reasons?
end    
assert(det(M) == det(N)*det(J),'N = inv(J)M hast to be integer valued');
FiledlVP = '';
FileWav = '';
hata = []; hatb = [];
d = size(M,2);
assert(det(M) == det(N)*det(J),'N = inv(J)M hast to be integer valued');
dM = patternDimension(M);
dN = patternDimension(N);
hM = generatingSetBasis(transpose(M));
epsilon = diag(snf(M)); epsilon = epsilon(d-dM+1:d);
assert( ( isvector(g) || isa(g,'function_handle')),... %neither vector nor function
    'g has to be a vector or value or a function handle');
if ischar(pp.File) %one file indicated -> load dlVP compute Wavelet
    FiledlVP = pp.File; %try loading
end
if iscellstr(pp.File) && (length(pp.File)==2)
    FiledlVP = pp.File{1};
    FileWav = pp.File{2};
end
if ~isempty(FiledlVP)
    debug('text',3,'Text',['Loading coefficients from file',FiledlVP]);
    if exist(FiledlVP,'file')
        vars = load(FiledlVP,'M','J','hata','orthonormalized');
        if ( (all(all(vars.M==M))) &&  (all(all(vars.J==J))) && (all(size(hata)==epsilon')) )
            hata = vars.hata;
            if pp.Orthonormalize && ~varsorthonormalized
                hata = orthogonalizeTranslatesInSpace(hata,M,J,'Validate',false);
            end
        end
        debug('text',3,'Text',['The specified file ''',FiledlVP,''' does not contain coefficients for M w.r.t J, will overwrite them.']);
    else
        debug('text',3,'Text',['The specified file ''',FiledlVP,''' does not exist yet; trying to write to it']);
    end
end
if ~isempty(FileWav)
    debug('text',3,'Text',['Loading coefficients from file',FileWav]);
    if exist(FileWav,'file')
        vars = load(FileWav,'M','J','hatb','orthonormalized');
        if ( (all(all(vars.M==M))) &&  (all(all(vars.J==J))) && (all(size(vars.hatb)==epsilon')) )
            hatb = vars.hatb;
            if pp.Orthonormalize && ~vars.orthonormalized
                hatb = orthogonalizeTranslatesInSpace(hatb,M,J,'Validate',false);
            end
        end
        debug('text',3,'Text',['The specified file ''',FileWav,''' does not contain coefficients for M wrt. J, will overwrite them.']);
    else
        debug('text',3,'Text',['The specified file ''',FileWav,''' does not exist yet; trying to write to it']);
    end
end
NTg = transpose(N)*generatingSetBasis(transpose(J));
InvNy = N\patternBasis(patternNormalForm(J));
hN = generatingSetBasis(transpose(N));
lambdag = round(generatingSetBasisDecomp(NTg,transpose(M),'Target','symmetric','Validate',false));
P = zeros(dM,dN);
for i=1:dN
    P(:,i) = generatingSetBasisDecomp(hN(:,i),transpose(M),'Target','symmetric','Validate',false);
end
P = round(P);
if (numel(hata)>0) && (numel(hatb)>0) %both successfully loaded
    return
elseif (numel(hata)==0) && (numel(hatb)>0) %Compute a from b
    hata = inf(epsilon');
    debug('text',3,'Text','Computing corresponding scaling function as orthogonal compplement of the loaded wavelet function');
    debug('time',3,'StartTimer','computeScalingForthtiming');
    summation = nestedFor(zeros(1,dM),epsilon'-1);
    while(summation.hasNext()) %Can this be done faster?
        ind = summation.next();
        indcp1 = num2cell(ind'+1);
        indshift = modM(ind'+lambdag,diag(epsilon),'Validate',false,'Target','unit','Index',true);
        indshiftcp1 = num2cell(indshift+1);
        hata(indcp1{:}) = exp(-2*pi*1i* (InvNy'*(hM*ind')))*hatb(indshiftcp1{:});
    end
    debug('time',3,'StopTimer','computeScalingForthtiming');
elseif (numel(hata)==0) && (numel(hatb)==0) 
    hata = inf(epsilon');
    summation = nestedFor(zeros(1,dM),epsilon'-1);
    debug('text',3,'Text','Computing de la Vallée Poussin scaling function');
    debug('time',3,'StartTimer','computeScalingFtiming');
    while(summation.hasNext()) %Can this be done faster?
        ind = summation.next();
        indcp1 = num2cell(ind'+1);
        x = transpose(N)\modM( hM * ind',transpose(M),'Validate',false,'Target','symmetric');
        hata(indcp1{:}) = 1/abs(det(J))* BnSum(J,g,x);
    end
    debug('time',3,'StopTimer','computeScalingFtiming');
end %now only hatb might still be 0, but the general wavelet is now computed
hatb = inf(epsilon');
debug('text',3,'Text','Computing corresponding Wavelet function');
debug('time',3,'StartTimer','computeWaveletFtiming');
summation = nestedFor(zeros(1,dM),epsilon'-1);
while(summation.hasNext()) %Can this be done faster?
    ind = summation.next();
    indcp1 = num2cell(ind'+1);
    indshift = modM(ind'+lambdag,diag(epsilon),'Validate',false,'Target','unit','Index',true);
    indshiftcp1 = num2cell(indshift+1);
    hatb(indcp1{:}) = exp(-2*pi*1i* (InvNy'*(hM*(ind'))) ) * hata(indshiftcp1{:});
end
debug('time',3,'StopTimer','computeWaveletFtiming');
if pp.Orthonormalize
    hata = orthTranslatesInSpace(hata,M,J,'Validate',false);
    hatb = orthTranslatesInSpace(hatb,M,J,'Validate',false);
end
% File savings
if ~isempty(FiledlVP)
    try
        orthonormalized = pp.Orthonormalize;
        save(FiledlVP,'M','J','hata','orthonormalized');
    catch err
        warning(['Could not save to file ''',FiledlVP,''', the following error occured: ',err.message]);
    end
end
if ~isempty(FileWav)
    try
        orthonormalized = pp.Orthonormalize;
        save(FileWav,'M','J','hatb','orthonormalized');
    catch err
        warning(['Could not save to file ''',FileWav,''', the following error occured: ',err.message]);
    end
end
end

function s = BnSum(J,g,x)
d = size(J,2);
s=0;
summation = nestedFor(-2*zeros(1,d),2*zeros(1,d));
while(summation.hasNext()) %Can this be done faster?
    p = summation.next()';
    if ~isa(g,'function_handle')
        s = s+pyramidFunction(g,x+transpose(J)*p);
    else %function handle by previous verification
        s = s+g(x+transpose(J)*p);
    end
end
end