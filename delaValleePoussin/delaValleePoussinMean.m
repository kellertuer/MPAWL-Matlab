function [ckphi, BSums] = delaValleePoussinMean(g,M,varargin)
% delaValleePoussin(g,M)
% Generate the de la Vallée Poussin mean based on the function g with
% respect to the translates of pattern(M), where g has to be a partition of
% unity in the d-dimensional space, nonnegative and positive on the unit
% cube. For simplicity g may also be a number or vector, which reduces to
% using the pyramidFunction
%
% INPUT
%    g : a function or vector characterizing the de la Vallée Poussin mean
%    M : matrix for the translates of the de la Vallée Poussin means
%
% OUTPUT
%   ckphi : Fourier coefficients of the de la Vallée Poussin means, where
%   all dimensions have odd number of entries, and c_0 is at the center
%   BSums : Corresponding Bracket sums
%
% OPTIONAL PARAMETERS
%   'Validate' : (true) whether or not to validate the input
%   'File'     : (string) or ({string,string}) whether or not to load (or
%                if not possible ty try to save) the coefficients. If a
%                second file is given, the same holds for the bracket sums.
%  'Orthonormalize': (true) whether or not to orthonormalize the translates
%  'Support'       : cube indicating the support, if g is a function. for
%                    the vector case, the support is determined
%                    automatically. If it is not given for g, the area
%                    M*unit cube is taken.
%
% ---
% MPAWL, R. Bergmann ~ 2014-09-18
p = inputParser;
addParamValue(p, 'Validate',true,@(x) islogical(x));
addParamValue(p, 'File',{});
addParamValue(p, 'Orthonormalize',true,@(x) islogical(x));
addParamValue(p, 'Support',[]);
parse(p, varargin{:});
pp = p.Results;
if (pp.Validate)
    isMatrixValid(M);
end
BSStr = '';
FileStr = '';
ckphi = [];
if ischar(pp.File)
    FileStr = pp.File; %try loading
end
if iscellstr(pp.File) && (length(pp.File)==2)
    BSStr = pp.File{2};
    FileStr = pp.File{1};
end
if ~isempty(FileStr)
    debug('text',3,'Text',['Loading coefficients from file',FileStr]);
    if exist(FileStr,'file')
        vars = load(FileStr,'M','ckphi');
        if (all(vars.M==M))
            ckphi = vars.ckphi;
            if pp.Orthonormalize
                torigin = (size(ckphi)-1)/2;
                ckBSq = bracketSums(ckphi,torigin,M,'Validate',false,'Compute','absolute Squares');
                ckphi = getFourierFromSpace(1/sqrt(abs(det(M))*ckBSq),ckphi,origin,M,'Validate',false);
            end
        end
        debug('text',3,'Text',['The specified file ''',FileStr,''' does not contain coefficients for M, will overwrite them.']);
    else
        debug('text',3,'Text',['The specified file ''',FileStr,''' does not exist yet; trying to write to it']);
    end
end
if ~isempty(BSStr)
    debug('text',3,'Text',['Loading Bracket sums from file',BSStr]);
    if exist(BSStr,'file')
        vars = load(BSStr,'M','BSums');
        if (all(vars.M==M))
            BSums = vars.BSums;
            if ~isempty(ckphi)
                return;
            end
        end
        debug('text',3,'Text',['The specified file ''',BSStr,''' does not contain Bracket sums for M, will overwrite them.']);
    else
        debug('text',3,'Text',['The specified file ''',BSStr,''' does not exist yet; trying to write to it']);
        if ~isempty(ckphi)
            torigin = (size(ckphi)-1)/2;
            BSums = bracketSums(ckphi,torigin,M,'Validate',false);
            try
                save(BSStr,'M','BSums');
            catch err
                warning(['Could not save to file ''',pp.File,''', the following error occured: ',err.message]);
            end
            return
        end
    end
end
d = size(M,2);
adM = abs(det(M));
if isvector(g)
    ind = 2*max(g,pp.Support);
    tmax = getMaxIndex(M,Target','symetric','Cube',ind);
    torigin = tmax+1;
    debug('text',3,'Text','Computing de la Vallée Poussin scaling function');
    ckphi = zeros(2*tmax+1);
    summation = nestedFor(zeros(1,d),2*tmax+1);
    while(summation.hasnext) %faster?
        ind = summation.next();
        indc = num2cell(ind);
        v = 1;
        p = M\(ind-torigin);
        for j=1:d
            v = v * pyramidFunction(g(j),p(j));
        end
        ckphi(indc{:}) = 1/adM*v;
    end
elseif isa(g,'function handle')
    ind = 2*(pp.Support);
    tmax = getMaxIndex(M,Target','symetric','Cube',ind);
    torigin = tmax+1;
    debug('text',3,'Text','Computing de la Vallée Poussin scaling function');
    ckphi = zeros(2*tmax+1);
    summation = nestedFor(zeros(1,d),2*tmax+1);
    while(summation.hasnext) %faster?
        ind = summation.next();
        indc = num2cell(ind);
        ckphi(indc{:}) = 1/adM*g(M\(ind-torigin));
    end
else
    error('Unknown input type g');
end
if pp.Orthonormalize
    torigin = (size(ckphi)-1)/2;
    ckBSq = bracketSums(ckphi,torigin,M,'Validate',false,'Compute','absolute Squares');
    ckphi = getFourierFromSpace(1/sqrt(abs(det(M))*ckBSq),ckphi,origin,M,'Validate',false);
end
% File savings
if (~isempty(BSStr) || (nargout==2)
    BSums = bracketSums(ckphi,torigin,M,'Validate',false);
    if ~isempty(BSStr)
        try
            save(BSStr,'M','BSums');
        catch err
            warning(['Could not save to file ''',BSStr,''', the following error occured: ',err.message]);
        end
    end
end
if ~isempty(FileStr)
    try
        save(FileStr,'M','ckphi');
    catch err
        warning(['Could not save to file ''',FileStr,''', the following error occured: ',err.message]);
    end
end
end