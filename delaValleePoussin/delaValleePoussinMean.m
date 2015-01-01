function [ckphi, BSums] = delaValleePoussinMean(pg,M,varargin)
% [ckphi, BSums] = delaValleePoussinMean(g,M)
% Generate the de la Vallée Poussin mean based on the function g with
% respect to the translates of pattern(M), where g has to be a partition of
% unity in the d-dimensional space, nonnegative everywhere and positive on
% the unit cube. For simplicity g may also be a number or vector, which
% reduces to using the pyramidFunction
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
%                    2*M*unit cube is taken.
%
% ---
% MPAWL, R. Bergmann ~ 2014-09-18
p = inputParser;
addParameter(p, 'Validate',true,@(x) islogical(x));
addParameter(p, 'File',{});
addParameter(p, 'Orthonormalize',true,@(x) islogical(x));
addParameter(p, 'Support',1);
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
    debug('text',3,'Text',['Loading coefficients from file ''',FileStr,''.']);
    if exist(FileStr,'file')
        vars = load(FileStr,'M','ckphi','orthonormalized');
        if (all(vars.M==M))
            ckphi = vars.ckphi;
            if pp.Orthonormalize && ~vars.orthonormalized
                torigin = (size(ckphi)-1)/2;
                ckBSq = bracketSums(ckphi,torigin,M,'Validate',false,'Compute','absolute Squares');
                ckphi = coeffsSpace2Fourier(M,1./(sqrt(ckBSq)),ckphi,torigin,'Validate',false);
            end
        end
        debug('text',3,'Text',['The specified file ''',FileStr,''' does not contain coefficients for M, they will be replaced.']);
    else
        debug('text',3,'Text',['The specified file ''',FileStr,''' does not exist yet.']);
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
        debug('text',3,'Text',['The specified file ''',BSStr,''' does not contain Bracket sums for M, they will be replaced.']);
    else
        debug('text',3,'Text',['The specified file ''',BSStr,''' does not exist yet.']);
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
debug('time',3,'StartTimer','constructing the de la Vallée Poussin mean');
d = size(M,2);
% adM = abs(det(M));
if isvector(pg) && length(pg)==1 && ~isa(pg,'function_handle') 
    g = pg.*ones(1,d);
else
    g = pg;
end
if isa(g,'function_handle')
    ind = 2*(pp.Support);
    tmax = getMaxIndex(transpose(M),'Target','symetric','Cube',ind)+1;
    torigin = tmax+1;
    debug('text',3,'Text','Computing de la Vallée Poussin scaling function');
    ckphi = zeros(2*tmax+1);
    griddims = cell(1,d);
    for i=1:dM
        griddims{i} = 1:(2*tmax(i)+1);
    end
    gridmeshes = cell(1,d);
    [gridmeshes{:}] = ndgrid(griddims{:});
    inds = zeros(dM,numel(ckphi));
    for i=1:d
        inds(i,:) = reshape(gridmeshes{i},1,[]);
    end
    indsc = cell(1,d);
    for i=1:d
        indsc{i} = inds(i,:);
    end
    ckphi(sub2ind(size(ckphi),indsc{:})) = g(transpose(M)\(inds-repmat(torigin,[1,numel(ckphi)])));
elseif isvector(g)
    ind = max(1+2*g,pp.Support*ones(size(g)));
    tmax = getMaxIndex(transpose(M),'Target','symetric','Cube',ind)+1; %+1 for security
    torigin = tmax+1;
    ckphi = zeros(2*tmax+1);
    debug('text',3,'Text','Computing de la Vallée Poussin scaling function');
    griddims = cell(1,d);
    for i=1:d
        griddims{i} = 1:(2*tmax(i)+1);
    end
    gridmeshes = cell(1,d);
    [gridmeshes{:}] = ndgrid(griddims{:});
    inds = zeros(d,numel(ckphi));
    for i=1:d
        inds(i,:) = reshape(gridmeshes{i},1,[]);
    end
    indsc = cell(1,d);
    for i=1:d
        indsc{i} = inds(i,:);
    end
    ckphi(sub2ind(size(ckphi),indsc{:})) = pyramidFunction(g,transpose(M)\(inds-repmat(torigin',[1,numel(ckphi)])));
else
    error('Unknown input type g');
end
debug('time',3,'StopTimer','constructing the de la Vallée Poussin mean');
if pp.Orthonormalize
    torigin = (size(ckphi)-1)/2;
    ckBSq = bracketSums(ckphi,torigin,M,'Validate',false,'Compute','absolute Squares');
    ckphi = coeffsSpace2Fourier(M,1./(sqrt(ckBSq)),ckphi,torigin,'Validate',false);
end
% File savings
if ~isempty(BSStr) || (nargout==2)
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
        orthonormalized = pp.Orthonormalize; %#ok<NASGU>
        save(FileStr,'M','ckphi','orthonormalized');
    catch err
        warning(['Could not save to file ''',FileStr,''', the following error occured: ',err.message]);
    end
end
end