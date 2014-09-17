function v = sample(M,f,varargin)
% v = sample(M,f) sample the function given by f on the pattern (M)
% v = sample(M,A) subsample the det(M)^d array A on the pattern (M)
%
% INPUT
%   M   : a regular integer matrix for the pattern(M) to sample on
%   f/A : a function handle or an det(M)^d, where d=size(M,1). in order
%         to sample f or subsample A.
%
% OUTPUT
%     v : array of sampled values with respect to the cycles of M.
%
% OPTIONAL PARAMETERS
%   'Validate' : (true) Whether to check the validtidy of input M and origin
%   'File' : ([]) specify a file string. If a file exists, this method
%            tries to load from that file. If that fails or the files does
%            not exist, the file is created (or overwritten!) with the sample values.
%   'Sampling Method': ('pointwise') Specify a method to sample, i.e. how
%                      to call f. 'pointwise' will call f for each point
%                      seperately. 'point row' will call f once with all
%                      points rowwise in a matrix, i.e. a matrix det(M)xd,
%                      similarly 'point col' will call f with a dxdet(M)
%                      matrix.
% ---
% MPAWL, R. Bergmann ~ 2014-09-16
p = inputParser;
addParamValue(p, 'Validate',true,@(x) islogical(x));
addParamValue(p, 'File','');
addParamValue(p, 'SamplingMethod','Pointwise');
parse(p, varargin{:});
pp = p.Results;
if (pp.Validate)
    isMatrixValid(M);
    isMatrixValid(J);
end
if ( numel(pp.File)>0) && (strcmp(pp.File(end-3:end),'.mat'))
    if exist(pp.File,'file')
       vars = load(pp.File,'M','v');
       if (all(vars.M==M))
           v = vars.v;
           return;
       end
       debug('text',3,'Text',['The specified file ''',pp.File,''' does not contain sample values for M, will overwrite them.']);
    else
        debug('text',3,'Text',['The specified file ''',pp.File,''' does not exist yet; trying to write to it']);
    end
end
d = size(M,1);
dM = patternDimension(M);
m = abs(det(M));
epsilon = diag(snf(M)); epsilon = epsilon(d-dM+1:d);
pMBasis = patternBasis(M,'Validate',false);
if ~isa(f, 'function_handle') % no function handle
    assert(all(size(f))==[m,m],['The data array has to be of size ',num2str(m),'x',num2str(m),'.']);
end
pointSet = 0;
if (strcmp(pp.SamplingMethod,'pointwise'));
    pointSet = zeros(m,d); %m rows, d cols -> create in row mode
end
summation = nestedFor(zeros(1,dM),epsilon-1);
v = zeros(epsilon);
debug('text',2,'Text','Starting to sample');
debug('time',3,'StartTimer','samplingTimer');
while summation.hasNext()
    ind = summation.next();
    indc = num2cell(ind);
    pt = 2*pi*modM(pMBasis*ind,eye(d),'Target','symmetric','Validate',false);
    if ~isa(f,'function handle')
        indA = num2cell( (modM(pMBasis*ind,eye(d),'Target','symmetric','Validate',false) + 0.5)*m);
        v(indc) = f(indA);
    else %function handle
        if (sum(size(pointSet))==2) % point wise sampling
            v(indc) = f(pt);
        else
            pointSet(sub2ind(epsilon,ind),:) = pt; %#ok<AGROW>
        end
    end
end
if (sum(size(pointSet))~=2)
    if (strcmp(pp.SamplingMethod,'point row'))
        v = f(pointSet);
    elseif (strcmp(pp.SamplingMethod,'point col'))
        v = f(pointSet');
    else
        error('unknown sampling type');
    end
    v = reshape(v,epsilon); %reshape to fit cycles after sampling
end
if (numel(pp.File) > 0)
    try
        save(pp.File,'M','v');
    catch err
        warning(['Could not save to file ''',pp.File,''', the following error occured: ',err.message]);
    end    
end

