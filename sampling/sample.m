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
%   'SamplingMethod': ('pointwise') Specify a method to sample, i.e. how
%                      to call f. 'pointwise' will call f for each point
%                      seperately. 'point row' will call f once with all
%                      points rowwise in a matrix, i.e. a matrix det(M)xd,
%                      similarly 'point col' will call f with a dxdet(M)
%                      matrix.
%   'Domain' : ('space') sample domain: sample in 'space' (on P(M)) or 
%               in 'frequency' (w.r.t. FFT, i.e. on G(M^T))
%   'Scale'  : (2*pi in space, 1 in frequency) factor by which the pattern
%              is scaled. Can also be a matrix, e.g. for rotationg all points
%   'Target' : ('symmetric') whether to produce a basis for the 'unit' cube or
%               the nearly 'symmetric' block. 
%   'SampleSize': ([1]) the size of a function evaluations
% ---
% MPAWL, R. Bergmann ~ 2014-09-16
p = inputParser;
addParameter(p, 'Validate',true,@(x) islogical(x));
addParameter(p, 'File','');
addParameter(p, 'SamplingMethod','pointwise');
addParameter(p, 'Target','symmetric');
addParameter(p, 'Domain','space');
addParameter(p, 'Scale',NaN);
addParameter(p, 'SampleSize',[1]);
parse(p, varargin{:});
pp = p.Results;

if (pp.Validate)
    isMatrixValid(M);
end
if numel(pp.File)>0
    if exist(pp.File,'file')
       vars = load(pp.File,'M','v');
       if (all(vars.M==M))
            debug('text',3,'Text',['Sampling values loaded from file ''',pp.File,'''.']);
           v = vars.v;
           return;
       end
       debug('text',3,'Text',['The specified file ''',pp.File,''' does not contain sample values for M, will overwrite them.']);
    else
        debug('text',3,'Text',['The specified file ''',pp.File,''' does not exist yet; trying to write to it']);
    end
end
if (isnan(pp.Scale))
    if strcmp(pp.Domain,'space') && isa(f, 'function_handle')
        pp.Scale = 2*pi;
    else
        pp.Scale = 1;
    end
end
if strcmp(pp.Domain,'frequency')
    M = transpose(M);
end
m = abs(det(M));
d = size(M,1);
Dim = patternDimension(M);
Size = patternSize(M);

if strcmp(pp.Domain,'space')
    Basis = patternBasis(M,'Target',pp.Target);
else
    Basis = generatingSetBasis(M,'Target',pp.Target);
end

if ~isa(f, 'function_handle') % no function handle
    assert(all(size(f)==[m,m]),['The data array has to be of size ',num2str(m),'x',num2str(m),'.']);
end

pointSet = zeros([d,Size]);

for i = 1:Dim
  b = repmat((0:(Size(i)-1)),[d,1]) .* repmat(Basis(:,i),[1,Size(i)]);

  s = ones(1,Dim+1);
  s(1) = d;
  s(i+1) = Size(i);

  b = reshape(b,s);

  s = [d,Size];
  s(i+1) = 1;
  s(1) = 1;

  pointSet = pointSet + repmat(b,s);
end

if strcmp(pp.Domain,'space')
    pointSet = pp.Scale*modM(pointSet,eye(d),'Target',pp.Target);
else
    pointSet = pp.Scale*modM(pointSet,M,'Target',pp.Target);
end
pointSet = reshape(pointSet,d,[]); %row of points

debug('text',2,'Text','Starting to sample');
debug('time',3,'StartTimer','samplingTimer');
if strcmp(pp.SamplingMethod,'pointwise')
  v = zeros([pp.SampleSize,prod(Size)]);
    for index = 1:size(pointSet,2)
        if ~isa(f,'function_handle')
            indA = num2cell( (modM(pointSet(:,index),eye(d),'Target','unit','Validate',false))*m+1);
            v(:,index) = f(indA{:});
        else %function handle
            v(:,index) = f(pointSet(:,index));
        end
    end
else
    if strcmp(pp.SamplingMethod,'point row') %det(M)xd
        v = f(pointSet');
    elseif strcmp(pp.SamplingMethod,'point col')
        v = f(pointSet);
    else
        error('unknown sampling type');
    end
end
debug('time',3,'StopTimer','samplingTimer');

v = reshape(v,[pp.SampleSize,Size]);

if (numel(pp.File) > 0)
    try
        save(pp.File,'M','v');
        debug('text',3,'Text',['Sampling values saved to file ''',pp.File,'''.']);
    catch err
        warning(['Could not save to file ''',pp.File,''', the following error occured: ',err.message]);
    end
end

end

