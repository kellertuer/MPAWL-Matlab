function T = exportDataOnPattern2mat(varargin)
% exportDataOnPattern2TikZ(mM,f) given data f with respect to the cycles of the
%   patternBasis of the matrix mM, this methods plots a unit plane with the
%   data items represented by parallelograms colord in grayscale of a
%   colormap.
%
% INPUT
%   mM : a matrix indicating the pattern 
%   f  : data on the pattern given as real values between 0 and 1 or they
%   are scaled to that.
%   
%   OPTIONAL PARAMETERS
%       'File'         : (String) Write the generated Code directly into
%                          a file (please provide without extension, you'll
%                          get a tex/tikz and a dat file
%
% OUTPUT
%       dataStr : TikZ code as string
% ---
% MPAWL ~ D. Merkert, 2016-02-17
ip = inputParser;
addRequired(ip,'mM');
addRequired(ip,'f');
addRequired(ip,'File');
parse(ip, varargin{:});
vars = ip.Results;
dM = patternDimension(vars.mM);
ptsB = patternBasis(vars.mM,'Validate',false,'Target','symmetric');
%
f = vars.f;
mM = vars.mM;
sM = patternSize(vars.mM);


fileID = H5F.create(vars.File,'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');

datatypeIDData = H5T.copy('H5T_NATIVE_DOUBLE');
datatypeIDMatrix = H5T.copy('H5T_NATIVE_DOUBLE');
datatypeIDDimension = H5T.copy('H5T_NATIVE_DOUBLE');
datatypeIDBasis = H5T.copy('H5T_NATIVE_DOUBLE');
datatypeIDSize = H5T.copy('H5T_NATIVE_DOUBLE');


dataspaceIDData = H5S.create_simple(length(size(f)), size(f), []);
dataspaceIDMatrix = H5S.create_simple(length(size(mM)), size(mM), []);
dataspaceIDDimension = H5S.create_simple(length(size(dM)), size(dM), []);
dataspaceIDBasis = H5S.create_simple(length(size(ptsB)), size(ptsB), []);
dataspaceIDSize = H5S.create_simple(length(size(sM)), size(sM), []);

datasetIDData = H5D.create(fileID,'data',datatypeIDData,dataspaceIDData,'H5P_DEFAULT');
datasetIDMatrix = H5D.create(fileID,'matrix',datatypeIDMatrix,dataspaceIDMatrix,'H5P_DEFAULT');
datasetIDDimension = H5D.create(fileID,'dimension',datatypeIDDimension,dataspaceIDDimension,'H5P_DEFAULT');
datasetIDBasis = H5D.create(fileID,'basis',datatypeIDBasis,dataspaceIDBasis,'H5P_DEFAULT');
datasetIDSize = H5D.create(fileID,'size',datatypeIDSize,dataspaceIDSize,'H5P_DEFAULT');

H5D.write(datasetIDData,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',permute(f,length(size(f)):-1:1));
H5D.write(datasetIDMatrix,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',permute(mM,length(size(mM)):-1:1));
H5D.write(datasetIDDimension,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',permute(dM,length(size(dM)):-1:1));
H5D.write(datasetIDBasis,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',permute(ptsB,length(size(ptsB)):-1:1));
H5D.write(datasetIDSize,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',permute(sM,length(size(sM)):-1:1));

end

