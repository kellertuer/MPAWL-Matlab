function [decompTree,spaces] = decomposeData2D(gs,Js,M,data,varargin)
% decomposeData2D(gs,Js,M,data) decompose data with respect to certain
% wavelet(s) and dilation matrix/matrices J(s) starting from M with data.
% The first g determines the scaling function to which the data refers.
%
% INPUT
%   gs   : (vector of function handles or strings) function(s) determining
%          the wavelets on the levels
%   Js   : (cells of cells of string(s)) dilation matrix(matrices) strings
%          determining the decomposition patterns (see dilationMatrix2D for
%          all strings). This has to be the same length as gs.
%   M    : pattern(M) the data was sampled on
%   data : sampling values to start the computation with, usually taken as
%          coefficients of the translates.
%
% OUTPUT
%   decompTree : (struct) For each decomposition Level this cell-structure contains
%                one entry scale and one wavelet conteining the coefficents
%                data posesses in the corresponding spaces.
%   spaces     : (struct) Fourier coefficients of all involved functions in order to
%                reconstruct images of wavelet parts, the size
%
% OPTIONAL PARAMETERS
%   'Orthonormalize' : (false) whether or not all translates are normalized
%   'ImageOutput'   : ('None') whether or not to display images of the
%                      corresponding decompositions, 'Wavelet' or 'Scaling'
%                      or 'Both' are the other possible values
%   'ImagePrefix'    : ('') if given, activates the image output into png
%                      files using the Strinf as a prefix
%   'Plotresolution' : (size(data) specify an image resolution larger than
%                      size(data) for the above mentioned images
%   'SpacePrefix'    : ('') if given, activates the saving of the space
%                     Fourier coefficients or tries to load them from there.
%   'Levels'         : (length(gs)) Number of Levels to decompose, which
%                      has to be less than length(gs). Only exception: If g
%                      is only one element and Js are cells of strings,
%                      then these are used on all Levels specified by this
%                      natural number.
% ---
% MPAWL, R. Bergmann ~2014-09-20

% ---
% MPAWL, R. Bergmann ~2014-09-20

p = inputParser;
addParamValue(p, 'ImageOutput','None');
addParamValue(p, 'ImagePrefix','');
addParamValue(p, 'Plotresolution',size(data));
addParamValue(p, 'SpacePrefix','');
addParamValue(p, 'Levels',length(gs));
parse(p, varargin{:});
pp = p.Results;

assert(length(gs)==length(Js),'The amoung of functions/vectors g has to equal the number of sets of matrices in J');

isMatrixValid(M);

d = size(M,1);
dM = patternDimension(M);
epsilon = diag(snf(M)); epsilon = epsilon(d-dM+1:d);

assert(size(data)==epsilon,'Data has to be the size of the patternDimensions of M');

assert( pp.Levels <= length(gs) || length(gs)==1,'Too many Levels');

decompTree(pp.Levels) = struct('scale',0,'wavelet',0);
spaces(pp.Levels) = struct('scaleck',0,'waveletck',0);

%
%

if nargout==1
[decompTree] = recDecomp2Ddata(pp.Levels,gs,Js,M,ckM,data,decompTree,spaces,pp.ImageOutput,pp.ImagePrefix,pp.Plotresolution,pp.SpacePrefix, pp.Orthonormalize);
else
[decompTree,spaces] = recDecomp2Ddata(pp.Levels,gs,Js,M,ckM,data,decompTree,spaces,pp.ImageOutput,pp.ImagePrefix,pp.Plotresolution,pp.SpacePrefix,pp.Orthonormalize);
end

end

%local recursion
function [decompTree, spaces] = recDecomp2Ddata(lvl,g,J,M,ckphi,data,dTree,spaces,ImgOut,ImgPre,PlotRes,SpacePre.Orth)
% pre

% (1) run through all J
lvlmats = J{1};
for i=1:length(lvlmats)
    N = inv(dilationMatrix2D(lvlmats{i}))*M;
    try
        if numel(SpacePre)>0
            StrScale = [SpacePre,'-'lvlmats{i}];
            StrWave = [SpacePre,'-',lvlmats{i}];
            [hatS,hatW] = delaValleePoussinSubspaces(g{1},M,J,'Validate',false,'File',{StrScake,StrWave},'Orthonormalize',Orth);
        else
            [hatS,hatW] = delaValleePoussinSubspaces(g{1},M,J,'Validate',false,'Orthonormalize',Orth);
        end
        
    catch
        warning(['The matrix N,'num2str(N),' is not valid, continuing with next dilationMatrix2D']);
    end
end
% (2) reduce level and call self again

end

