function [decompTree,spaces] = decomposeData2D(gs,Js,M,hatdata,varargin)
% decomposeData2D(gs,Js,M,data) decompose data with respect to certain
% wavelet(s) and dilation matrix/matrices J(s) starting from M with data.
% The first g determines the scaling function to which the data refers.
%
% INPUT
%   gs      : (vector of function handles or strings) function(s)
%             determining the scaling functions on the levels
%             0,...,length(Js)
%   Js      : (cells of cells of string(s)) dilation matrix(matrices)
%             strings determining the decomposition patterns (se
%             dilationMatrix2D for all strings). This has to be
%             length(gs)-1, or both equal to length 1.
%   M       : pattern(M) the data was sampled on
%   hatdata : Fourier-Transform of sampled data given with rspect to the
%          translates of the g{1}-delaValléePoussin mean.
%
% OUTPUT
%   decompTree : (struct-tree) For each decomposition Level this cell-structure contains
%                one entry scale and one wavelet conteining the coefficents
%                data posesses in the corresponding spaces.
%   spaces     : (struct-tree) Fourier coefficients of all involved functions of the spaces in order to
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
% MPAWL, R. Bergmann ~2014-09-28, last edit: 2014-09-29

p = inputParser;
addParamValue(p, 'ImageOutput','None');
addParamValue(p, 'ImagePrefix','');
addParamValue(p, 'Plotresolution',size(data));
addParamValue(p, 'SpacePrefix','');
addParamValue(p, 'Levels',length(gs));
parse(p, varargin{:});
pp = p.Results;

assert((length(gs)==length(Js)+1)||(length(gs)==1&&length(Js)==1,'The amoung of functions/vectors g has to be one more than sets of matrices in J');

isMatrixValid(M);

d = size(M,1);
dM = patternDimension(M);
epsilon = diag(snf(M)); epsilon = epsilon(d-dM+1:d);

assert(size(hatdata)==epsilon,'Data has to be the size of the patternDimensions of M');

assert( pp.Levels-1 <= length(gs) || length(gs)==1,'Too many Levels');

decompTree(pp.Levels) = struct('scale',0,'wavelet',0);
spaces(pp.Levels) = struct('scaleck',0,'waveletck',0);

%
%
if ~strcmp(pp.ImageOutput,'None')
    if numel(pp.SpacePrefix)>0
        fileMStr = [SpacePrefix,'-M.dat'];
    else
        fileMStr='';
    end
    ckM = deLaValleePoussinMean(gs{1},M,'Validate',false,'Orthonormalize',pp.Orthonormalize,File,fileMStr);
else
    ckM = [];
end
cellnames = Js{1};
cellnames{end+1} = 'Scale';
decompTree = cell2struct(cell(size(cellnames)),cellnames,2);
if nargout==1
    [decompTree] = recDecomp2Ddata(pp.Levels,gs{2:end},Js,M,ckM,data,decompTree,spaces,pp.ImageOutput,pp.ImagePrefix,pp.Plotresolution,pp.SpacePrefix, pp.Orthonormalize);
else
    spaces(pp.Levels) = cell2struct(cell(size(cellnames)),cellnames,2);
    [decompTree,spaces] = recDecomp2Ddata(pp.Levels,gs{2:end},Js,M,ckM,data,decompTree,spaces,pp.ImageOutput,pp.ImagePrefix,pp.Plotresolution,pp.SpacePrefix,pp.Orthonormalize);
end

end

%local recursion
function [decompTree, decompSpaces] = recDecomp2Ddata(lvl,g,J,M,ckphiM,hatdata,dTree,spaces,ImgOut,ImgPre,PlotRes,SpacePre.Orth)
% (1) run through all J
lvlmats = J{1};
decompTree = dTree;
decompSpaces = spaces;
for i=1:length(lvlmats)
    N = inv(dilationMatrix2D(lvlmats{i}))*M; %#ok<MINV>
    if numel(ckphiM)>0
        lOrigin = (size(ckphiM)+1)/2;
    end
    try
        if numel(SpacePre)>0
            StrScale = [SpacePre,'-',lvlmats{i}];
            StrWave = [SpacePre,'-',lvlmats{i}];
            [hatS,hatW] = delaValleePoussinSubspaces(g{1},M,dilationMatrix2D(lvlmats{i}),'Validate',false,'File',{StrScale,StrWave},'Orthonormalize',Orth);
        else
            [hatS,hatW] = delaValleePoussinSubspaces(g{1},M,dilationMatrix2D(lvlmats{i}),'Validate',false,'Orthonormalize',Orth);
        end
        constructionSuc = true;
    catch
        warning(['The matrix N,',num2str(N),' is not valid, continuing with next dilationMatrix2D']);
        constructionSuc = false;
    end
    if constructionSuc %succesfully constructed/loaded Wavelets/Scaling functions -> decompose and recursively continue
        [hatDataSubS,hatDataSubW] = patternFWT(M,dilationMatrix2D(lvlmats{i}),hatdata,hatS,hatW,'Validate',false);
        % Output Images
        if strcmp(ImgOut,'Wavelet') || strcmp(ImgOut,'Both') || numel(spaces)>0
            ckWav = coeffsSpace2Fourier(M,hatW,ckphiM,lOrigin);
        end
        if strcmp(ImgOut,'Wavelet') || strcmp(ImgOut,'Both')
            ckDataPart = coeffsSpace2Fourier(N,hatDataSubW,ckWav,lOrigin);
            figure;
            img = discretePlotFourierSeries(PlotRes,ckDataPart);
            imagesc(img); colormap(rwb);
            title(['Wavelet Part ',ImgPre,'-',lvlmats{i}]);
            if numel(ImgPre)>0 %Save image?
                imgmin = min(min(img));
                imgmax = max(max(img));
                map = colormap(rwb(256));
                frs = 255*(img-imgmin)/(imgmax-imgmin);
                imwrite(frs,map,[ImgPre,'-',lvlmats{i},'-WaveletPart.png'],'png');
            end
        end
        if (strcmp(ImgOut,'Scale') || strcmp(ImgOut,'Both')) || numel(ckphiM)>0 || numel(spaces)>0 %any need for ckScale?
            ckScal = coeffsSpace2Fourier(M,hatW,ckphiM,lOrigin);
        end
        if strcmp(ImgOut,'Scale') || strcmp(ImgOut,'Both')
            ckDataPart = coeffsSpace2Fourier(N,hatDataSubS,ckScal,lOrigin);
            figure;
            img = discretePlotFourierSeries(PlotRes,ckDataPart);
            imagesc(img); colormap(rwb);
            title(['Scaling Part ',ImgPre,'-',lvlmats{i}]);
            if numel(ImgPre)>0 %Save image?
                imgmin = min(min(img));
                imgmax = max(max(img));
                map = colormap(rwb(256));
                frs = 255*(img-imgmin)/(imgmax-imgmin);
                imwrite(frs,map,[ImgPre,'-',lvlmats{i},'-ScalePart.png'],'png');
            end
        end
        % Recursion
        if (lvl>1)
            % Many ifs for each of the optional parameters
            % Fourier Coefficients of Mother space
            if numel(ckphiM)>0
                ckN = ckScal;
            else
                ckN = [];
            end
            % Extend DecompTree
            % ImgPre
            if numel(ImgPre)>0
                ImgNPre = [ImgPre,'-',lvlmats{i}];
            else
                ImgNPre = '';
            end
            newcells = J{2};
            newcells{end+1} = 'Scale';
            newcells{end+1} = 'Wavelet';
            if numel(ImgPre)>0
                SpaceNPre = [SpacePre,'-',lvlmats{i}];
            else
                SpaceNPre = '';
            end
            if numel(spaces) ==0
                subTree = recDecomp2Ddata(lvl-1,g{2:end},J{2:end},N,ckN,hatDataSubS,cell2struct(cell(size(J{2})),J{2},2),[],ImgNPre,PlotRes,SpaceNPre,Orth);
            else
                [subTree,subSpaces] = recDecomp2Ddata(lvl-1,g{2:end},J{2:end},N,ckN,hatDataSubS,cell2struct(cell(size(Js{2})),Js{2},2), cell2struct(cell(size(Js{2})),Js{2},2),ImgNPre,PlotRes,SpaceNPre,Orth);
            end
            decompTree(lvlmats{i}) = subTree;
            if (numel(spaces)>0)
                decompSpaces(lvlmats{i}) = subSpaces;
            end
        else % lvl = 1 - end of recursion, generate reasonable end-structs
            decompTree(lvlmats{i}) = struct('Scale',hatDataSubS,'Wavelet',hatDataSubW);
            if numel(spaces) > 0
                decompSpaces(lvlmats{i}) = struct('Scale',ckScale,'Wavelet',ckWavelet);
            end
        end
    else %if the construction was not succesfull, the tree is just initialized as empty
        decompTree(lvlmats{i}) = struct('Scale',[],'Wavelet',[]);
        if numel(spaces) > 0
            decompSpaces(lvlmats{i}) = struct('Scale',[],'Wavelet',[]);
        end
    end %end of succesfull construction if.
end %end running over all Js of this level (levelmats)
end