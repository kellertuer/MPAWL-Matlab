function T = exportDataOnPattern2vtk(varargin)
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
addParameter(ip,'File', '');
parse(ip, varargin{:});
vars = ip.Results;
dS = 8; %number of digits
d = size(vars.mM,2);
assert(d==3,'This method only works for lattices in 3-D');
dM = patternDimension(vars.mM);
ptsB = patternBasis(vars.mM,'Validate',false,'Target','symmetric');
%
if ~isempty(vars.File) %Export given
    exportFile = fopen(vars.File,'w','n','UTF-8');
end
dataStr = ['# vtk DataFile Version 3.0\n',...
    'Unstructured Grid Example\n',...
    'ASCII\n',...
    'DATASET UNSTRUCTURED_GRID\n'...
];
polyX = [-1/2,  1/2,  1/2, -1/2, -1/2,  1/2, 1/2, -1/2];
polyY = [-1/2, -1/2, -1/2, -1/2,  1/2,  1/2, 1/2,  1/2];
polyZ = [-1/2, -1/2,  1/2,  1/2, -1/2, -1/2, 1/2,  1/2];
polyPts = (vars.mM\([polyX;polyY;polyZ]));
nF = nestedFor(zeros(1,dM),patternSize(vars.mM)-1);
PatchPtsX = zeros(size(polyPts,2),abs(det(vars.mM)));
PatchPtsY = zeros(size(PatchPtsX));
PatchPtsZ = zeros(size(PatchPtsX));
colors = zeros(size(PatchPtsX,2),1);
i=0;
j=0;
adM = abs(det(vars.mM));
while nF.hasNext()
    ind = nF.next().';
    i=i+1;
    tempPts = polyPts + repmat(modM(ptsB*ind,eye(3),'Target','symmetric'),[1,size(polyPts,2)]);
    PatchPtsX(:,i) = tempPts(1,:)';
    PatchPtsY(:,i) = tempPts(2,:)';
    PatchPtsZ(:,i) = tempPts(3,:)';
    indC = num2cell(ind+1);
    colors(i) = vars.f(indC{:});
    if any(abs(tempPts(:))>1/2) % outside -> shift but whereto?
        aInd = any(tempPts<-1/2,2);
        if aInd(1)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)' + 1;
            PatchPtsY(:,adM+j) = tempPts(2,:)';
            PatchPtsZ(:,adM+j) = tempPts(3,:)';
            colors(adM+j) = vars.f(indC{:});    
        end
        if aInd(2)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)';
            PatchPtsY(:,adM+j) = tempPts(2,:)' + 1;
            PatchPtsZ(:,adM+j) = tempPts(3,:)';
            colors(adM+j) = vars.f(indC{:});    
        end
        if aInd(3)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)';
            PatchPtsY(:,adM+j) = tempPts(2,:)';
            PatchPtsZ(:,adM+j) = tempPts(3,:)' + 1;
            colors(adM+j) = vars.f(indC{:});    
        end
        if aInd(1) && aInd(2)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)' + 1;
            PatchPtsY(:,adM+j) = tempPts(2,:)' + 1;
            PatchPtsZ(:,adM+j) = tempPts(3,:)';
            colors(adM+j) = vars.f(indC{:});    
        end
        if aInd(1) && aInd(3)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)' + 1;
            PatchPtsY(:,adM+j) = tempPts(2,:)';
            PatchPtsZ(:,adM+j) = tempPts(3,:)' + 1;
            colors(adM+j) = vars.f(indC{:});    
        end
        if aInd(2) && aInd(3)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)';
            PatchPtsY(:,adM+j) = tempPts(2,:)' + 1;
            PatchPtsZ(:,adM+j) = tempPts(3,:)' + 1;
            colors(adM+j) = vars.f(indC{:});    
        end
        if aInd(1) && aInd(2) && aInd(3)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)' + 1;
            PatchPtsY(:,adM+j) = tempPts(2,:)' + 1;
            PatchPtsZ(:,adM+j) = tempPts(3,:)' + 1;
            colors(adM+j) = vars.f(indC{:});    
        end
        aInd = any(tempPts>1/2,2);
        if aInd(1)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)' - 1;
            PatchPtsY(:,adM+j) = tempPts(2,:)';
            PatchPtsZ(:,adM+j) = tempPts(3,:)';
            colors(adM+j) = vars.f(indC{:});    
        end
        if aInd(2)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)';
            PatchPtsY(:,adM+j) = tempPts(2,:)' - 1;
            PatchPtsZ(:,adM+j) = tempPts(3,:)';
            colors(adM+j) = vars.f(indC{:});    
        end
        if aInd(3)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)';
            PatchPtsY(:,adM+j) = tempPts(2,:)';
            PatchPtsZ(:,adM+j) = tempPts(3,:)' - 1;
            colors(adM+j) = vars.f(indC{:});    
        end
        if aInd(1) && aInd(2)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)' - 1;
            PatchPtsY(:,adM+j) = tempPts(2,:)' - 1;
            PatchPtsZ(:,adM+j) = tempPts(3,:)';
            colors(adM+j) = vars.f(indC{:});    
        end
        if aInd(1) && aInd(3)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)' - 1;
            PatchPtsY(:,adM+j) = tempPts(2,:)';
            PatchPtsZ(:,adM+j) = tempPts(3,:)' - 1;
            colors(adM+j) = vars.f(indC{:});    
        end
        if aInd(2) && aInd(3)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)';
            PatchPtsY(:,adM+j) = tempPts(2,:)' - 1;
            PatchPtsZ(:,adM+j) = tempPts(3,:)' - 1;
            colors(adM+j) = vars.f(indC{:});    
        end
        if aInd(1) && aInd(2) && aInd(3)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)' - 1;
            PatchPtsY(:,adM+j) = tempPts(2,:)' - 1;
            PatchPtsZ(:,adM+j) = tempPts(3,:)' - 1;
            colors(adM+j) = vars.f(indC{:});    
        end
    end
end
dataStr = [dataStr, ...
    'POINTS ',num2str(prod(size(PatchPtsX))),' double\n'];
for i=1:size(PatchPtsX,2)
  for j=1:size(PatchPtsX,1)
    dataStr = [dataStr,num2str(PatchPtsX(j,i),dS),' ',num2str(PatchPtsY(j,i),dS),' ',num2str(PatchPtsZ(j,i),dS),'\n'];
  end
end
dataStr = [dataStr,'CELLS ',num2str(size(PatchPtsX,2)),' ',num2str(9*size(PatchPtsX,2)),'\n'];
for i=1:size(PatchPtsX,2)
    dataStr = [dataStr,'8',...
      ' ',num2str((i-1)*8),...
      ' ',num2str((i-1)*8+1),...
      ' ',num2str((i-1)*8+2),...
      ' ',num2str((i-1)*8+3),...
      ' ',num2str((i-1)*8+4),...
      ' ',num2str((i-1)*8+5),...
      ' ',num2str((i-1)*8+6),...
      ' ',num2str((i-1)*8+7),...
      '\n'];
end
dataStr = [dataStr,'CELL_TYPES ',num2str(size(PatchPtsX,2)),'\n'];
for i=1:size(PatchPtsX,2)
    dataStr = [dataStr,'12\n'];
end
dataStr = [dataStr,'CELL_DATA ',num2str(size(PatchPtsX,2)),'\n'];
dataStr = [dataStr,'SCALARS scalars double\n'];
dataStr = [dataStr,'LOOKUP_TABLE default\n'];
for i=1:size(PatchPtsX,2)
    dataStr = [dataStr,num2str(colors(i)),'\n'];
end

if ~isempty(vars.File) %Export given
    fprintf(exportFile,dataStr);
end
if ~isempty(vars.File) %Export given
    fclose(exportFile);
end
if nargout==0
    clear dataStr;
end
end

