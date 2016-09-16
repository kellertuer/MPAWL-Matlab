function dataStr = exportDataOnPattern2TikZ(varargin)
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
%       'Colormap'     : (String, 'parula') colormap to volor data items
%                       (see pgfplots documentation Ch. 5 Sec. 2 colormaps)
%       'Boundary'     : (false) whether or not to draw the boundary of
%                            every parallelotope
%       'File'         : (String) Write the generated Code directly into
%                          a file
%       'ExportHeader' : (false) whether to export the Header of the .asy-File
%                    or just the drawing commands.
%
% OUTPUT
%       dataStr : TikZ code as string
% ---
% MPAWL ~ R. Bergmann, 2015-11-08
ip = inputParser;
addRequired(ip,'mM');
addRequired(ip,'f');
addParameter(ip,'Colormap','parula');
addParameter(ip,'Boundary',false); % frames on elements?
addParameter(ip,'ExportHeader',false);
addParameter(ip,'File', '');
parse(ip, varargin{:});
vars = ip.Results;
dS = 2; %number of digits
d = size(vars.mM,2);
assert(d==2,'This method only works for lattices in the plane');
dM = patternDimension(vars.mM);
ptsB = patternBasis(vars.mM,'Validate',false,'Target','symmetric');
%
if ~isempty(vars.File) %Export given
    exportFile = fopen(vars.File,'w','n','UTF-8');
end
dataStr = '';
if vars.ExportHeader
    dataStr = ['\\documentclass{standalone}\n',...
        '\\usepackage{pgfplots}\n',...
        '\\usepgfplotslibrary{colormaps}\n',...
        '\\pgfplotsset{compat=1.12}\n',...
        '\\tikzset{BoxC/.style={/utils/exec={\\pgfplotscolormapdefinemappedcolor{#1}},%%\n',...
        '\tdraw=mapped color!80!black,very thin, fill=mapped color!80!white}}\n',...
        '\\begin{document}%%\n'];
end
dataStr = [dataStr,'\t\\begin{tikzpicture}\n',...
    '\t\t\\begin{axis}[colormap/',vars.Colormap,',\n,',...
    '\t\t\txmin=-0.5, xmax=0.5, ymin=-0.5, ymax=0.5,\n',...
    '\t\t\taxis equal,clip=true,xtick=\\empty, ytick=\\empty,hide axis]\n',...
    '\\pgfplotsextra{\\clip (axis cs:-0.5,-0.5) rectangle (axis cs:0.5,0.5);}\n'];
if ~isempty(vars.File) %Export given
    fprintf(exportFile,dataStr);
end
polyX = [-1/2, 1/2, 1/2, -1/2];
polyY = [-1/2, -1/2, 1/2, 1/2];
polyPts = (vars.mM\([polyX;polyY]));
nF = nestedFor(zeros(1,dM),patternSize(vars.mM)-1);
PatchPtsX = zeros(size(polyPts,2),abs(det(vars.mM)));
PatchPtsY = zeros(size(PatchPtsX));
colors = zeros(size(PatchPtsX,2),1);
i=0;
j=0;
adM = abs(det(vars.mM));
while nF.hasNext()
    ind = nF.next().';
    i=i+1;
    tempPts = polyPts + repmat(modM(ptsB*ind,eye(2),'Target','symmetric'),[1,size(polyPts,2)]);
    PatchPtsX(:,i) = tempPts(1,:)';
    PatchPtsY(:,i) = tempPts(2,:)';
    indC = num2cell(ind+1);
    colors(i) = vars.f(indC{:});
    if any(abs(tempPts(:))>1/2) % outside -> shift but whereto?
        aInd = any(tempPts<-1/2,2);
        if aInd(1)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)' + 1;
            PatchPtsY(:,adM+j) = tempPts(2,:)';
            colors(adM+j) = vars.f(indC{:});    
        end
        if aInd(2)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)';
            PatchPtsY(:,adM+j) = tempPts(2,:)' + 1;
            colors(adM+j) = vars.f(indC{:});    
        end
        if aInd(1) && aInd(2)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)' + 1;
            PatchPtsY(:,adM+j) = tempPts(2,:)' + 1;
            colors(adM+j) = vars.f(indC{:});    
        end
        aInd = any(tempPts>1/2,2);
        if aInd(1)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)' - 1;
            PatchPtsY(:,adM+j) = tempPts(2,:)';
            colors(adM+j) = vars.f(indC{:});    
        end
        if aInd(2)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)';
            PatchPtsY(:,adM+j) = tempPts(2,:)' - 1;
            colors(adM+j) = vars.f(indC{:});    
        end
        if aInd(1) && aInd(2)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)' - 1;
            PatchPtsY(:,adM+j) = tempPts(2,:)' - 1;
            colors(adM+j) = vars.f(indC{:});    
        end
    end
end
% scale colors (1000 due to colormaps from pgfplots
colors = round(1000*(colors-min(colors(:)))./(max(colors(:))-min(colors(:))));
for i=1:size(PatchPtsX,2)
    line = ['\t\t\t\\draw[BoxC=',num2str(colors(i)),']'];
  for j=1:size(PatchPtsX,1)
    if j>1
        line = [line,' -- '];
    end
    line = [line,'(axis cs:',num2str(PatchPtsX(j,i),dS),',',num2str(PatchPtsY(j,i),dS),')'];
  end
  line = [line,' -- cycle;\n'];
  if ~isempty(vars.File) %Export given
      fprintf(exportFile,line);
  end
  if nargout>0
      dataStr = sprintf([dataStr,line]);
  end
end
% Export end of file lines
endLine = '\t\t\\end{axis}\n\t\\end{tikzpicture}\n';
if vars.ExportHeader
    endLine = [endLine,'\\end{document}'];
end
if ~isempty(vars.File) %Export given
    fprintf(exportFile,endLine);
end
if nargout>0
    dataStr = sprintf([dataStr,endLine]);
end
if ~isempty(vars.File) %Export given
    fclose(exportFile);
end
if nargout==0
    clear dataStr;
end
end

