function T = exportDataOnPattern2PgfPlots(varargin)
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
%       'DataRange'     : (2-vector, [-inf,inf]) data range for the colormap
%                         [-inf,inf] is used to indicate the usual data range given by the data vector itself
%       'Boundary'     : (false) whether or not to draw the boundary of
%                            every parallelotope
%       'File'         : (String) Write the generated Code directly into
%                          a file (please provide without extension, you'll
%                          get a tex/tikz and a dat file
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
addParameter(ip,'Colormap','viridis');
addParameter(ip,'Boundary',false); % frames on elements?
addParameter(ip,'imgSize','8cm'); % frames on elements?
addParameter(ip,'DataRange',[-inf,inf]);
addParameter(ip,'ExportHeader',false);
addParameter(ip,'File', '');
parse(ip, varargin{:});
vars = ip.Results;
d = size(vars.mM,2);
assert(d==2,'This method only works for lattices in the plane');
dM = patternDimension(vars.mM);
ptsB = patternBasis(vars.mM,'Validate',false,'Target','symmetric');
if ~isempty(vars.File) %Export given
    if vars.ExportHeader
        exportFile = fopen([vars.File,'.tex'],'w','n','UTF-8');
    else
        exportFile = fopen([vars.File,'.tikz'],'w','n','UTF-8');
    end
end
dataStr = '';
if vars.ExportHeader
dataStr = ['\\documentclass{standalone}\n',...
        '\\usepackage{pgfplots}\n',...
        '\\usepgfplotslibrary{colormaps}\n',...
        '\\pgfplotsset{compat=1.12}\n',...
        '\\input{viridispgfplots}\n',...
        '\\begin{document}%%\n'];
end    
dataStr = [dataStr,'\t\\begin{tikzpicture}\n',...
    '\t\t\\def\\axisdim{',vars.imgSize,'}\n'];
%
polyX = [-1/2, 1/2, 1/2, -1/2];
polyY = [-1/2, -1/2, 1/2, 1/2];
polyPts = (vars.mM\([polyX;polyY]));
%         \pgfdeclareplotmark{thiscube}{
%          \pgfpathmoveto{\pgfpointdiff{\pgfpointxy{0.23438}{-0.5}}{\pgfpointxy{0}{0}}}
%          \pgfpathlineto{\pgfpointdiff{\pgfpointxy{0.26562}{-0.5}}{\pgfpointxy{0}{0}}}
%          \pgfpathlineto{\pgfpointdiff{\pgfpointxy{-0.23438}{0.5}}{\pgfpointxy{0}{0}}}
%          \pgfpathlineto{\pgfpointdiff{\pgfpointxy{-0.26562}{0.5}}{\pgfpointxy{0}{0}}}
%          \pgfpathclose
%          \pgfusepath{fill}
%        }

lineStr = '\n';
for j=1:size(polyPts,2)
    if j>1
        lineStr = [lineStr,'\t\t\t\\pgfpathlineto']; %#ok<AGROW>
    else
        lineStr = [lineStr,'\t\t\t\\pgfpathmoveto']; %#ok<AGROW>
    end
    lineStr = [lineStr,'{\\pgfpointdiff{\\pgfpointxy{',num2str(polyPts(1,j)),'}{',num2str(polyPts(2,j)),'}}{\\pgfpointxy{0}{0}}}\n']; %#ok<AGROW>
end
    lineStr = [lineStr,'\t\t\t\\pgfpathclose\n\t\t\t\\pgfusepath{fill,stroke}\n\t\t'];
    dataStr = [dataStr,...
    '\t\t\\pgfdeclareplotmark{thiscube}{',lineStr,'}\n',...
    '\t\t\\clip (.5*\\axisdim,.5*\\axisdim) rectangle (1.5*\\axisdim,1.5*\\axisdim);\n',...
    '\t\t\\begin{axis}[colormap/',vars.Colormap,',width=2*\\axisdim,height=2*\\axisdim,scale only axis=true,\n',...
    '\t\t\txmin=-1, xmax=1, ymin=-1, ymax=1,\n',...
    '\t\t\taxis equal,clip=true,xtick=\\empty, ytick=\\empty,hide axis,\n'];
    if (vars.DataRange(1)>-inf) || (vars.DataRange(2)<inf)
        dataStr = [dataStr,'\t\t'];
        if (vars.DataRange(1)>-inf)
            dataStr = [dataStr,'point meta min=',num2str(vars.DataRange(1)),','];
        end
        if (vars.DataRange(1)>-inf) || (vars.DataRange(2)<inf)
            dataStr = [dataStr,'point meta max=',num2str(vars.DataRange(2)),','];
        end
        dataStr = [dataStr,'\n'];
    end
    dataStr = [dataStr,']',...
    '\t\t\\addplot[scatter,scatter src=explicit,,mark=thiscube,only marks,mark options={draw opacity=.5,line width=0.01mm}]\n',...
    '\t\t\ttable[x=X,y=Y,meta=color] {',vars.File,'.dat};\n',...
    '\t\t\\end{axis}\n',...
    '\t\\end{tikzpicture}\n'];
if vars.ExportHeader
    dataStr = [dataStr,'\\end{document}'];
end
if ~isempty(vars.File) %Export given
    fprintf(exportFile,dataStr);
    fclose(exportFile);
end
nF = nestedFor(zeros(1,dM),patternSize(vars.mM)-1);
PatchPtsX = zeros(abs(det(vars.mM)),1);
PatchPtsY = zeros(size(PatchPtsX));
colors = zeros(size(PatchPtsX,2),1);
i=0;
j=0;
adM = abs(det(vars.mM));
while nF.hasNext()
    ind = nF.next().';
    i=i+1;
    tempPt =  modM(ptsB*ind,eye(2),'Target','symmetric');   
    PatchPtsX(i) = tempPt(1);
    PatchPtsY(i) = tempPt(2);
    indC = num2cell(ind+1);
    colors(i) = vars.f(indC{:});
    tempPts = polyPts + repmat(modM(ptsB*ind,eye(2),'Target','symmetric'),[1,size(polyPts,2)]);
    if any(abs(tempPts(:))>1/2) % outside -> shift but whereto?
        aIndBelow = any(tempPts<-1/2,2);
        if aIndBelow(1)
            j=j+1;
            PatchPtsX(adM+j) = tempPt(1)+1;
            PatchPtsY(adM+j) = tempPt(2);
            colors(adM+j) = vars.f(indC{:});    
        end
        if aIndBelow(2)
            j=j+1;
            PatchPtsX(adM+j) = tempPt(1);
            PatchPtsY(adM+j) = tempPt(2) + 1;
            colors(adM+j) = vars.f(indC{:});    
        end
        if aIndBelow(1) && aIndBelow(2)
            j=j+1;
            PatchPtsX(adM+j) = tempPt(1) + 1;
            PatchPtsY(adM+j) = tempPt(2) + 1;
            colors(adM+j) = vars.f(indC{:});    
        end
        aIndAbove = any(tempPts>1/2,2);
        if aIndAbove(1)
            j=j+1;
            PatchPtsX(adM+j) = tempPt(1) - 1;
            PatchPtsY(adM+j) = tempPt(2);
            colors(adM+j) = vars.f(indC{:});    
        end
        if aIndAbove(2)
            j=j+1;
            PatchPtsX(adM+j) = tempPt(1);
            PatchPtsY(adM+j) = tempPt(2) - 1;
            colors(adM+j) = vars.f(indC{:});    
        end
        if aIndAbove(1) && aIndAbove(2)
            j=j+1;
            PatchPtsX(adM+j) = tempPt(1) - 1;
            PatchPtsY(adM+j) = tempPt(2) - 1;
            colors(adM+j) = vars.f(indC{:});    
        end
    end
end
% scale colors (1000 due to colormaps from pgfplots
%colors = round(1000*(colors-min(colors(:)))./(max(colors(:))-min(colors(:))));
% scale to 0,1
%colors = (colors-min(colors(:)))./(max(colors(:))-min(colors(:)));
%
T = table(PatchPtsX(:),PatchPtsY(:),colors(:),'VariableNames',{'X' 'Y' 'color'});
if ~isempty(vars.File) %Export given
    writetable(T,[vars.File,'.dat'],'Delimiter','\t');
end
end