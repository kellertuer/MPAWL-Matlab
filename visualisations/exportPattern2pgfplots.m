function dataStr = exportPattern2pgfplots(varargin)
% exportPattern2pgfplots(mM,f) given data f with respect to the cycles of the
%   patternBasis of the matrix mM, this methods plots a unit plane with the
%   data items represented by parallelograms colord in grayscale of a
%   colormap.
%
% INPUT
%   mM : a matrix indicating the pattern 
%   
%   OPTIONAL PARAMETERS
%       'File'         : (String) Write the generated Code directly into
%                          a file
%       'ExportHeader' : (false) whether to export the Header of the .asy-File
%                    or just the drawing commands.
%
% OUTPUT
%       dataStr : string containing the TikZ code
% ---
% MPAWL ~ R. Bergmann, 2015-11-08
ip = inputParser;
addRequired(ip,'mM');
addParameter(ip,'ExportHeader',false);
addParameter(ip,'File', '');
parse(ip, varargin{:});
vars = ip.Results;

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
        '\\pgfplotsset{compat=1.12}\n',...
        '\\begin{document}%%\n'];
end
dataStr = [dataStr,'\t\\begin{tikzpicture}\n',...
    '\t\t\\begin{axis}[%%\n',...
    '\t\t\txmin=-0.5, xmax=0.5, ymin=-0.5, ymax=0.5,\n',...
    '\t\t\taxis equal image,axis lines=middle,\n\t\tclip=false,\n',...
    '\t\t\txtick={-0.5,-0.25,0,0.25,0.5},\n',...
    '\t\t\txticklabels={$-\\frac{1}{2}$,$-\\frac{1}{4}$,$0$,$\\frac{1}{4}$,$\\frac{1}{2}$},\n',...
    '\t\t\tytick={-0.5,-0.25,0,0.25,0.5},\n',...
    '\t\t\tyticklabels={$-\\frac{1}{2}$,$-\\frac{1}{4}$,$0$,$\\frac{1}{4}$,$\\frac{1}{2}$},\n',...
    '\t\t\tonly marks,mark=*,mark options={scale=1.5, fill=blue!75!black}]\n',...
    '\\addplot+[black,mark options={scale=1.5, black!75}] plot coordinates {'];
    pts = pattern(patternNormalForm(vars.mM));
    for j=1:size(pts,2)
        dataStr=[dataStr,'(',num2str(pts(1,j)),',',num2str(pts(2,j)),')'];
        if j<size(pts,2)
            dataStr = [dataStr,' '];
        end
    end
    dataStr = [dataStr,'};\n'];
endLine = '\t\t\\end{axis}\n\t\\end{tikzpicture}\n';
if vars.ExportHeader
    endLine = [endLine,'\\end{document}'];
end
dataStr = [dataStr,endLine];
if ~isempty(vars.File) %Export given
    fprintf(exportFile,dataStr);
end
if ~isempty(vars.File) %Export given
    fclose(exportFile);
end
if nargout>0
    dataStr = sprintf(dataStr);
else
    clear dataStr;
end
end

