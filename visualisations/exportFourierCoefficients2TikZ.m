function pixelvals = exportFourierCoefficients2TikZ(varargin)
% exportFourierSeries2pgfplots(ckf,samplepoints,File) - export a time space
%   plot of the Fourier series given by the Fourier series of ckf and the
%   origin
%
% INPUT
%   ckf   : a matrix indicating the pattern
%
%   OPTIONAL PARAMETERS
%   'Origin'       : (size(ckf)-1)/2) index indicating the origin in the ckf
%   'File'         : (string) indicating the exportfiles (excluding .tex/.dat)
%   'Data'         : (String) export data to a file of different name
%   'ExportHeader' : (false) whether to export the Header of the .tex-File
%                    or just the drawing commands.
%
% OUTPUT
%       dataStr : string containing the TikZ code
% ---
% MPAWL ~ R. Bergmann, 2015-11-08
ip = inputParser;
addRequired(ip,'ckf');
addParameter(ip,'File','');
addParameter(ip,'origin',NaN);
addParameter(ip,'ExportHeader',false);
addParameter(ip,'Data', '');
parse(ip, varargin{:});
vars = ip.Results;
if isnan(vars.origin)
    origin = (size(vars.ckf)+1)/2;
else
    origin = vars.origin;
end
dims = size(vars.ckf);
[k2,k1] = meshgrid(1:dims(1),1:dims(2)); %k2 constant in y, k1 constant in x-dir
k2 = k2 - origin(2);
k1 = k1 - origin(1);

%
% write output
dataFile = '';
if ~isempty(vars.Data)
    dataFile = vars.Data;
elseif ~isempty(vars.File)
    dataFile = [vars.File,'.dat'];
end
%
if ~isempty(dataFile)
    exportDataFile = fopen(dataFile,'w','n','UTF-8');
% fprintf(exportDataFile,'\tx\ty\tz\n');
for i=1:dims(1)
    for j=1:dims(2)
        if (vars.ckf(i,j)~=0)
        fprintf(exportDataFile,'\t%5.2f\t%5.2f\t%.7e\n',k1(i,j),k2(i,j),vars.ckf(i,j));
        end
    end
end
fclose(exportDataFile);
end
fxmin = min(k1(:)); fxmax = max(k1(:));
fymin = min(k2(:)); fymax = max(k2(:));
fzmin = min(vars.ckf(:)); fzmax = max(vars.ckf(:));
% generate TeX
if ~isempty(vars.File) %Export given
    exportFile = fopen([vars.File,'.tex'],'w','n','UTF-8');
    dataStr = '';
    if vars.ExportHeader
        dataStr = ['%%!TEX TS-options = --shell-escape\n\\documentclass{standalone}\n',...
            '\\usepackage{pgfplots,amsmath}\n',...
            '\\pgfplotsset{compat=1.12}\n',...
            '\\begin{document}%%\n'];
    end
    dataStr = [dataStr,'\t\\pgfplotsset{tick label style={font=\\small}}\n',...
        '\t\\begin{tikzpicture}\n',...
        '\t\t\\begin{axis}[\n\t\t\tview={-45}{35},\n',...
        '\t\t\txmin=',num2str(fxmin),', xmax=',num2str(fxmax),',\n',...
        '\t\t\tymin=',num2str(fymin),', ymax=',num2str(fymax),',\n',...
        '\t\t\tzmin=',num2str(fzmin),', zmax=',num2str(fzmax),',\n',...
		'\t\t\tmesh/ordering=y varies,\n',...
		'\t\t\tenlargelimits={abs=.0125},\n',...
        '\t\t\tonly marks,mark=*,mark options={scale=1.25}]\n',...
    '\\addplot3+[black,mark options={scale=1.5, blue!75!black}] file {',dataFile,'};\n',...
        '\t\t\\end{axis}\n',...
        '\t\\end{tikzpicture}\n'];
    if vars.ExportHeader
        dataStr = [dataStr,'\\end{document}'];
    end
    fprintf(exportFile,dataStr);
    fclose(exportFile);
end
if nargout > 0
    pixelvals = cat(3,k1,k2,vars.ckf);
end
end
