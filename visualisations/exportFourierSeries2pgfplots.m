function pixelvals = exportFourierSeries2pgfplots(varargin)
% exportFourierSeries2pgfplots(ckf,samplepoints,File) - export a time space
%   plot of the Fourier series given by the Fourier series of ckf and the
%   origin
%
% INPUT
%   ckf   : a matrix indicating the pattern
%   dims : number of points in x and y direction
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
addRequired(ip,'dims');
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
ckf = vars.ckf(:);
[k2,k1] = meshgrid(1:dims(1),1:dims(2)); %k2 constant in y, k1 constant in x-dir
k2 = k2(:) - origin(2);
k1 = k1(:) - origin(1);

[x2m,x1m] = meshgrid( ((1:vars.dims(1))-1)/vars.dims(1), ((1:vars.dims(2))-1)/vars.dims(2));
x2m = x2m-0.5;
x1m = x1m-0.5;
x1 = x1m(:)'; x2 = x2m(:)';

% inner product of c_k*exp(2*pi*x^Tk) and summed over k
v = real(sum(...
    repmat(ckf,[1,size(x1,2)]).*...
    exp(-2*pi*1i*...
        (repmat(k1,[1,size(x1,2)]).*repmat(x1,[size(k1,1),1])...
        + repmat(k2,[1,size(x2,2)]).*repmat(x2,[size(k2,1),1]))...
    ),1));
v = reshape(v,vars.dims);
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
for i=1:vars.dims(1)
    for j=1:vars.dims(2)
        fprintf(exportDataFile,'\t%.7e\t%.7e\t%.7e\n',x1m(i,j),x2m(i,j),v(i,j));
    end
    fprintf(exportDataFile,'\n');
end
fclose(exportDataFile);
end
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
        '\t\t\txtick={-0.5,-0.25,0,0.25,0.5},\n',...
		'\t\t\txticklabels={$-\pi$,$-\tfrac{\pi}{2}$,$0$,$-\tfrac{\pi}{2}$,$\pi$},\n',...
		'\t\t\tytick={-0.5,-0.25,0,0.25,0.5},\n',...
		'\t\t\tyticklabels={$-\pi$,$-\tfrac{\pi}{2}$,$0$,$-\tfrac{\pi}{2}$,$\pi$},\n',...
        '\t\t\txmin=-.5,xmax=.5,\n',...
        '\t\t\tymin=-.5,ymax=.5,\n',...
        '\t\t\tzmin=0, zmax=1,\n',...
		'\t\t\tmesh/ordering=y varies,\n',...
		'\t\t\tenlargelimits={abs=.0125}\n\t\t]\n',...
		'\t\t\t\\addplot3[smooth,surf, color=white, opacity=0.5,\n',...
		'\t\t\t\tfill opacity=0.5, faceted color=blue!50!black\n',...
		'\t\t\t] file {',dataFile,'};\n',...
        '\t\t\\end{axis}\n',...
        '\t\\end{tikzpicture}\n'];
    if vars.ExportHeader
        dataStr = [dataStr,'\\end{document}'];
    end
    fprintf(exportFile,dataStr);
    fclose(exportFile);
end
if nargout > 0
    pixelvals = cat(3,x1m,x2m,v);
end
end
