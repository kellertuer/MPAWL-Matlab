function M = dilationMatrix2D(str)
% dilationMatrix2D(str) return a matrix of the 2D cases based on a string
% from the set {'X','Y','D','Xp','Xm','Yp','Ym'}
%
% INPUT
%   str : a string representing a matrix
%
% OUTPUT
%   M   : corresponding matrix M
%
% ---
% MPAWL, R. Bergmann ~ 2014-09-17
switch str
    case 'X'
        M = [2,0;0,1];
    case 'Y'
        M = [1,0;0,2];
    case 'D'
        M = [1,-1;1,1];
    case 'Xp'
        M = [2,0;1,1];
    case 'Xm'
        M = [2,0;-1,1];
    case 'Yp'
        M = [1,1;0,2];
    case 'Ym'
        M = [1,-1;0,2];
    otherwise
        error('Unknown 2D dilation matrix');
end

