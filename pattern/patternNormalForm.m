function pMf = patternNormalForm(M,varargin)
% patternNormalform(M) computes the pattern normal form of M
%
%   The pattern normal form of the matrix, i.e. a matrix obtained by
%   computing the gaussian elimination (in integers) to obtain an upper
%   triangular matrix, where all elements above the diagonal are smaller
%   than the diagonal element
%   
%   INPUT
%       M: a regular integral dxd dimensional matrix
%
%   OUTPUT
%       pMf: The matrix in pattern normal form
%
% ---
% MPAWL 1.0, R. Bergmann ~ 2013-09-17 ~ last update: 2013-11-15
    p = inputParser;
    addParamValue(p, 'Validate',true,@(x) islogical(x));
    parse(p, varargin{:});
    optionals = p.Results;

    if (optionals.Validate)
        isMatrixValid(M);
    end
    pMf = M;
    d = size(M,1);
    % form upper triangular matrix
    for col = 1:(d-1)
        for row = col+1:d
            pMf = gcdonrows(pMf,col,row);
        end
    end
    % get diagonal positive
    for row=1:d
        if pMf(row,row)<0 
            pMf(row,:) = -pMf(row,:);
        end
    end
    % get upper nonzero values of a colum to lie between 0 and
    for col=2:d
        for row=col-1:-1:1
            f=0;
            if (pMf(row,col)<0) || (pMf(row,col)>=pMf(col,col))
                f = -floor(pMf(row,col)/pMf(col,col));
            end
            pMf(row,:) = pMf(row,:) + f*pMf(col,:);
        end
    end
end

% local functions
function rM = gcdonrows(M,ci,ri)
% gcdonRows(M,rowindex,colindex) perform a GCD on complete rows
% ---
% MPAWL 1.0 written on 2013-09-11 by Ronny Bergmann

% Just transcribed from Mathematica...
% TODO: Can that be optimized in MatLab by using Matrix-Vector notation?

rM = M; % computation matrix
if rM(ri, ci) ~= 0
    % Modify by row addition, such that ci,ci is nonnegative
    if rM(ci,ci) < 0
        rM(ci,:) = rM(ci,:) - fix(rM(ci,ci)/rM(ri,ci))*rM(ri,:);
    end
    % get it positive
    if rM(ci,ci) == 0
        rM(ci,:) = rM(ci,:) + sign(rM(ri,ci))*rM(ri,:);
    end
    % get the second entry in that column also positive
    if rM(ri,ci) < 0
        f = ceil(rM(ri,ci)/rM(ci,ci));
        if f==0
            f=f+1;
        end
        rM(ri,:) = rM(ri,:) - f*rM(ci,:);
    end
    % Euclidean algorithm on rows in order to get (ri,ci) to zero
    while rM(ri,ci) ~= 0
        if abs(rM(ci,ci)) > abs(rM(ri,ci))
            f = floor(rM(ci,ci)/rM(ri,ci));
            if mod(rM(ci,ci),rM(ri,ci))==0
                f = f-sign(rM(ri,ci))*sign(rM(ci,ci));
            end
            rM(ci,:) = rM(ci,:)-f*rM(ri,:);
        else
            f = floor(rM(ri,ci)/rM(ci,ci));
            rM(ri,:) = rM(ri,:) - f*rM(ci,:);
        end
    end
end
end