function fh = plotPattern(M)
%plotPattern(M) - plot the pattern of M on 
%   the symmetric unit cube
%
% INPUT
%   M: a regular integer matrix
%
% OUTPUT
%   fh : the fihure handle
    pts = pattern(patternNormalForm(M));
    scatter(pts(1,:),pts(2,:),'filled');
    axis([-.5 .5 -.5 .5])
    if nargout > 0
        fh = fig;
    end
end

