function cmap = rwb(n)
% rwb(n) produce the red-white-blue colormap
% INPUT
%    n : (256) colormap length, if not provided, 256 is used
%
% OUTPUT
%   cmap : the colormap
%
% ---
% MPAWL, R. Bergmann, 2014-09-29

if nargin==0
    n = 256;
end
cmap = [0:1/(ceil(n/2)-1):1, ones(1,floor(n/2));...
        0:1/(ceil(n/2)-1):1, (1-mod(n,2)*1/(floor(n/2))):-1/(ceil(n/2)-1):0; ...
        ones(1,ceil(n/2)), (1-mod(n,2)*1/(floor(n/2))):-1/(ceil(n/2)-1):0]';
% cmap = [ones(1,128), 1:-1/127:0; 0:1/127:1, 1:-1/127:0; 0:1/127:1,ones(1,128)]';
end

