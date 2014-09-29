function pixelsimg = discretePlotFourierSeries( resolution, coefficients)
% discretePlotFourierSeries(resolution, coefficients(, origin))
% produce an image or array fom given Fourier coefficients into pixel image
% of size resolution
%
% INPUT
%   resolution   : resolution of the image to be computed
%   coefficients : Fourier coefficients
%   origin       : (optional) index of the zero coefficients. It is
%                  computed as (size(coefficients-1)/2, is not given.
%
% OUTPUT
%   pixelsimg    : resulting pixel image (or higher dimensional data).
%
% OPTIONAL PARAMETERS
% ---
% MPAWL, R. Bergmann ~2014-09-28

origin = (size(coefficients)-1)/2;

% lengths okay
assert(( (length(resolution)==1)&&(isvector(coefficients)) )...
        ||( length(resolution)==length(size(coefficients)) ),...
    'The coefficients are of different dimension then the resolution');
% resolution large enough
assert(any(size(coefficients)<=resolution),...
    ['The resolution (',num2str(resolution),...
     ' is too small to capture all fequencies of the coefficients (',...
        num2str(size(coefficients)),').']);
    
    debug('text',3,'Text',['Generating an image of size ',resolution,'.']);

    diff = resolution - size(coefficients);
    diffeven = mod(diff,2);
    resorigin = zeros(size(resolution));
    
    % for even entries - share equally, for odd, one more to the left
    diff(diffeven) = diff/2;
    resorigin(diffeven) = origin(diffeven) + diff(diffeven);

    diff(~diffeven) = floor(diff/2);
    resorigin(~diffeven) = origin(~diffeven) + diff(~diffeven) + 1;

    paddright = diff;
    paddleft = diff; paddleft(~diffeven) = paddleft(~diffeven) + 1;
    
    imgck = padarray(coefficients,paddleft,0,'pre');
    imgck = padarray(imgck,paddright,0,'post');
    
    % special version of ifftshift - shift the origin to 1,1
    
    pixelsimg = circshift( fftn(circshift(imgck,-resorigin)), resorigin);
end

