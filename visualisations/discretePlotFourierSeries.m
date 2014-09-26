function pixelsimg = discretePlotFourierSeries( resolution, coeffs,varargin)
% discretePlotFourierSeries(resolution, coeffs(, origin))
% produce an image or array fom given Fourier coefficients into pixel image
% of size resolution
%
% INPUT
%   resolution : resolution of the image to be computed
%   coefficients : Fourier coefficients
%   origin : (optional) index of the zero coefficients. It is computed as
%            (size(coefficients-1)/2, is not given.
%

p = inputParser;
addOptional('origin',(size(coeffs)-1)/2);
parse(p, varargin{:});
pp = p.Results;

end

