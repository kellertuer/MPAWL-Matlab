function y = pyramidFunction(palpha,x)
% y = pyramidFunction(alpha,x)
% evaluate the d-dimensional analog of the pyramid-function, which is a
% basis for the one-dimensional de la Vallée Poussin scaling function. For
% generality, this function is evaluated w.r.t. the unit cube, i.e. from
% -0.5-alpha to 0.5+alpha, where alpha is a vector or a number (meaning it
% works in any dimension the same) and x is a point or a set of points
%
% INPUT
%   alpha : describing the length of the linear decay, if it is a number,
%           each dimension has the same decay, and each dimension has to be
%           nonnegative and less or equal to 0.5
%       x : a point or a matrix consisting of a point per column
%
% OUTPUT
%       y : result(s) of the pyramid function at the point(s) x
%
% ---
% MPAWL, R. Bergmann ~2014-09-18

if isrow(palpha)
    alpha = palpha';
else
    alpha = palpha;
end
alpha = ones(size(x,1),1).*alpha;
y = ones(1,size(x,2));
for i=1:size(x,1)
y = y.*pyramidfunction1D(alpha(i),x(i,:));
end
end

%local 1D
function y = pyramidfunction1D(alpha,x) %alpha a number, x a row
y = zeros(size(x));
if (alpha==0)
    y(abs(x)<1/2) = 1 * ones(size(y(abs(x)<1/2)));
    y(abs(x)==1/2) = 1/2 * ones(size(y(abs(x)==1/2)));
else
    y( abs(x) < 1/2-alpha) = 1 * ones(size( y( abs(x) < 1/2-alpha) ));
    y( (abs(x)>=1/2-alpha) & (abs(x) <= 1/2+alpha)) = (0.5+alpha-abs(x( (abs(x)>=1/2-alpha) & (abs(x) <= 1/2+alpha) )))./(2*alpha);
    %rest zero per se
end
end