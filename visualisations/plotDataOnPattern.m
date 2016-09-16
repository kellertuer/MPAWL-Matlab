function fh = plotDataOnPattern(varargin)
% plotDataOnPattern(mM,f) given data f with respect to the cycles of the
%   patternBasis of the matrix mM, this methods plots a unit plane with the
%   data items represented by parallelograms colord in grayscale of a
%   colormap.
%
% INPUT
%   mM : a matrix indicating the pattern 
%   f  : data on the pattern given as real values between 0 and 1 or they
%   are scaled to that.
%   
%   OPTIONAL PARAMETERS
%       'Colormap'       : (gray) colormap to volor data items
%       'DataRange'      : ([-inf,inf]) cut data Range to certain range
%       'Boundary'       : (false) whether or not to draw the boundary of
%                            every parallelotope
%
% OUTPUT
%       fh : figure handle of the resulting figure
% ---
% MPAWL ~ R. Bergmann, 2015-11-08
p = inputParser;
addRequired(p,'mM');
addRequired(p,'f');
addParameter(p,'DataRange',[-inf,inf]);
addParameter(p,'Colormap',gray(256));
addParameter(p,'Boundary',false); % frames on elements?
parse(p, varargin{:});
vars = p.Results;

d = size(vars.mM,2);
assert(d==2,'This method only works for lattices in the plane');
dM = patternDimension(vars.mM);
ptsB = patternBasis(vars.mM,'Validate',false,'Target','symmetric');
polyX = [-1/2, 1/2, 1/2, -1/2];
polyY = [-1/2, -1/2, 1/2, 1/2];
polyPts = (vars.mM\([polyX;polyY]));
nF = nestedFor(zeros(1,dM),patternSize(vars.mM)-1);
PatchPtsX = zeros(size(polyPts,2),round(abs(det(vars.mM))));
PatchPtsY = zeros(size(PatchPtsX));
colors = zeros(size(PatchPtsX,2),1);
i=0;
j=0;
adM = int32(abs(det(vars.mM)));
f = vars.f;
if isfinite(vars.DataRange(1)) % lower bound given
    f = max(vars.f,vars.DataRange(1));
end
if isfinite(vars.DataRange(2)) % upper bound given
    f = min(f,vars.DataRange(2));
end
while nF.hasNext()
    ind = nF.next().';
    i=i+1;
    tempPts = polyPts + repmat(modM(ptsB*ind,eye(2),'Target','symmetric'),[1,size(polyPts,2)]);
    PatchPtsX(:,i) = tempPts(1,:)';
    PatchPtsY(:,i) = tempPts(2,:)';
    indC = num2cell(ind+1);
    colors(i) = f(indC{:});
    if any(abs(tempPts(:))>1/2) % outside -> shift but whereto?
        aInd = any(tempPts<-1/2,2);
        if aInd(1)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)' + 1;
            PatchPtsY(:,adM+j) = tempPts(2,:)';
            colors(adM+j) = f(indC{:});    
        end
        if aInd(2)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)';
            PatchPtsY(:,adM+j) = tempPts(2,:)' + 1;
            colors(adM+j) = f(indC{:});    
        end
        if aInd(1) && aInd(2)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)' + 1;
            PatchPtsY(:,adM+j) = tempPts(2,:)' + 1;
            colors(adM+j) = vars.f(indC{:});    
        end
        aInd = any(tempPts>1/2,2);
        if aInd(1)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)' - 1;
            PatchPtsY(:,adM+j) = tempPts(2,:)';
            colors(adM+j) = f(indC{:});    
        end
        if aInd(2)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)';
            PatchPtsY(:,adM+j) = tempPts(2,:)' - 1;
            colors(adM+j) = f(indC{:});    
        end
        if aInd(1) && aInd(2)
            j=j+1;
            PatchPtsX(:,adM+j) = tempPts(1,:)' - 1;
            PatchPtsY(:,adM+j) = tempPts(2,:)' - 1;
            colors(adM+j) = f(indC{:});    
        end
    end
end
if vars.Boundary
    patch(PatchPtsX,PatchPtsY,colors,'Clipping','On');
else
    patch(PatchPtsX,PatchPtsY,colors,'Clipping','On','EdgeColor','None');
end
axis image
axis([-1/2 1/2 -1/2 1/2]);
colormap(vars.Colormap)
    caxis([min(f(:)),max(f(:))]);
axis off
if nargout>0
    fh = gcf;
end
end

