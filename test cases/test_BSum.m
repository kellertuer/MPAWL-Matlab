% Transcribed from Example 3 of the Mathematica implementation
%%
M = [32,4; -1,8];
assert(abs(det(M))==260,'The determinant is not correct, did you change the matrix?')
patternBasis(M);
IndMax = getMaxIndex(transpose(M));
origin = IndMax + 1;

ckDM = zeros((2*IndMax+1));
summ = nestedFor(-IndMax,IndMax);
while summ.hasNext() % Loop for the Dirichlet-Case
    v = summ.next();
    t = max(abs(transpose(inv(M))*v'));
    tInd = num2cell(v+origin);
        % M(sub2ind(size(M),num2cell(t){:}))
        % this is even worse than Mathematica handling indices!
        % It even needs to store this cell! WTF
    if (t<=1/2)
        ckDM(sub2ind(size(ckDM),tInd{:})) = 1;
    else
        ckDM(sub2ind(size(ckDM),tInd{:})) = 0;
    end
end

ckDM

dMBS = bracketSums(ckDM,origin,M);

dMBS;
assert(dMBS(131)==2,'This value should be 2');

ckDMIP = coeffsSpace2Fourier(1./(abs(det(M))*dMBS),ckDM,origin,M)

t = bracketSums(ckDMIP,origin,M);
assert(all(t == ones(size(t))*1/260),...
    'for the IP they all should be 1/260 = 0.0038');
