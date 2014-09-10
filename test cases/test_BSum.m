% Transcribed from Example 3 of the Mathematica implementation
%%
M = [32,4; -1,8];
abs(det(M))
IndMax = getMaxIndex(transpose(M))

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

dMBS