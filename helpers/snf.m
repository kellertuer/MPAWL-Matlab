% [U,S,V] = snf(A)
% S = snf(A)
%
% Computes the Smith normal form S of a given integer matrix A and their
% corresponding transformation matrices, where S = UAV holds.
%
% INPUT
%   A : An integer matrix of dimension m times n
%
% OUTPUT
%   U : An m times m transformation matrix 
%   S : The Smith normal form Diagonal matrix consisting of the elementary
%       divisors
%   V : An n times n transformation matrix 
%
% If only one or no return value is specified, the matrix S is returned.
% 
% MPAWL 1.0, D. Merkert ~ 2014-07-31 ~ last edit: 2014-08-19 (R. Bergmann)

% TODO: The code could be improved by (1) adding more comments and (2)
% reduce the possibilities of too large values in between (maybe by
% adapting the code from A. Pascoletti.

function [U,S,V] = snf(A)

assert(sum(sum(abs(A))) ~= 0, 'MPAWL:snf',...
    'The input matrix must not be the zero matrix');

m = size(A,1);
n = size(A,2);

S = A;
U = eye(m,m);
V = eye(n,n);
index = 1;

while (index < m && index < n)   
    % Obtain the minimal value (t) and one set of indices with t=S(u,v)
    [t,u,v] = matrixAbsMin(S,index);
    d = gcdm(S,index);
    while (d ~= t)
        if (S(u,v) < 0)
            [U,S,V] = multRow(S,u,-1,U,V);
        end
        if (existsNonDivisableEntryInColOrRow(S,index,u,v) == 1)
            [U,S,V] = reduceGCDCase1(S,index,u,v,U,V);
        else
            [U,S,V] = reduceToZero(S,index,u,v,U,V);
            [t,u,v] = matrixAbsMin(index,S);
            d = gcdm(S,index);
            
            if (d ~= t)
                [U,S,V] = reduceGCDCase2(S,index,u,v,U,V);
                [U,S,V] = reduceGCDCase1(S,index,u,v,U,V);
            end
        end
        [t,u,v] = matrixAbsMin(S,index);
        d = gcdm(S,index);
    end % end while d ~= t    
    [U,S,V] = putMinimalElementToCorrectPosition(S,index,u,v,U,V);
    if (S(index,index) < 0)
        [U,S,V] = multRow(S,index,-1,U,V);
    end
    [U,S,V] = reduceToZero(S,index,index,index,U,V);
    index = index+1;
end % end while indices
% Debug check
% verifySmithNormalForm(U,S,V,A);
if nargout < 2 % return just S into the one variable or to ans.
    U = S;
end
end

% Local helpers to simplify the above code

function res = gcdm(M,index)
% checkForNaN(M);
res = M(index,index);
for i = index:size(M,1)
    for j = index:size(M,2)
        res = gcd(res,M(i,j));
    end
end
end

function res = existsNonDivisableEntryInColOrRow(S,index,u,v)
% checkForNaN(S);
res = 0;
for j=index:size(S,2)
    if (mod(S(u,j),S(u,v)) ~= 0)
        res = 1;
        return;
    end
end
for i=index:size(S,1)
    if (mod(S(i,v),S(u,v)) ~= 0)
        res = 1;
        return;
    end
end
end

function [U,S,V] = reduceGCDCase1(S,index,u,v,U,V)
% checkForNaN(S);
auv = S(u,v);
for j=index:size(S,2)
    alm = S(u,j);
    if (mod(alm,auv) ~= 0)
%        r = rem(alm,auv);
        q = fix(alm./auv);
        [U,S,V] = subCol(S,j,v,q,U,V);
    end
end
for i=index:size(S,1)
    alm = S(i,v);
    if (mod(alm,auv) ~= 0)
%        r = rem(alm,auv);
        q = fix(alm./auv);
        [U,S,V] = subRow(S,i,u,q,U,V);
    end
end
end

function [U,S,V] = reduceGCDCase2(S,index,u,v,U,V)
% checkForNaN(S);
auv = S(u,v);
for i = index:size(S,1)
    for j = index:size(S,2)
        aij = S(i,j);
        if (mod(aij,auv) ~= 0)
            [U,S,V] = subRow(S,u,i,-1,U,V);
            return;
        end
    end
end
end

function [U,S,V] = reduceToZero(S,index,u,v,U,V)
% checkForNaN(S);
auv = S(u,v);
for j=(index+1):size(S,2)
    alm = S(u,j);
    z = alm/auv;
    [U,S,V] = subCol(S,j,v,z,U,V);
end
for i=(index+1):size(S,1)
    alm = S(i,v);
    z = alm/auv;
    [U,S,V] = subRow(S,i,u,z,U,V);
end
end

function [U,S,V] = putMinimalElementToCorrectPosition(S,index,u,v,U,V)
[U,S,V] = swapRow(S,index,u,U,V);
[U,S,V] = swapCol(S,index,v,U,V);
end

function [t,u,v] = matrixAbsMin(M,index)
t = inf;
u = NaN;
v = NaN;

for i=index:size(M,1)
    for j=index:size(M,2)
        if (M(i,j)~=0)
            if (abs(M(i,j)) < t)
                t = abs(M(i,j));
                u = i;
                v = j;
            end
        end
    end
end
end

function [U,S,V] = swapRow(S,i,j,U,V)
tmp = S(i,:);
S(i,:) = S(j,:);
S(j,:) = tmp;

tmp = U(i,:);
U(i,:) = U(j,:);
U(j,:) = tmp;
end

function [U,S,V] = swapCol(S,i,j,U,V)
tmp = S(:,i);
S(:,i) = S(:,j);
S(:,j) = tmp;

tmp = V(:,i);
V(:,i) = V(:,j);
V(:,j) = tmp;
end

function [U,S,V] = multRow(S,i,f,U,V)
S(i,:) = f*S(i,:);
U(i,:) = f*U(i,:);
end

function [U,S,V] = multCol(S,j,f,U,V)
S(:,j) = f*S(:,j);
V(:,j) = f*V(:,j);
end

function [U,S,V] = subRow(S,i,j,f,U,V)
S(i,:) = S(i,:) - S(j,:) .* f;
U(i,:) = U(i,:) - U(j,:) .* f;
end

function [U,S,V] = subCol(S,i,j,f,U,V)
S(:,i) = S(:,i) - S(:,j) .* f;
V(:,i) = V(:,i) - V(:,j) .* f;
end

% function checkForNaN(S)
% if (sum(sum(isnan(S))) ~= 0)
%     error('You have a NaN');
% end
% if (sum(sum(isinf(S))) ~= 0)
%     error('You have a Inf');
% end
% if (norm(S-round(S)) > 0)
%     error('You have non-integers');
% end
% end

% function verifySmithNormalForm(U,S,V,A)
% checkForNaN(U);
% checkForNaN(S);
% checkForNaN(V);
% checkForNaN(A);
% if (sum(sum(U*A*V-S)) ~= 0)
%     error('Matrices do not fulfill U*A*V-S=0');
% end
% old = 1;
% for i=1:(min(size(S,1),size(S,2)))
%     if (mod(S(i,i),old) ~= 0)
%         error('Divisability property not fulfilled');
%     end
%     old = S(i,i);
% end
% end
