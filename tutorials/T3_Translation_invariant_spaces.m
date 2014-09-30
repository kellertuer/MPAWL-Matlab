% Tutorial 3:
% ---
% Translation invariant spaces
%
% This tutorial demonstrates how to construct and work with functions,
% whose translates generate an m=|det(M)| dimensional space.
%
% ---
% MPAWL, R. Bergmann ~ 2014-09-30
clc
format compact
setDebugLevel(3);
disp('--- Tutorial 3: Translation invariant spaces ---');
disp('We again use the matrix');
M = [32,4;-1,8] %#ok<*NOPTS>

disp('Having a pattern of cardinality');
adM = abs(det(M))

disp('[Repetition from Tutorial 1] We take a closer look at generatingSet(transpose(M))');
Mt = transpose(M)
pts2 = generatingSet(Mt);
figure(1);
plot(pts2(1,:),pts2(2,:),'bo','MarkerSize',9','MarkerFaceColor',[.5,.5,1]);
title('Plot generatingSet(Mt) on the symmetric interval');
axis tight
axis square

disp('there is a maximal index, that characterizes the minimal surrounding rectange of G(M^T), using getMaxIndex');
IndMax = getMaxIndex(transpose(M))

disp('used by the corresponding Dirichlet kernel, which is defined in frequency:');
ckDM = dirichletKernel(M,'Orthonormalize',false)

disp('But its translates may also be orthonormalized');
ckDMorth = dirichletKernel(M)
ckBSq = bracketSums(ckDMorth,(size(ckDMorth)+1)/2,M,'Validate',false,'Compute','absolute Squares')

disp('Or even used as an interpolating function by emploing the usual Bracket sums');
ckBSIP = bracketSums(ckDM,(size(ckDMorth)+1)/2,M,'Validate',false,'Compute','Bracket')
ckDMIP = coeffsSpace2Fourier(M,1./(abs(det(M))*ckBSIP),ckDM,(size(ckDMorth)+1)/2)

disp('The Bracket sums of this last kernel are');
ckBSDMIP = bracketSums(ckDMIP,(size(ckDMorth)+1)/2,M,'Validate',false,'Compute','Bracket')

disp('Which translates as time coefficients to (again transposed for space reasons');
patternFFT(M,ckBSDMIP)'
