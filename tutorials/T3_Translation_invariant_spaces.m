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
start = pwd;
cd(fileparts(which(mfilename)));
run('../initMPAWL.m') %Initialize Library
setDebugLevel(3);
setDebugLevel('time',3);
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
ckDMIP = coeffsSpace2Fourier(M,1./(ckBSIP),ckDM,(size(ckDMorth)+1)/2)

disp('The Bracket sums of this last kernel are');
ckBSDMIP = bracketSums(ckDMIP,(size(ckDMorth)+1)/2,M,'Validate',false,'Compute','Bracket')

disp('Which translates as time coefficients to (again transposed for space reasons');
patternIFFT(M,ckBSDMIP)'


disp('his can also be visualized, but we have to take the real part only, the imaginary part is of order 10^(-14) only though)');
figure(2);
imgDM = 1/abs(det(M))*real(FourierSeries2Img(256*[1,1],ckDM));
%correct axis
range = max(max(abs(imgDM)));
imagesc(imgDM,[-range,range]); %unify range of the image in order to have white=middle=zero
colormap rwb
title('Plot of the modified dirichlet Kernel DM(x,y), which is already interpolating');

%%
figure(3);
[X,Y] = meshgrid(-pi:2*pi/256:pi-2*pi/256,-pi:2*pi/256:pi-2*pi/256);
[Xl,Yl] = meshgrid(-pi:2*pi/64:pi-2*pi/64,-pi:2*pi/64:pi-2*pi/64);
imgDMl = 1/abs(det(M))*real(discretePlotFourierSeries(64*[1,1],ckDM));
rangel = max(max(abs(imgDMl)));
hold on
surf(X,Y,imgDM,'FaceAlpha',.8,'FaceColor','interp','EdgeColor','k','EdgeAlpha',0.2);
colormap rwb
caxis([-rangel,rangel]);
% mesh(Xl,Yl,imgDMl+0.01,'EdgeColor','k','FaceAlpha',0); hidden=off
hold off
%%
figure(4);
imgDMo = 1/abs(det(M))*real(FourierSeries2Img(256*[1,1],ckDMorth));
%correct axis
range = max(max(abs(imgDMo)));
imagesc(imgDMo,[-range,range]); %unify range of the image in order to have white=middle=zero
colormap rwb
title('Plot of the modified orthonormalized dirichlet Kernel DM(x,y)');

figure(5);
imgDiff = 1/abs(det(M))*real(FourierSeries2Img(256*[1,1],ckDMorth-ckDM));
%correct axis
range = max(max(abs(imgDiff)));
imagesc(imgDiff,[-range,range]); %unify range of the image in order to have white=middle=zero
colormap rwb
title('Difference of these two functions.');
