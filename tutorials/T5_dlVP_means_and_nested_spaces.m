% Tutorial 5:
% ---
% de la Vallée POussin means and nested spaces
%
% What the Dirichlet kernels were used to in the last Tutorial, can also be
% done with the de la Vallée Poussin means, and furthermore, this tutorial
% will introduce the construction of nested spaces
%
% ---
% MPAWL, R. Bergmann ~ 2014-09-30
clc
format compact
setDebugLevel(3);
disp('--- Tutorial 5: de la Vallée POussin means and nested spaces ---');
disp('While the dirichlet kernel is just a modified characteristic function,');
disp('the de la Vallée Poussin kernel is based on a function g, for example the pyramid function');

figure(1);
plot(-1:1/16:1,pyramidFunction(1/8,-1:1/16:1));
title('A 1D pyramid function with alpha = 1/8');

disp('Or 2D even with different steepness');
[X,Y] = meshgrid(-1:1/32:1,-1:1/32:1);
Z = reshape(  pyramidFunction([1/8,1/4],[X(:)';Y(:)']), size(X));
figure(2);
surf(X,Y,Z)
title('a 2D pyramid function with alpha = 1/8[1,2]');

disp('This can be used for a de la Vallée Poussin mean, e.g. for the simple pattern');
M = 32*[1,0;0,1] %#ok<*NOPTS>

ckphiM = delaValleePoussinMean(@(x)(pyramidFunction([1/14,1/14],x)),M,'Support',0.5+1/14);
% Short hand : delaValleePoussinMean([1/14,1/14],M);
max1 = (size(ckphiM')+1)/2-1;
[Xc,Yc] = meshgrid(-max1(1):max1(1),-max1(2):max1(2));
figure(3);
scatter3(Xc(:),Yc(:),ckphiM(:),'bo')
axis tight
title('Fourier coefficients of a de la Vallée Poussin mean');

disp('Or using the shorthand with a vector g and a nonsquared matrix N (or even a number that gets vectorized.');
N = [16,12;0,16]
ckphiN = delaValleePoussinMean(1/14,N);
max2 = (size(ckphiN')+1)/2-1;
[Xc2,Yc2] = meshgrid(-max2(1):max2(1),-max2(2):max2(2));
figure(4);
scatter3(Xc2(:),Yc2(:),ckphiN(:),'bo')
title('special matrix N');
disp('corresponding function');

%%
figure(5);
n=128;
imgphiN = 1/abs(det(N))*real(discretePlotFourierSeries(n*[1,1],ckphiN));
%correct axis
rangeN = max(max(abs(imgphiN)));
[Xn,Yn] = meshgrid(-pi:2*pi/n:(pi-2*pi/n),-pi:2*pi/n:(pi-2*pi/n));
surf(Xn,Yn,imgphiN,'FaceAlpha',.8,'FaceColor','interp','EdgeColor','k','EdgeAlpha',0.2);
colormap rwb
caxis([-rangeN,rangeN]);
title('de la Vallee Poussin kernel');

ckDN = dirichletKernel(N);
max3 = (size(ckDN')+1)/2-1;
figure(6);
n=128;
imgDN = 1/abs(det(N))*real(discretePlotFourierSeries(n*[1,1],ckDN));
%correct axis
rangeDN = max(max(abs(imgDN)));
[Xn,Yn] = meshgrid(-pi:2*pi/n:(pi-2*pi/n),-pi:2*pi/n:(pi-2*pi/n));
surf(Xn,Yn,imgDN,'FaceAlpha',.8,'FaceColor','interp','EdgeColor','k','EdgeAlpha',0.2);
colormap rwb
caxis([-rangeN,rangeN]);
title('Corresponding dirichlet Kernel');

%%
disp(' (b) Subspaces');
disp(' Lets divide the Figure 5 mean into two spaces with respect to X');

dilX = dilationMatrix2D('X')
%%
disp('Both coefficient sets are with respect to the ckphiN Fourier coefficients above. Origin, i.e. c_0 is always the center for these.');
[coeffS, coeffW] = delaValleePoussinSubspaces(1/14,N,dilX);

disp('To reconstruct, use coeffsSpace2Fourier(N,coeffs[S/W],ckphiN,origin)');

figure(7);
ckphiSubS = coeffsSpace2Fourier(N,coeffS,ckphiN,(size(ckphiN)+1)/2);
imgScale = 2/abs(det(N))*real(discretePlotFourierSeries(n*[1,1],ckphiSubS));
%correct axis
rangeNS = max(max(abs(imgScale)));
surf(Xn,Yn,imgScale,'FaceAlpha',.8,'FaceColor','interp','EdgeColor','k','EdgeAlpha',0.2);
axis tight
colormap rwb
caxis([-rangeNS,rangeNS]);
title('Scaling subspace function');
 
figure(8);
ckphiSubW = coeffsSpace2Fourier(N,coeffW,ckphiN,(size(ckphiN)+1)/2);
imgWave = 2/abs(det(N))*real(discretePlotFourierSeries(n*[1,1],ckphiSubW));
%correct axis
rangeNW = max(max(abs(imgWave)));
surf(Xn,Yn,imgWave,'FaceAlpha',.8,'FaceColor','interp','EdgeColor','k','EdgeAlpha',0.2);
axis tight
colormap rwb
caxis([-rangeNW,rangeNW]);
title('Wavelet subspace function');