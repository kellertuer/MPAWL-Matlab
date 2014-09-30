% Tutorial 4:
% ---
% Sampling and Interpolation
%
% Before decomposing a certain function in any TI space we have to obtain
% sampling values and perform a change of basis from the interpolating
% basis
%
% ---
% MPAWL, R. Bergmann ~ 2014-09-30
clc
format compact
setDebugLevel(3);
disp('--- Tutorial 4: Sampling and Interpolation ---');
disp(' (a) Sampling a box spline ');
disp('Let''s assume the following Xi to define a Boxspline, each column having multiplicity 1');
Xi = pi*[1,0,0.125,0,0.125; 0,1,0,0.125,-0.125] %#ok<*NOPTS>
ct = sum(Xi,2)/2; %center point
nu = ones(length(Xi),1);

n = 512;
[X,Y] = meshgrid(-pi:2*pi/n:pi-2*pi/n, -pi:2*pi/n:pi-2*pi/n);
debug('time',3,'StartTimer','Box Spline sampling on a grid');
Z = reshape(box_eval(Xi,nu,[X(:),Y(:)]+ones(length(X(:)),1)*ct'),size(X));
bsr = max(max(abs(Z)));
debug('time',3,'StopTimer','Box Spline sampling on a grid');
figure(1);
surf(X,Y,Z,'FaceColor','interp','EdgeAlpha',0.1);
axis tight
colormap rwb
caxis([-bsr,bsr]);

disp('Let''s sample this function using sample(M,fun,...) on a simple pattern: pixels');
M = 128*[1,0;0,1]
data = sample(M,@(x)(box_eval(Xi,nu,x+ones(length(x),1)*ct')),...
    'SamplingMethod','point row','File','tutorials/T4-files/sampleBoxSpline.mat');

disp('Using the Dirichlet Kernel from the last Tutorial including its Bracket Sums,');
disp('we can perform a change of basis from the interpolating basis (points from above) to the Dirichlet kernel translates');

[ckdM,dMBS] = dirichletKernel(M,'File',{'tutorials/T4-files/ckDM.mat','tutorials/T4-files/ckDM-BS.mat'});
hata = changeBasis(M,data,dMBS,'Input','time','Output','Fourier');

disp('And reconstruct its Fourier coefficients');

ckBoxSpline = coeffsSpace2Fourier(M,hata,ckdM,(size(ckdM)+1)/2);

disp('and compare the result with the above, denser sampled first data');

diffImg = real(discretePlotFourierSeries(n*[1,1],ckBoxSpline))-Z;
dIr = max(max(abs(diffImg)));
figure(2);
imagesc(diffImg,[-dIr,dIr]); 
colormap rwb
title('Illustration of the interpolation error');

disp(' (b) Sampling a linear function ');

data2 = sample(M,@(x)(sum(x)),...
    'File','tutorials/T4-files/samplelinear.mat');

hata2 = changeBasis(M,data2,dMBS,'Input','time','Output','Fourier');

disp('And reconstruct its Fourier coefficients');

ckLinear = coeffsSpace2Fourier(M,hata2,ckdM,(size(ckdM)+1)/2);

disp('and plotting the image resembles the Gibbs phenomenon at the borders');

resImg = real(discretePlotFourierSeries(n*[1,1],ckLinear));
rIr = max(max(abs(resImg)));
figure(3);
imagesc(resImg,[-rIr,rIr]); 
colormap rwb
title('Sampled version of the linear function f(x,y) = x+y');
