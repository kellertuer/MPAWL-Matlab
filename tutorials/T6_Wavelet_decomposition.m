% Tutorial 6:
% ---
% Wavelet decomposition
%
% We illustrate the decomposition of a sampled box spline into wavelet
% parts
%
% ---
% MPAWL, R. Bergmann ~ 2014-10-02
clc
format compact
start = pwd;
cd(fileparts(which(mfilename)));
run('../initMPAWL.m') %Initialize Library
setDebugLevel(2);
setDebugLevel('time',3);
disp('--- Tutorial 6: Wavelet decomposition ---');
disp(' (a) Sampling a box spline [see Tutorial 4]');
Xi = pi*[1,0,0.125,0,0.125; 0,1,0,0.125,-0.125] %#ok<*NOPTS>
ct = sum(Xi,2)/2; %center point
nu = ones(length(Xi),1);

n = 1024;
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
title('Box Spline we will sample');
%%
M = n/2*[1,0;0,1]
data = sample(M,@(x)(box_eval(Xi,nu,x+ones(length(x),1)*ct')),...
    'SamplingMethod','point row','File','T6-files/sampleBoxSpline.mat');

[ckdM,dMBS] = dirichletKernel(M,'File',{'T6-files/ckDM.mat','T6-files/ckDM-BS.mat'});
origin = (size(ckdM)+1)/2;

cdata = changeBasis(M,data,dMBS,'Output','Fourier');

ckdata = coeffsSpace2Fourier(M,cdata,ckdM,origin);

figure(2);
imgSamples = real(FourierSeries2Img([n,n],ckdata) );
rIr = max(max(abs(imgSamples)));
imagesc(imgSamples,[-rIr,rIr]); 
colormap rwb
title('Sampled version of the Box spline');
%%
matStr = 'X';
[CoeffsDXs,CoeffsDXw] = dirichletKernelSubspaces(M,dilationMatrix2D(matStr),...
    'File',{'T6-files/coeffsDXS.mat','T6-files/coeffsDXW.mat'});

ckDXs = coeffsSpace2Fourier(M,CoeffsDXs,ckdM,origin);
ckDXw = coeffsSpace2Fourier(M,CoeffsDXw,ckdM,origin);

[cdataXS,cdataXW] = patternFWT(M,dilationMatrix2D(matStr),cdata,CoeffsDXs,CoeffsDXw);
%%
ckDataXS = coeffsSpace2Fourier(round(inv(dilationMatrix2D(matStr))*M),cdataXS,ckDXs,origin);
figure(3);
imgSamplesS = real( FourierSeries2Img([n,n],ckDataXS) );
rIr = max(max(abs(imgSamplesS)));
imagesc(imgSamplesS,[-rIr,rIr]); 
colormap rwb
title('Scaling space part');
axis tight
axis square

ckDataXW = coeffsSpace2Fourier(round(inv(dilationMatrix2D(matStr))*M),cdataXW,ckDXw,origin);
figure(4);
imgSamplesW = real( FourierSeries2Img([n,n],ckDataXW) );
rIr = max(max(abs(imgSamplesW)));
imagesc(imgSamplesW,[-rIr,rIr]); 
colormap rwb
title('Wavelet space part');
axis tight
axis square
%%
disp(' (b) De la Vallée Poussin means and a decomposition Tree');
disp(' The localization of the edges is not that good using the dirichlet kernel. Let''s try de la Vallée Poussin');

t = 1/1024;

[ckdlVPM,dlVPMBS] = delaValleePoussinMean(t,M,'File',{'T6-files/ckdlVPM.mat','T6-files/ckdlVPM-BS.mat'});
origin2 = (size(ckdlVPM)+1)/2;
cdata2 = changeBasis(M,data,dlVPMBS,'Output','Fourier');

disp('But we will use the multilevel approach this time, even only for one level: decomposeData2D');
decomposeData2D(t,{{'X'}},M,cdata2,'ImageOutput','Both')
%%
disp('But it is also possible to perform more than one decomposition and to save/load coefficients and images');
decomposeData2D(t*ones(1,6),{{'D'},{'X'},{'X'},{'X'},{'Y'}},M,cdata2,'ImageOutput','Both','ImagePrefix','T6-files/img','SpacePrefix','T6-files/space')