%
% Tutorial 7
%   A Tree of Wavelet Decompositions
%
% ---
% MPAWL 1.0 ~ R. Bergmann 2015-01-01
clc
format compact
startFolder = pwd;
cd(fileparts(which(mfilename)));
run('../initMPAWL.m') %Initialize Library
setDebugLevel(2);
setDebugLevel('time',3);
tutFolder = 'T7-files/';
disp('--- Tutorial 7: A Tree of Wavelet Decompositions ---');
% Define a Testfunction
f = @(x) (abs(x) <= 1) .* ( (x+1).^2.*(x-1).^2 );
fp = @(x) (abs(x) <= 1) .* (4* x.*(x.^2-1) );
fpp = @(x) (abs(x) <= 1) .* (12*x.^2-4);

x = linspace(0,pi,100);
xs = 8/7/pi*x;
figure(1); plot(x,f(xs),x,fp(xs),x,fpp(xs));
title('A function f and its first and second derivative.');
xlim([0,pi]); legend('f','f''','f''''');
disp(['The given function f has a jump in its second derivative. ',char(10),...
    'it will be used as the radial part of the data fr this tutorial, i.e. ',char(10),...
    'the data has a jump on a circle of radius 7/8pi in its second order',char(10),...
    'derivatives, wth outward pointing direction.',char(10)]);

fr = @(X,Y) f(8/7/pi*abs(sqrt(X.^2+Y.^2)));

[X,Y] = meshgrid(-pi:pi/20:pi);
Z = fr(X,Y);
figure(2); surf(X,Y,Z);xlim([-pi,pi]);ylim([-pi pi]);
disp(['Analogous to Tutorial 6, we now sample and decomposie thei function,',char(10),...
    'but with different matrix decompositions arranged in a tree.',char(10)]);

M = 256*eye(2);
help decomposeData2D
data = sample(M,@(x)( fr(x(1),x(2)) ),'SamplingMethod','pointwise','File',[tutFolder,'sampledFct.mat']);

[ckMphi, ckMphiBS] = delaValleePoussinMean(1/14, M, 'File',{[tutFolder,'dlVP_ckS.mat'],[tutFolder,'dlVP_ckS_BS.mat']});

cdata = changeBasis(M,data,ckMphiBS,'Input','time','Output','Fourier');

decomposeData2D(1/14*ones(1,4),{{'Y','Ym'},{'Y','Ym'},{'X'}},M,cdata,...
    'Orthonormalize',true,'ImageOutput','Wavelet','ImagePrefix',[tutFolder,'img'],...
    'SpacePrefix',[tutFolder,'dlVP-']);

cd startFolder