% Tutorial 2:
% ---
% The Fourier transform on arbitrary patterns
%
% This tutorial introduces the Fourier-Transform on arbitrary patterns
%
% ---
% MPAWL, R. Bergmann ~ 2014-09-29
clc
format compact
clc
format compact
start = pwd;
cd(fileparts(which(mfilename)));
run('../initMPAWL.m') %Initialize Library
setDebugLevel(3);
setDebugLevel('time',3);
disp('--- Tutorial 2: The Fourier Transform on arbitrary patterns');
disp(' Imagine the pattern given by');
M = [16,4;0,16] %#ok<*NOPTS>

disp('[Repetition from Tutorial 1] Its pattern dimension is ');
dM = patternDimension(M)

disp('its basis');
hM = patternBasis(M)

disp('And its cycle length (elementary divisors)');
elDivs = diag(snf(M))

disp('and its pattern by');
pts = pattern(patternNormalForm(patternNormalForm(M)));
figure(1);
plot(pts(1,:),pts(2,:),'bo','MarkerSize',9','MarkerFaceColor',[.5,.5,1]);
title('Plot pattern(M) on the symmetric interval');
axis tight
axis square

disp('[/Repetition from Tutorial 1]');
disp('Having a matrix of sampling values of size of the elDivs can be seen as');
disp('sampling a function on the pattern points, here the delta peak.');
b = zeros(elDivs');
b(1,1) = 1

disp('The Fourier transform now needs both M and b: patternFFT(M,b):');
hatb = patternFFT(M,b)

disp('it also works, if you provide a vector. It is interpreted as b(:),')
disp('i.e. reshaped to the matrix before performing the (here 2D) fft and then reshaped to a vector again for the result.');
disp('For display reasons, the transpose of the result is shown');
hatbvec = patternFFT(M,b(:))'


disp('where the Fourier transform is not scaled (neither by 1/abs(det(M)) not by sqareroot of ~).');
disp('where here both circles are ordered with respect to the generatingSetBasis of Mt, i.e. as multiples of');
hMt = generatingSetBasis(transpose(M))

disp('And in fact, the inverse patternIFFT(M,b) reproduces');
b2 = patternIFFT(M,hatb)

disp('And of course similarily for the vector (again shown the transpose of the result for space reasons).');
b2vec = patternIFFT(M,hatbvec')'

disp('Which of course reshaped also resembles b, where reshape is done with the elementary divisors as dimensions.');
b3 = reshape(b2vec',elDivs')