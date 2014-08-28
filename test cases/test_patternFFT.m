setDebugLevel('text',3');
setDebugLevel('time',3');
M = [16,4; 0,16];
diag(snf(M))
b = zeros(4,64);
b(1,1) = 1.;

hatb = patternFFT(M,b)
bs = patternFFT(M,hatb)

hatb2 = patternFFT(M,b(:))
bs2 = patternFFT(M,hatb2)