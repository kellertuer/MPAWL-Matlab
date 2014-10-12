setDebugLevel(3);
M = [2,0;0,2];
hata = [0.001, 0.002; 0.003, 0.005];

ckphi = [0.01, 0.02, 0.03; 0.1, 0.2, 0.3; 1 2 3];

origin = [2,2];

ck1 = coeffsSpace2Fourier_old(M,hata,ckphi,origin)

ck2 = coeffsSpace2Fourier(M,hata,ckphi,origin)