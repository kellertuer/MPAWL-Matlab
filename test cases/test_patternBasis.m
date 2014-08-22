%
% Verify Validity of Pattern Basis
%
% ---
% MPAWL 1.0, R. Bermnann ~ 2014-08-20

disp('The first matrix has pattern dimension two and elementary divisors 4 and 64');
M1 = [16,4;0,16]
patternDimension(M1) % 2
diag(snf(M1)) % 4,64
patternBasis(M1) % [0, 1/64; 1/4, -1/16]
generatingSetBasis(transpose(M1)) %[4,1;1,16]