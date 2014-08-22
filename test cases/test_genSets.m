%
% Verify Validity of generating Set functions
%
% ---
% MPAWL 1.0, R. Bermnann ~ 2014-08-22

disp('The first matrix has pattern dimension two and elementary divsors 1 and 260');
M1 = [32,4;-1,8];
dM1 = patternDimension(M1)
assert(dM1==1,'The pattern dimension should be 1');
elem = diag(snf(M1))
assert(all(elem==[1,260]'),'The elementary divisors are not [1,260] as expected');
V = generatingSetBasis(transpose(M1))
assert(all(V==[0;1]),'The basis should consist of only one element namely [0,1]''');

v = [3,0]';
k = generatingSetBasisDecomp([3,0],transpose(M1));
v2 = modM(V*k,transpose(M1))
assert(all(v==v2),'The decomposition is wrong, it does not repoduce the original vector');
