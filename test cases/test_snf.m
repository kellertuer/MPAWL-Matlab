%
%
% Test cases for the Smith normal form algorithm
%
% MPAWL 1.0, R. Bergmann ~ 2014-08-19
disp('Starting test of the smith normal form function snf(A).');
disp('M1 - stays as it is:');
M1 = [5,0; 0,5]
[U1,S1,V1] = snf(M1)

assert(all(all(M1==S1)),'MPAWL:test_snf','The first test failed, M1 should be S1');

try 
disp('M2 - The zero matrix should produce an error.');
M2 = [0,0,0;0,0,0;0,0,0]
goterror = false;
snf(M2)
catch err
    disp('As Expected, there was an error for the zero matrix');
    disp(getReport(err));
    goterror=true;
end
assert(goterror,'MPAWL:test_snf','The zero matrix should produce an error');
% Has elementary divisors 2 and 64.
M3 = [16,0;14,8];
snf(M3)
assert(ans(1,1)==2 && ans(2,2)==64,'MPAWL:test_snf','The elementary divisors are not correct');