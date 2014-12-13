% Tutorial 1:
% ---
% Patterns and generating Sets
%
% This tutorial introduces the basic functions for patterns and generating
% Sets and illustrates their usage.
%
% ---
% MPAWL, R. Bergmann ~ 2014-09-29
clc
format compact
start = pwd;
cd(fileparts(which(mfilename)));
run('../initMPAWL.m') %Initialize Library
setDebugLevel(3);
setDebugLevel('time',3);

disp('--- Tutorial 1: Patterns and generating Sets ---');
disp(' (a) pattern ');
disp('We start wird a simple matrix');

M = [32,4;-1,8] %#ok<*NOPTS>

disp(['and use the command pattern(M) to create the |det(M)|=',...
    num2str(abs(det(M))),' points of the pattern.']);
disp('In order to do so, we need the matrix in patternNormalForm(M):');

pnfM = patternNormalForm(M)
pts = pattern(patternNormalForm(M));

figure(1);
plot(pts(1,:),pts(2,:),'bo','MarkerSize',9','MarkerFaceColor',[.5,.5,1]);
title('Plot pattern(M) on the symmetric interval');
axis tight
axis square

disp('How many vectors (and their integer multiples) are needed for this matrix? patternDimension(M) tells us:')
dM = patternDimension(M)

disp('and which vector is that? patternBasis(M) returns them as columns of a matrix');
hM = patternBasis(M)

disp('Actually both are bases on the Smith normal form of a matrix, snf(M):');

snfM = snf(M)

disp('having any point, i.e. an integer offset of a pattern point, we can shift it back using modM(k,M) by using M=eye(2).');
disp('and the additional Option ''Target'' set to ''symmetric''');
originalp = pts(:,4)
offpt = originalp + [5,1]'
pt = modM(offpt,eye(2),'Target','symmetric')
disp('Similarly, a pattern could also be created for the unit cube, using pattert(M,''Target'',''unit'')');
disp('Of course, the basis vector above creates any pattern, but that one on the unit cube only,');
disp('when using modM(c*hM,eye(2)), c=0,...,259, is used');

disp(' (b) generating Set');
disp('In the TI spaces, the generating set of M transpose is usually used');
disp('for generatingSet(M) it is important to really use M, _not_ its patternNormalForm.');

Mt = transpose(M)
pts2 = generatingSet(Mt);
figure(2);
plot(pts2(1,:),pts2(2,:),'bo','MarkerSize',9','MarkerFaceColor',[.5,.5,1]);
title('Plot generatingSet(Mt) on the symmetric interval');
axis tight
axis square
disp('The vector reconstructing the set (with respect to modM(c*kMt,Mt)) is')
kMt = generatingSetBasis(Mt);

disp('Here, the other way around is more important in most cases:');
disp('which multiple recreates (up to modulo) a point? use generatingSetBasisDecomp(p2,Mt) with');

% for i=1:length(pts2)
p2 = [3;3] %pts2(:,i)
c = generatingSetBasisDecomp(p2,Mt)

disp('and verify by computing modM(c*kMt,Mt) for the symmetric case.');

p2r = modM(c*kMt,Mt,'Target','symmetric')

disp('All values should be integer, so modM also possesses the possibility to only produce those, by setting ''Index'' to true');
p2r2 = modM(c*kMt,Mt,'Target','symmetric','Index',true)
disp('--- End of Tutorial 1');
% end
cd(start)