echo off
%%
%% Demo-script for the Box-Spline Evaluation (Leif Kobbelt - 11/07/96)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, home

echo on

%% This is a script file to demonstrate the use of the box_eval() function
%% which allows the stable and efficient evaluation of box-splines.

%% A box-spline is a piecewise polynomial function from R^s to R.
%% It is defined by a matrix of k mutually distinct (s-dim) column vectors.
%% The first example uses s = k = 2 ...

X  = [ 1  0 ;
       0  1 ];

%% The multiple occurrence of identical column is logged in a vector nu
%% holding the number of times each direction is repeated.

nu = [3;3];

pause  %% press any key to continue

%% The box-spline evaluation can be applied to a whole set of sample
%% points. Here we construct a uniform grid over an appropriate interval.

[xx,yy] = meshgrid(((1:20)-2)/6,((1:20)-2)/6);
 p      = [xx(:) yy(:)];

%% To also demonstrate the efficiency of the algorithm, we measure the
%% execution time.

pause  %% press any key to continue

tic
b = box_eval(X,nu,p);
toc

surf(reshape(b,20,20))

%% The surface you see right now is the graph of a biquadratic tensor
%% product Box-spline function.

pause  %% press any key to see more examples

clc, home

%% The next example has k = 4 distinct directions in R^s = R^2
%% each occurring only once.

 X      = [ 1  0   1  1 ;
            0  1  -1  1 ];
 nu     = [1;1;1;1];
[xx,yy] = meshgrid(((1:20)-2)/6,((1:20)-8)/6);
 p      = [xx(:) yy(:)];

tic
b = box_eval(X,nu,p);
toc

surf(reshape(b,20,20))

%% The resulting function is called the Zwart-Powell-element

pause  %% press any key for one more example

clc, home

%% The last example uses k = 3 directions in R^s = R^2. Each has
%% the multiplicity two.

 X      = [ 1  0  1 ;
            0  1  1 ];
 nu     = [2;2;2];
[xx,yy] = meshgrid(((1:20)-2)/5,((1:20)-2)/5);
 p      = [xx(:) yy(:)];

tic
b = box_eval(X,nu,p);
toc

surf(reshape(b,20,20))

%% This last function is obtained by repeating the directions defining
%% the Courant-element twice.

