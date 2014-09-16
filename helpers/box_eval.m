function b = box_eval(X,nu,p)
%%
%%  Efficient and stable evaluation of box-splines (Leif Kobbelt -  5/15/96
%%                                                               - 11/07/96)
%%
%%        [Vectorized version: Spline is evaluated
%%                    at all locations in p simultaneously]
%%
%%                      b = box_eval(X,nu,p)
%%
%%  b  : Returned vector of function values at locations p
%%  X  : Matrix of distinct directions defining the box-spline
%%       (columns = s-dim direction vectors)
%%  nu : Vector holding the multiplicities of the columns in X
%%  p  : Vector of points where the spline is to be evaluated
%%       (rows = s-dim point vectors)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  global BoxEv_I   % Identity matrix ... eye(BoxEv_k)
  global BoxEv_J   % ones(length(BoxEv_p),1) --> matrix - BoxEv_J * row_vector
  global BoxEv_k   % Number of rows in BoxEv_X
  global BoxEv_N   % Hash table of normal vectors
  global BoxEv_s   % Dimension of domain space
  global BoxEv_u   % Hashing function
  global BoxEv_p   % Vector of points where the spline is to be evaluated
  global BoxEv_X   % Matrix of distinct directions (rows) defining the box-spline

  BoxEv_p = p;
  BoxEv_X = X';

  [BoxEv_k,BoxEv_s] = size(BoxEv_X);
  [      n,BoxEv_s] = size(p);

  BoxEv_I = eye(BoxEv_k);
  BoxEv_J = ones(n,1);

%% Hashing function

  BoxEv_u = 2.^(0:BoxEv_k-1);

%% Solve normal equation for least norm representation of p

  Y = BoxEv_X(find(nu),:);
  Y = (Y'*Y)\Y';

%% Compute hash table of normal vectors

  BoxEv_N = box_norm(BoxEv_s-1,BoxEv_k,zeros(1,BoxEv_k));

%% Do recursion ...

  b = box_rec(nu,zeros(1,BoxEv_k),Y,p*Y);

%% Garbage Collection

clear global BoxEv_I BoxEv_J BoxEv_k BoxEv_N BoxEv_s BoxEv_u BoxEv_p BoxEv_X
