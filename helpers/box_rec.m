function b = box_rec(n,m,Y,t)
%%
%%  Recursive Evaluation of Box-Splines (Leif Kobbelt -  5/15/96
%%                                                    - 11/07/96)
%%
%%       called by box_eval()  ...  b = box_rec(n,m,Y,t)
%%
%%  b : Returned vector of function values at locations BoxEv_p
%%  n : Vector holding the multiplicities of the rows in BoxEv_X
%%  m : Current position in the recursion tree (for delayed translation)
%%  Y : Matrix to compute the least norm representation of BoxEv_p
%%      with respect to the remaining directions BoxEv_X(find(n),:)
%%  t : Least norm representation of BoxEv_p
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

  if (sum(n)>BoxEv_s)

%% Recursion case ...

    b = 0;
    j = 1;

%% Sum over the remaining directions in BoxEv_X ...

    for i = 1:BoxEv_k

%% Update multiplicity of directions and position in recursion tree

      nn = n-BoxEv_I(:,i);
      mm = m+BoxEv_I(i,:);

%% Recursive calls

      if (n(i)>1)
        b = b+      t(:,j) .*box_rec(nn,m ,Y,t                         );
        b = b+(n(i)-t(:,j)).*box_rec(nn,mm,Y,t-BoxEv_J*(BoxEv_X(i,:)*Y));
        j = j+1;
      elseif (n(i)>0)

%% Update least norm representation

        Z = BoxEv_X(find(nn),:);
        if (rank(Z) == BoxEv_s)
          Z = (Z'*Z)\Z';
          b = b+      t(:,j) .*box_rec(nn,m ,Z,(BoxEv_p-BoxEv_J*(m *BoxEv_X))*Z);
          b = b+(n(i)-t(:,j)).*box_rec(nn,mm,Z,(BoxEv_p-BoxEv_J*(mm*BoxEv_X))*Z);
        end
        j = j+1;
      end
    end

%% Normalization

    b = b/(sum(n)-BoxEv_s);
  else

%% Base case ... compute characteristic function

    b = 1;

%% Delayed translations

    v = find(n);
    z = BoxEv_p-BoxEv_J*(m*BoxEv_X);

%% Check against all hyperplanes

    for i = 1:BoxEv_s

%% Lookup normal vector to current hyperplane

      NN = BoxEv_N(:,1+BoxEv_u*(n-BoxEv_I(:,v(i))));

%% Box is half-open!!!

      p  = BoxEv_X(v(i),:)*NN;
      q  = z*NN;
      b  = min(b,1-((p>0&q<0)|(p<0&q>=0)));
      q  = (BoxEv_p-BoxEv_J*((m+BoxEv_I(v(i),:))*BoxEv_X))*NN;
      b  = min(b,1-((p>0&q>=0)|(p<0&q<0)));
    end

%% Normalization

    b = b/abs(det(BoxEv_X(v(1:BoxEv_s),:)));
  end
