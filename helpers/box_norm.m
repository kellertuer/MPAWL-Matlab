function N = box_norm(t,k,M)
%%
%%  Precomputation of normal vectors (Leif Kobbelt-  5/15/96
%%                                                  11/07/96)
%%
%%     called by box_eval()  ...  N = box_norm(t,k,M)
%%
%%  N : Returned hash table of normal vectors
%%  t : Number of rows to be selected before base case is reached
%%  k : Next row of BoxEv_X to be considered for selection
%%  M : Bitvector indicating the selected rows of BoxEv_X
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  global BoxEv_I   % Identity matrix ... eye(BoxEv_k)
  global BoxEv_s   % Dimension of domain space
  global BoxEv_X   % Matrix of distinct directions (rows) defining the box-spline

%% Allocate space according to hashing function

  N = zeros(BoxEv_s,2^k);

%% Cut leafless branches

  if (k>=t)
    if (t>0)

%% Left  half of the hash table: row X(k,:) not selected
%% Right half of the hash table: row X(k,:)     selected

      N = [box_norm(t,k-1,M),box_norm(t-1,k-1,M+BoxEv_I(k,:))];

    else

%% Normal vector is orthogonal to all selected rows ...

      N(:,1) = null(BoxEv_X(find(M),:));
    end
  end
