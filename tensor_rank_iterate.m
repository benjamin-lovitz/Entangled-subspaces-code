%%  TENSOR_RANK_ITERATE    Numerically estimates an optimization over states with tensor rank <= k
%   This function has one required argument:
%     X: a square positive semidefinite matrix
%     K: a positive integer
%     DIM: a vector with positive integer entries
%
%   val = tensor_rank_iterate(X,K,DIM) is a lower bound of the maximum value
%   of <v|X|v>, where |v> ranges over pure states with tensor rank <= K.
%   X and |v> both live in a tensor product space whose local dimensions
%   are specified by the vector DIM.
%
%   This function has one optional argument:
%     TOL (default 10^-10)
%   
%   [VAL,VP] = tensor_rank_iterate(X,K,DIM,TOL) computes the same lower
%   bound described above, stopping the algorithm once we make less
%   progress on the lower bound than TOL per iteration. The optional
%   output VP is a cell array that contains the local vectors of
%   a state |v> with tensor rank <= k that attains <v|X|v> = VAL.
%
%   requires: QETLAB (www.qetlab.com)
%             
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   last updated: September 22, 2020

function [val,vp] = tensor_rank_iterate(X,k,dim,varargin)

dX = length(X); % total dimension (product of local dimensions)
p = length(dim); % number of parties

% Set optional argument defaults: k=1, dim=sqrt(length(X)), tol=10^-5, v0
% is a random initial vector (set to -1 for now).
[tol] = opt_args({ 10^(-10) },varargin{:});

% Some of this preparation is unnecessary when k = 1, but it's cheap so we
% don't really care.
kperm = perm_inv([[1:2:2*p-1],[2:2:2*p]]);
kdim = [k*ones(1,p),dim];

psi = GHZState(k,p);
psiI = PermuteSystems(kron(psi,speye(dX)),kperm,kdim,1);
PSI = psiI*psiI';
PSIX = PermuteSystems(kron(psi*psi',X),kperm,kdim,0);

% Randomly generate a starting vector v.
for j = p:-1:1 % pre-allocate for speed
    vp{j} = randn(dim(j)*k,1) + 1i*randn(dim(j)*k,1);
    vp{j} = vp{j}/norm(vp{j});
end
v = psiI'*Tensor(vp);
v = v/norm(v);

% Preparation is done; now do the actual iteration.
it_err = 2*tol+1;
val = real(v'*X*v);
while it_err > tol
    it_err = 0;
    it_perm = randperm(p);
            
    % Loop through the p parties.
    for j = p:-1:1
        % Fix p-1 of the parties and optimize over the other party.
        for m = p:-1:1
            if(it_perm(j) == m)
                xp{m} = speye(dim(m)*k);
            else
                xp{m} = vp{m};
            end
        end
        V0 = Tensor(xp);
        V1 = V0'*PSI*V0;

        try
            [vp{it_perm(j)},nval] = eigs(V0'*PSIX*V0,V1,1,'LR');
        catch err
            % ARPACK errors? Quit out and return the best lower bound
            % found so far.
            if(strcmpi(err.identifier,'MATLAB:eigs:ARPACKroutineError'))
                return
            else
                rethrow(err);
            end
        end
        vp{it_perm(j)} = vp{it_perm(j)}/norm(vp{it_perm(j)});

        it_err = it_err + abs(val-nval);
        val = real(nval);
        v = psiI'*Tensor(vp);
        v = v/norm(v);
    end
end