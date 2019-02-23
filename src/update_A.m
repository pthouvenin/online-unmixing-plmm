function A = update_A(Y,M,A,dM,muA,delta)
% Abundance update (PALM).
%%
% Code : Pierre-Antoine Thouvenin, May 16th 2015.
%%
%-------------------------------------------------------------------------%
% Input: 
% > Y       hypespectral image (L|N);
% > M       endmember matrix (L|R);
% > A       initial abundances (L|N);
% > dM      variability matrix (L|R);
% > muA     additional regularization parameter (\mu * |A|²/2 added to 
%           the "local" objective function) (1|1);
% > delta   lower bound of the Lipschitz's constant (1|1).
%
% Output:
% < A   updated abundances  (R|N).
%-------------------------------------------------------------------------%
%%
R = size(A,1);
% Lipschitz' constant computation
Md = M + dM;
mu = max([1.1*norm(Md'*Md + muA*eye(R),'fro'),delta]);
% Update
gradh = Md'*(Md*A - Y) + muA*A;
A = A - gradh/mu;
% Projection onto the simplex
A = max(bsxfun(@minus,A,max(bsxfun(@rdivide,cumsum(sort(A,1,'descend'),1)-1,(1:R)'),[],1)),0);

end

