function A = update_A_tsmooth(Y,M,A,dM,muA,delta,alpha,A_prev)
% Abundance update (PALM) (temporal smoothness prior).
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
% > delta   lower bound of the Lipschitz's constant (1|1);
% > alpha   regularization parameter (temporal smoothness) (1|1);
% > A_prev  abundances at time t-1 (temporal smoothness) (L|N).
%
% Output:
% < A   updated abundances  (R|N).
%-------------------------------------------------------------------------%
%%
R = size(A,1);

% Update
Md = M + dM;
mu = Md'*Md + (alpha+muA)*eye(R);
mu = max([1.1*norm(mu,'fro'),delta]);
gradh = Md'*(Md*A - Y) + muA*A + alpha*(A - A_prev);
A = A - gradh/mu;
% Projection onto the simplex
A = max(bsxfun(@minus,A,max(bsxfun(@rdivide,cumsum(sort(A,1,'descend'),1)-1,(1:R)'),[],1)),0); 

end
