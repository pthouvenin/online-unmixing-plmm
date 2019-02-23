function M = update_M(gradM,lambdaM,M,C,D,t,N_iter)
% Generic endmember update (PALM, which reduces to projected gradient steps 
% in the present case).
%%
% Code : Pierre-Antoine Thouvenin, November 17th 2015.
%%
%-------------------------------------------------------------------------%
% Input: 
% > gradM      gradient (w.r.t. M) of the global objective function  (L|R);
% > lambdaM    anonymous function to compute the Lipschitz constant of gradM (1|1);
% > M          initial endmember matrix (L|R);
% > C,D        auxiliary variables (described in the main algorithm) (L|R);
% > t          index of the current iteration (1|1);
% > N_iter     number of iterations ;
%
% Output:
% < M   updated endmember matrix (L|R).
%-------------------------------------------------------------------------%
%%
%--------------------------------------------------------------
% Update (projected gradient descent)
%--------------------------------------------------------------
lambda = lambdaM(C,t);
for n = 1:N_iter
    M = M - gradM(M,C,D,t)/lambda;
    % Projections (M >= 0)
    M = (M>=0).*M;
end

end
