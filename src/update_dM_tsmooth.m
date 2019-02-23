function dM = update_dM_tsmooth(Y,M,A,dM,mu,delta,sigma2,kappa2,E,gamma,dM_prev,t,nIterDykstra)
% Variability update (PALM) (temporal smoothness prior).
%%
% Code : Pierre-Antoine Thouvenin, May 7th 2015.
%%
%-------------------------------------------------------------------------%
% Input: 
% > Y             hypespectral image (L|N);
% > M             endmember matrix (L|R);
% > A             abundances (L|N);
% > dM            initial variability terms (L|R);
% > mu            additional regularization parameter (\mu * |dM|²/2 added to 
%                 the "local" objective function) (1|1);
% > delta         lower bound for the Lipschitz's coefficients;
% > sigma2        upper bound on the variability energy 
%                 ( |dM|² \leq \sigma_2) (1|1);
% > kappa2        upper bound on the energy of the average variability 
%                 ( |E[dM]|² \leq \kappa_2) (1|1);
% > E        	  auxiliary variable (described in the main algorithm) (L|R)
% > gamma         regularization parameter (temporal smoothness) (1|1);
% > dM_prev       variability term at time t-1 (temporal smoothness) (L|R);
% > t		      index of the current iteration (1|1);
% > nIterDykstra  number of iterations for the Dykstra algorithm (1|1).
%
% Output:
% < dM    updated variability matrix (L|R).
%-------------------------------------------------------------------------%
%%
[L,R] = size(dM);
Proj1 = @(y) proj_fb(y,zeros(L,R),sqrt(sigma2));
Proj2 = @(y) proj_fb(y,-E,t*sqrt(kappa2));

% Update
nu = A*(A') + (gamma+mu)*eye(R);
nu = max([1.1*sqrt(sum(nu(:).^2)),delta]);
gradh = ((M + dM)*A - Y)*(A') + mu*dM + gamma*(dM - dM_prev);
U = dM - gradh/nu;
% Approximate projection onto { dM | |dMt|² <= sigma² and |\sum_t dMt|² <= kappa²}
dM = dykstra(U,Proj1,Proj2,nIterDykstra);

end
