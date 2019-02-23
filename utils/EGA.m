function R = EGA(Y)
% Estimation of the endmember number R using an eigen-gap approach.
%-------------------------------------------------------------------------%
% Input:
% > Y     hyperspectral image (L|N).
%
% Output:
% < R     estimated endmember number.
%-------------------------------------------------------------------------%
%% Reference
% [1] A. Halimi, P. Honeine, M. Kharouf, C. Richard and J.-Y. Tourneret,
% "Estimating the Intrinsic Dimension of Hyperspectral Images Using an 
% Eigen-Gap Approach," submitted. Preprint available on arXiv.
%-------------------------------------------------------------------------%
%% Code
% Pierre-Antoine Thouvenin, December 5th 2015.
%-------------------------------------------------------------------------%

%% Remarks
% - problème liés au mauvais conditionnement des données pour l'inversion !
% - code très long à exécuter en raison de la présence des deux boucles for !
% - il y a certainement un problème dans l'implémentation : résultats très
% peu convaincants jusqu'à présent

% to be optimized with ARMADILLO
[L,N] = size(Y);

% Constants
c = L/N;
beta_c = (1 + sqrt(c))*((1 + sqrt(1/c))^(1/3));
psi_N = 4*sqrt(2*log(log(N)));
dN = psi_N*beta_c/(N^(2/3));

% Data covariance matrix
Ry = Y*(Y')/N;

% Ry eigenvectors (decreasing order)
[V,lambdaRy] = eig(Ry,'vector'); % vp déjà classées, il suffit d'inverser l'ordre
[lambdaRy,id] = sort(lambdaRy,'descend');
V = V(:,id);

% Noise covariance matrix estimate
E = zeros(N,L);

for l = 1:L
    yl = Y;
    yl(l,:) = [];
    E(:,l) = Y(l,:)' - (yl')*((yl*(yl'))\(yl*Y(l,:)')); % problème inversion !! 
    % trouver un autre moyen dans le cas où l'inversion est impossible
end
Sig = E'*E/N;

% Rs eigenvectors (decreasing order)
Rs = Ry - Sig;
[W,lambdaRs] = eig(Rs,'vector');
[~,id] = sort(lambdaRs,'descend');
W = W(:,id);

% Computation of \sigma_k^2
% sigma2 = zeros(L,1);
% for l = 1:L
%     sigma2(l) = (V(:,l)'*Sig*W(:,l))/(V(:,l)'*W(:,l));    
% end
sigma2 = (sum(V.*(Sig*W),2)./sum(V.*W,2));     

% Estimated endmember number
Delta = -diff(lambdaRy./sigma2); % re-vérifier !
K = find(Delta < dN,1,'first');
R = K + 1;

end
