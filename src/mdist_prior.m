function G = mdist_prior(R)
% Compute the matrix for the mutual distance prior.
%
% Input
% > R  number of endmembers
%
% Output
% < G  mutual distance prior matrix (\Psi = ||MG||²/2)
%-------------------------------------------------------------------------%
%%
% Code : Pierre-Antoine Thouvenin, October 29th 2015.
%%
G = zeros(R,R*R);
for l = 0:R-1
    el = zeros(R,1);
    el(l+1) = 1;
    G(:,1+l*R:(l+1)*R) = -eye(R) + el*ones(1,R);
end
G = sparse(G);
end
