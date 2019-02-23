function [M2,id] = reorder(M1,M2)
% Reorder the endmembers contained in M2.
%%
% Code : Pierre-Antoine Thouvenin, February 16th 2015.
%%
%-------------------------------------------------------------------------%
% Inputs:
% > M1     endmembers 1 (reference order)
% > M2     endmembers 2
% Outputs:
% < M2        reordered endmembers.
%-------------------------------------------------------------------------%
%%
K = size(M1,2);

%--------------------------------------------------------------
% Spectral angle and l2-norm computation for data reordering
%-------------------------------------------------------------- 
s = zeros(K,1);
r = zeros(K,1);
SAM = zeros(K,1);
GMSE_Mk = zeros(K,1);
id1 = zeros(K,1);
id2 = zeros(K,1);

for k = 1:K
    for l = 1:K
        s(l) = 180*acos( (M1(:,k).')*M2(:,l) /(norm(M1(:,k))*norm(M2(:,l))) )/pi;
        r(l) = norm(M1(:,k) - M2(:,l),2);
    end
    [SAM(k), id1(k)] = min(s);
    [GMSE_Mk(k),id2(k)] = min(r);
end

% keyboard

id1 = unique(id1,'stable');
id2 = unique(id2,'stable');
id = id1;
if length(id1) < K
    if length(id2) < K
        id = [id1;setdiff(1:K,id1)'];
    else
        id = id2;
    end
end
M2 = M2(:,id);