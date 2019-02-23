function [U,P,y_bar,Y_proj] = PCA(Y,R)
% Principal component analysis (adapted from find_endm)
%%
% Code : Pierre-Antoine Thouvenin, March 20th 2015.
%%
%-------------------------------------------------------------------------%
% Input: 
% > Y     hypespectral image (L|N);
% > R     endmember number.
%
% Output:
% < U      inverse projector (L|R-1)
% < P      projector (R-1|L)
% < y_bar  data mean (L|1);
% < Y_proj projected data (PCA space) (R-1|N).
%-------------------------------------------------------------------------%
%%
%--------------------------------------------------------------
% PCA
%--------------------------------------------------------------
y_bar = mean(Y,2);
Rmat = bsxfun(@minus,Y,y_bar);
Rmat = Rmat*Rmat';   % empirical covariance matrix

OPTIONS.disp = 0;    % diagnostic information display level
OPTIONS.maxit = 600; % maximum number of iterations
[V,D] = eigs(Rmat,R-1,'LM',OPTIONS) ; % first R-1 eigenvectors

% projector
P = V';
% P = D^(-1/2)*(V'); % permet de remédier aux problèmes numériques pouvant apparaître si le simplex est trop "étiré" 
                     % dans une direction (différence d'amplitude importante entre les valeurs propres)
% inverse projector
U = pinv(P);  % pinv(P)*D^(1/2);

% projecting
Y_proj = P*bsxfun(@minus,Y,y_bar);  % projection of the data on the K principal axes.

end

