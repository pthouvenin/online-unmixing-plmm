function initialize_A_dM(FileName,M,T,dm)
% Initialize abundances and variability matrices.
%%
% Code : Pierre-Antoine Thouvenin, July 29th 2015.
%%
%-------------------------------------------------------------------------%
% Input: 
% > FileName    name of the .mat file containing the multitemporal hyperspectral image (contains a cell(1|T), 
% named Y, s.t. Y{t} is the HS image at time t) (file saved in matlab format '-v7.3');
% > M           endmember matrix (L|R);
% > T           number of images composing the sequence;
% > dm          scalar value to initialize the variability matrices (default value: 0).
%-------------------------------------------------------------------------%
[L,K] = size(M);
A = cell(1,T);
dM = cell(1,T);
save('results.mat','A','dM','-v7.3');             % create a file containing the cell A and dM
results = matfile('results.mat','Writable',true); % allows the cells A and dM to be progressively updated without keeping all the data into memory
data = matfile(FileName);
clear A dM;

% keyboard

for t = 1:T       
    % SUNSAL algorithm
    A1 = sunsal(M,cell2mat(data.Y(1,t)),'POSITIVITY','yes','ADDONE','yes');
    results.A(1,t) = {max(bsxfun(@minus,A1,max(bsxfun(@rdivide,cumsum(sort(A1,1,'descend'),1)-1,(1:K)'),[],1)),0)}; % projection onto the unit simplex
    results.dM(1,t) = {dm*ones(L,K)};
end

end
