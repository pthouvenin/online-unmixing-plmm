function M = initialize_M(FileName,R,T)
%%
% Code : Pierre-Antoine Thouvenin, July 29th 2015.
%%
%-------------------------------------------------------------------------%
% Input: 
% > FileName    name of the .mat file containing the multitemporal hyperspectral image (contains a cell(1|T));
% > R           number of endmembers;
% > T           number of images composing the sequence;
%-------------------------------------------------------------------------%

Y_hull = [];
data = matfile(FileName);

for t = 1:T
    Y = cell2mat(data.Y(1,t));
    
    % PCA of the current image
    [U,V,Y_bar,Y_proj] = PCA(Y,R);
    
    % Convex hull (current image)
    D1 = Y_proj';
    S = convhulln(D1);        % convhull : operates row by row, only for 2D and 3D
    Y2 = Y(:,unique(S(:,1))); % only the points from the convex envelope of each image are kept in memory  
    
    % Updated convex hull
    Y_hull = [Y_hull,Y2];
end
M = vca(Y_hull,'Endmembers',R,'verbose','off');