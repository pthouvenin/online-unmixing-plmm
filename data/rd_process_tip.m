%% Create .mat files (matrix cube) for the image time series
clc; clear all; close all
format compact;
%--------------------------------------------------------------------------------------------%
%% RAPPEL
% - lexicographical ordering: [H,W,L] -> [L,H*W] Y = (reshape(permute(data,[2 1 3]),H*W,L))'; 
%                             [L,H*W] -> [H,W,L] data = permute(reshape(Y',W,H,L),[2 1 3]);
%
% - MATLAB ordering (column-wise): [H,W,L] -> [L,H*W] Y = reshape(data,H*W,L)';
%                                  [L,H*W] -> [H,W,L] data = reshape(Y',H,W,L);
%--------------------------------------------------------------------------------------------%
%%  
% Parameter initialization
T = 6;
R = 3;
H = 150;
W = 110;
N = H*W;

mask = [3:102,115:150,185:221]; % spectral mask [series 2 box 9 (178)] 
L = numel(mask);

indices = outlier_patch(63:66,77:80,H,W,0);
id = true(1,N);
id(indices) = false;

% Extract time series parameters (from the first image)
FileName = strcat('Im1','.hdr');
[W, H, L_true, interleave, offset, byte_order, data_type, wavelength, wavelength_unit] = extract_parameters(cd,FileName); 
data_type = strcat(data_type,'=>','double'); 
N = H*W;
Y = cell(1,T);

for t = 1:T
	FileName = strcat('Im',num2str(t));
	data = multibandread(fullfile(cd,FileName),[H,W,L_true],data_type,offset,interleave,byte_order,{'Column',1:W},{'Row',1:H},{'Band',1:L_true}); %(PathName,FileName)
    Y{t} = reshape(data(:,:,mask),N,L)'/10000; % column-major ordering (not lexicographical ordering)
    if t == 5 % image with outliers
        y = reshape(data(:,:,mask),N,L)'/10000;
        id = setdiff(1:N,indices);
        Y{t} = y(:,id);
    end
end

save('rd_tip.mat','Y','H','W','L','wavelength','wavelength_unit','mask','indices','-v7.3');
disp('... DONE');
disp('---------------------------------------------------------------------------');