%%
%-------------------------------------------------------------------------%
%          Online Unmixing of Multitemporal Hyperspectral                 %
%           Images Accounting for Spectral Variability                    %
%-------------------------------------------------------------------------%
%% File data
% File : main.m
% Author : P.A. Thouvenin (18/02/2015)
% Last modified : 28/06/2016
clc, clear all, close all, format compact;
addpath src; addpath utils; addpath data;
%-------------------------------------------------------------------------%
%% REF.
% P.-A. Thouvenin, N. Dobigeon and J.-Y. Tourneret,
% "Online unmixing of multitemporal hyperspectral images accounting for 
% spectral variability", IEEE Trans. Image Processing, vol. 25, no. 9, 
% pp. 3979-3990 Sep. 2016.
%-------------------------------------------------------------------------%
%% REMARKS
% v2 : corrected version, slightly optimized instructions.
%-------------------------------------------------------------------------%
%%
%--------------------------------------------------------------
% General parameters
%--------------------------------------------------------------
% Endmember number
R = 3;
% Number of images
T = 5; %10
% FileName (name of the .mat data file, saved with the '-v7.3' option to allow partial loading of matrices)
FileName = 'rd_tip.mat';
v = 'rd';
% Endmember initialization
% M0 = initialize_M(FileName,R,T); % use VCA on a single image if needed.
% save('data/M1','M0')
% return;

%--------------------------------------------------------------
% Parameters
%--------------------------------------------------------------
% Initial endmember matrix
load('M1','M0'); 
% Algorithm parameters
AlgoP = struct('delta',1e-10,'nIter',50,'nEpoch',20,'nIterPalm',50,'nIterDykstra',50,'rho',0.98);
% Regularization type
TypeP = struct('A','none','M','mdist','dM','tsp'); % abundance temporal smoothness removed, since outliers have been removed from the image at t = 5 
AddReg = struct('A',0,'M',0,'dM',0);
% Regularization constants
sigma2 = 1;
kappa2 = 1e-1;
HyP = struct('A',0,'M',1e-4,'dM',0,'sigma2',sigma2,'kappa2',kappa2); % real data
dm = 0;

%--------------------------------------------------------------
% Algorithm
%--------------------------------------------------------------
disp('Online unmixing: in progress...')
[M,f,RE,aSAM_Y,time] = online_unmixing(FileName,M0,dm,T,HyP,AlgoP,TypeP,AddReg); % M0
disp('-> Operation complete.')
disp('---------------------------------------------------------------------------');

%--------------------------------------------------------------
% Save results 
%--------------------------------------------------------------
save(['unmixing_results_',v],'M','f','HyP','AlgoP','AddReg','TypeP','time','M0','RE','aSAM_Y');

% % Slight modification to account for the outlier patch removed in the 5th
% % image
results = matfile('results.mat','Writable',true);
load('rd_tip.mat','indices');
N = 150*110;
A1 = zeros(R,N);
id = true(1,N);
id(indices) = false;
A1(:,id) = cell2mat(results.A(1,5));
results.A(1,5) = {A1};

%--------------------------------------------------------------
% Plot endmembers / corrupted endmembers
%--------------------------------------------------------------
results = matfile('results.mat');
dM = cell2mat(results.dM);
for r = 1:R
    figure;
    plot(M(:,r),'r-','Linewidth',1.5);
    for t = 1:T
        hold on;
        plot(M(:,r) + dM(:,r + (t-1)*R),'b-.')
    end
end
