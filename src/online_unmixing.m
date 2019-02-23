function [M,f,RE,aSAM_Y,time] = online_unmixing(FileName,M,dm,T,HyP,AlgoP,TypeP,AddReg,varargin)
% Online unmixing of multi-temporal hyperspectral images (MHSIs) accounting
% for temporal variability. 
%-------------------------------------------------------------------------%
%% Code
% Pierre-Antoine Thouvenin, December 5th 2015.
%-------------------------------------------------------------------------%
% Input:
% > FileName  name of the .mat file containing the multitemporal 
%             hyperspectral image (should contain a cell(1|T) named Y,
%             s.t. Y{t} = HS image at time t, saved in matlab format '-v7.3');
% > M         initial endmember matrix (L|R);
% > dm        scalar value to initialize the variability matrices 
%             (default value: 0);
% > T         number of hyperspeectral images composing the sequence;
%
% > HyP   structure containing the algorithm regularization constants:
%		>> HyP.A, HyP.M, HyP.dM : regularization constants 
%       >> HyP.sigma2   upper bound on the variability energy
%       >> HyP.kappaa2  upper bound on the energy of the average variability
%
% > AlgoP  structure containing the algorithm parameters
%       >> delta: lower bound on the Lipshitz' constants
%		>> nIter: number of iterations for the endmember update step
%		>> nEpoch: number of epochs to learn the mixture parameters
%		>> nIterPalm: number of iterations of the inner PALM algorithm
%		>> nIterDykstra: number of iterations of the inner Dykstra algorithm
%		>> rho: constant forgetting factor
%
% > TypeP structure mentionning the selected regularizations
%       >> TypeP.A: 'none' or 'tsp' (temporal smoothness prior)
%       >> TypeP.dM: 'none' or 'tsp' (temporal smoothness prior)
%       >> TypeP.M: 'none', 'mdist' (mutual distance) or 'dist' 
%       (distance to a reference endmember matrix, implying varargin = M0)
%
% > AddReg struture containing the values of optional regularization parameters
%		>> HyP.A, HyP.M, HyP.dM : additional regularization constants 
%         (optional quadratic penalizations)
%
% Output:
% < M        estimated endmember matrix (L|R)
% < f        value of the global objective function (1|nEpoch)
% < RE       quadratic reconstruction errors (nEpoch|T)
% < aSAM_Y   average angle between the reconstructed an the original data (in degrees) (nEpoch|T)
% < time     time spend to perform each epoch (1|Niter)
%-------------------------------------------------------------------------%
%%
%--------------------------------------------------------------
% Local variables
%--------------------------------------------------------------
[L,R] = size(M);

% Endmember update selection
switch TypeP.M   
    case 'none'
        lambdaM = @(C,t) max([1.1*norm(C/t + AddReg.M*speye(R),'fro'),AlgoP.delta]);
        gradM  = @(M,C,D,t) (M*C + D)/t + AddReg.M*M;
        termM = @(M) 0.5*AddReg.M*sum(M(:).^2);
        
    case 'mdist'
        G = mdist_prior(R);
        lambdaM = @(C,t) max([1.1*norm(C/t + HyP.M*G*(G') + AddReg.M*speye(R),'fro'),AlgoP.delta]);
        gradM   = @(M,C,D,t) (M*C + D)/t + AddReg.M*M + HyP.M*M*G*(G');
        termM   = @(M) 0.5*( AddReg.M*sum(M(:).^2) + HyP.M*sum(sum((M*G).^2)) ); 
        
    case 'dist'
        lambdaM = @(C,t) max([1.1*norm(C/t + (HyP.M + AddReg.M)*speye(R),'fro'),AlgoP.delta]);
        gradM   = @(M,C,D,t) (M*C + D)/t + AddReg.M*M + HyP.M*(M - varargin);
        termM = @(M) 0.5*( AddReg.M*sum(M(:).^2) + HyP.M*sum(sum((M - varargin).^2)) );
        
    otherwise
        warning('Unrecognized penalization type (TypeP.M): default option taken (''none'')');
        lambdaM = @(C,t) max([1.1*norm(C/t + AddReg.M*speye(R),'fro'),AlgoP.delta]);
        gradM  = @(M,C,D,t) (M*C + D)/t + AddReg.M*M;
        termM = @(M) 0.5*AddReg.M*sum(M(:).^2);
end

% Variable preallocation
f       = zeros(AlgoP.nEpoch,1); % objective function
RE      = zeros(AlgoP.nEpoch,T); % image by image reconstruction errors
aSAM_Y  = zeros(AlgoP.nEpoch,T);
time    = zeros(AlgoP.nEpoch,1);

C = zeros(R);
D = zeros(L,R);
E = zeros(L,R);
c = 0;

% Abundance and variability initialization
initialize_A_dM(FileName,M,T,dm);
t_tilde = 1;

% keyboard

%--------------------------------------------------------------
% Algorithm
%--------------------------------------------------------------
% Data loading
data = matfile(FileName); % FileName: to be created prior to the unmixing process
results = matfile('results.mat','Writable',true); % created in initialize_A_dM

tic

for q = 1:AlgoP.nEpoch
    
    for t = randperm(T)
        
        Y = cell2mat(data.Y(1,t));
        
        % Abundance update selection
        switch TypeP.A
            case 'none'      
                upA = @(Y,A,dM) update_A(Y,M,A,dM,AddReg.A,AlgoP.delta);
                termA = @(A) 0;
            case 'tsp'
                if t > 1
                    A_prev = cell2mat(results.A(1,t-1));
                    upA = @(Y,A,dM) update_A_tsmooth(Y,M,A,dM,AddReg.A,AlgoP.delta,HyP.A,A_prev);
                    termA = @(A) 0.5*HyP.A*sum((A(:) - A_prev(:)).^2);
                else
                    upA = @(Y,A,dM) update_A(Y,M,A,dM,AddReg.A,AlgoP.delta);
                    termA = @(A) 0;
                end
            otherwise
                warning('Unrecognized penalization : default selected (''none'')');
                upA = @(Y,A,dM) update_A(Y,M,A,dM,AddReg.A,AlgoP.delta);
        end 

        % Variability update selection
        switch TypeP.dM
            case 'none'
                up_dM = @(Y,A,dM,t) update_dM(Y,M,A,dM,AddReg.dM,AlgoP.delta,HyP.sigma2,HyP.kappa2,E,t,AlgoP.nIterDykstra);
                termdM = @(dM) 0;
            case 'tsp'
                if t > 1
                    dM_prev = cell2mat(results.dM(1,t-1));
                    up_dM = @(Y,A,dM,t) update_dM_tsmooth(Y,M,A,dM,AddReg.dM,AlgoP.delta,HyP.sigma2,HyP.kappa2,E,HyP.dM,dM_prev,t,AlgoP.nIterDykstra);
                    termdM = @(dM) 0.5*HyP.dM*sum(sum((dM - dM_prev).^2));
                else
                    up_dM = @(Y,A,dM,t) update_dM(Y,M,A,dM,AddReg.dM,AlgoP.delta,HyP.sigma2,HyP.kappa2,E,t,AlgoP.nIterDykstra);
                    termdM = @(dM) 0;
                end
            otherwise
                up_dM = @(Y,A,dM,t) update_dM(Y,M,A,dM,AddReg.dM,AlgoP.delta,HyP.sigma2,HyP.kappa2,E,t,AlgoP.nIterDykstra);
                termdM = @(dM) 0;
        end
        
        % Initialize results
        A = cell2mat(results.A(1,t));
        dM = cell2mat(results.dM(1,t));
        
        % Abundance and variability estimation (PALM)
        for k = 1:AlgoP.nIterPalm % (PALM iterations)   
		
            % Abundance update
            A = upA(Y,A,dM);
			
			% Variability update
			dM = up_dM(Y,A,dM,t_tilde);
        
        end 
        
        % Save results
        results.A(1,t) = {A};
        results.dM(1,t) = {dM};
        
        % Update auxiliary matrices
        f_factor = AlgoP.rho;
        c = f_factor*c + 0.5*sum(sum((Y - dM*A).^2));
        C = f_factor*C + A*(A');
        D = f_factor*D + (dM*A - Y)*(A');
        E = f_factor*E + dM;
        const = (1-AlgoP.rho^t_tilde)/(1-AlgoP.rho); 
        
        % Endmember update (projected gradient descent)
        M = update_M(gradM,lambdaM,M,C,D,const,AlgoP.nIter);
            
        % Error computations
        Y_hat = (M + dM)*A;
        RE(q,t) = sum((Y(:) - Y_hat(:)).^2)/numel(Y_hat);        
        aSAM_Y(q,t) = sum(abs(180*acos( sum(Y.*Y_hat,1)./(sqrt(sum(Y.^2,1)).*sqrt(sum(Y_hat.^2,1))) )/pi ))/size(Y,2);    
    
        f(q) = f(q) + 0.5*norm(Y - (M + dM)*A,'fro')^2 + termA(A) + termM(M) + termdM(dM);
        t_tilde = t_tilde + 1;
    end
    time(q) = toc;
end
