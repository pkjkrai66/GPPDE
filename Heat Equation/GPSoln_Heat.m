clear; close all;

addpath('./')
addpath('../')
addpath('../BIC and Likelihood Functions')
addpath('../Kernel Functions')
addpath('../Supporting Codes')

clear me mydir mydir2

% =========================================================================
% ~~~~~~~~~~~~~~~~~~~~~~ Generating Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                        ~~~~~~~~~~~~~~~
% Conditions.TimeVector = 0:.2:12;
% Conditions.TopBoundry = 50*ones(size(Conditions.TimeVector));
% Conditions.BottomBoundry = 100*ones(size(Conditions.TimeVector));
% Conditions.HeightVect = 0:1:10;
% Conditions.InitiHead = 0*ones(size(Conditions.HeightVect));
% % -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x
% 
% z.length     = 10;
% z.discreteZ  = 0:.1:z.length;
% z.deltaZ     = 0.1;
% 
% t.discreteT       = 0:.2:10;
% t.deltaT  = 0.2;
% k = 0.835;
% 
% [temperature, deltaZ] = HeatEqnImplicitMethod(Conditions, z, t, k);
% Temp = temperature(:, 2:end);
% Temp_stdev = 1.5;
% temp_Observation = Temp + normrnd(0, Temp_stdev, size(Temp));
% 
% tVector = t.discreteT(2:end);
% zVector = z.discreteZ; % Lengths at which GPs to be fitted
load('HeatSyntheticData_Errored')
TTT = length(tVector);
ZZZ = length(zVector);

[tGrid, zGrid] = meshgrid(tVector, zVector);
Ytrain = temp_Observation(:);
Xtrain = [zGrid(:), tGrid(:)];
X = Xtrain;
TZ = size(Xtrain, 1);
meany = mean(Ytrain);

figure;
H.scatterdata = scatter3(zGrid(:), tGrid(:), Ytrain);
%{
% figure(); 
% axes; 
% for i = 1:TTT
%     cla
%     plot(zVector, Temp(:,i), '-o')
%     hold on
%     plot(zVector, errored_Observation(:,i), 'xr', 'MarkerSize', 8, 'linewidth', 3)
%     plot(zVector, errored_Observation2(:,i), 'xg', 'MarkerSize', 8, 'linewidth', 2)
%     ylim([min(min(errored_Observation)), max(max(errored_Observation))]);
%     pause
% end
%}

f1 = figure();
subplot(1, 2, 1)
H.af1 = surf(zGrid, tGrid, temp_Observation);
% zstr = sprintf('temperature (t; in $\degree$C)', char(176)); zlabel(zstr);
xlabel('{\it{x}}; in (cm)'); ylabel('{\it{t}}; in (s)'); 
zlabel('temperature ({\it{T}}; in °C)')
title('temp. observation')
% -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-
%% ~~~~~~~~~~~~~~~~ Estimation of GP Hyperparameters ~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~ Gaussian Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                          ~~~~~~~~~~~~~~~~
%                  ~~  Default Matlab Function  ~~

% kparams0 = [750 5 .5];
%
% % Ini value for kernel parameters
% % SyntheticSMData_1: [750 5 .5] --> [6.5851242, 0.02544194, 0.00263082]
%
% gprMdl = fitrgp(Xtrain,Ytrain,'KernelFunction','squaredexponential',...
%                 'KernelParameters',kparams0([1,2]),'Sigma',kparams0(3));
% ypred = predict(gprMdl, Xtrain);
% estHyp = [gprMdl.KernelInformation.KernelParameters', gprMdl.Sigma];
% GPestmd_Y_MatDefault = reshape(ypred, [ZZZ, TTT]);
% figure; surf(zGrid, tGrid, GPestmd_Y_MatDefault)
% title('From Default GP function')
% estHyperpata = estHyp
% -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x

% -------------------------------------------------------------------------
%           Optimizing hyperparameters by minimizing 
%           negative-log-marginal-likelihood and it's gradient 
%           Optimization Function used: `fminunc`

kparams0 = [2 100 1];

GPoptions = optimoptions('fminunc','Algorithm','trust-region',...
                        'SpecifyObjectiveGradient',true);
fun = @(hyppar)negative_LogMarginalLikelihood(hyppar,Xtrain,Ytrain,meany);
estHyperpara = fminunc(fun, kparams0, GPoptions);
% estHyperpara = [1.88219353405592 99.9545901720395 1.55246473008603];

% -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x
% -x-x-x-x-x-x-x- OVER: Estimation of GP Hyperparameters -x-x-x-x-x-x-x

%%
tVector2  = tVector;
zVector2 = min(zVector):.5:max(zVector);
[tGrid2, zGrid2] = meshgrid(tVector2, zVector2);
X2 = [zGrid2(:), tGrid2(:)];
ZZZ2 = length(zVector2);
TTT2 = length(tVector2);
                % -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-
% -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x
clear XA TZA X2A X1A CxyA GPestmd_YA GramMatrix 
%% ~~~~~~~~~~~~~~ Implimentation of GP Hyperparameters ~~~~~~~~~~~~~~~~~~~
% 			   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gppar_temp = estHyperpara(1:2);
gpnoise_temp = estHyperpara(3);

Ctrain = RBFcov(Xtrain, Xtrain, gppar_temp);
GPConst = (Ctrain + gpnoise_temp^2*eye(TZ))\(Ytrain-meany);

Cxy = RBFcov(Xtrain, X, gppar_temp);
Cxy2 = RBFcov(Xtrain, X2, gppar_temp);
GPestmd_Y = meany + Cxy'*GPConst;
GPestmd_Y2 = meany + Cxy2'*GPConst;
GPestmd_Y_Mat = reshape(GPestmd_Y, [ZZZ, TTT]);
GPestmd_Y2_Mat = reshape(GPestmd_Y2, [ZZZ2, TTT2]);

% Calculating Performance Parameters ~~~~~~~~~~~~~~~
logMarginalLiklihood = logMarginalLikelihood(estHyperpara,Xtrain,Ytrain,meany);
R_square = CoeffOfDeter(Ytrain, GPestmd_Y); % CoeffOfDeter(Observed, Predicted)
RMSE = RootMeanSqrErr(Ytrain, GPestmd_Y);   % RootMeanSqrErr(Observed, Predicted)

% Plotting ~~~~~~~~~~~~~~~
figure(f1);
subplot(1, 2, 2)
H.p1a1f = surf(zGrid, tGrid, GPestmd_Y_Mat, 'linestyle', 'none');
zlabel('temperature ({\it{T}}; in °C)')
xlabel('{\it{x}}; in (cm)'); ylabel('{\it{t}}; in (s)');
title('GP estimated observation')

% ---------- Predictive Covariance --------------------
PredCovDiag = diag(diag(Ctrain)) + gpnoise_temp^2*eye(TZ);
PredCov = PredCovDiag - Cxy'*((Ctrain + gpnoise_temp^2*eye(TZ))\Cxy);
PredVarianceVect = diag(PredCov);
VarianceMatrix = reshape(PredVarianceVect, [ZZZ, TTT]);

% plotting Variance
figure; imagesc(tVector, zVector, VarianceMatrix)
ylabel('{\it{x}} (cm)')
xlabel('{\it{t}} (s)')
myVline([.25*max(zVector) .75*max(zVector)], 'color', 'black' )
myHline([.25*max(tVector) .75*max(tVector)], 'color', 'black' )
% ======== Calculating Derivatives of Theta (Gaussian Process) ============
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DiffWRT_Dim_z = 1;          % The derivative is to be carried out 
                        % w.r.t to 1st variable z
% -------------- First Derivative wrt z -----------------------

dC_z = dRBFcov(Xtrain, X, gppar_temp, DiffWRT_Dim_z);
delY_delz = dC_z'*GPConst;
delY_delz_GPMat = reshape(delY_delz, [ZZZ, TTT]);

dC_z2 = dRBFcov(Xtrain, X2, gppar_temp, DiffWRT_Dim_z);
delY_delz2 = dC_z2'*GPConst;
delY_delz_GPMat2 = reshape(delY_delz2, [ZZZ2, TTT2]);

%   CovMat_delY_delz{ci} = feval(@dxdzGPcov, X, X, gppar, DiffWRT_Dim_z)...
%                             -dC_z*((C + noise_param^2*eye(TZ))\dC_z');
%   VarianceMatrix_delY_delz{ci} = reshape(diag(CovMat_delSM_delz{ci}), ...
%                                             [ZZZ, TTT]);

% ------------ Second Derivative wrt z -------------------------

ddC_z = dxdxRBFcov(Xtrain, X, gppar_temp, DiffWRT_Dim_z);
del2Y_delz2 = ddC_z'*GPConst;
del2Y_delz2_GPMat = reshape(del2Y_delz2, [ZZZ, TTT]);

ddC_z2 = dxdxRBFcov(Xtrain, X2,gppar_temp, DiffWRT_Dim_z);
del2Y_delz22 = ddC_z2'*GPConst;
del2Y_delz2_GPMat2 = reshape(del2Y_delz22, [ZZZ2, TTT2]);

%   CovMat_del2Y_delz2{ci} = feval(@d4C_dx2dz2, X, X, gppar, DiffWRT_Dim_z,1)...
%                             -ddC_z*((C + noise_param^2*eye(TZ))\ddC_z');
%   VatianceMatrix_del2Y_delz2{ci} = reshape(diag(CovMat_del2SM_delz2{ci}), ...
%                                             [ZZZ, TTT]);
% ----------- Derivative of Theta wrt time (Gaussian Process) =============

DiffWRT_Dim_t = 2;         % The derivative is to be carried out
                       % w.r.t to 2nd variable t
dC_t = dRBFcov(Xtrain, X, gppar_temp, DiffWRT_Dim_t);
delY_deltime = dC_t'*GPConst;
delY_delt_GPMat = reshape(delY_deltime, [ZZZ, TTT]);

dC_t2 = dRBFcov(Xtrain, X2, gppar_temp, DiffWRT_Dim_t);
delY_delt2 = dC_t2'*GPConst;
delY_delt_GPMat2 = reshape(delY_delt2, [ZZZ2, TTT2]);

%   CovMat_delY_delt{ci} = feval(@dxdzGPcov, X, X, gppar, DiffWRT_Dim_t)...
%                             -dC_t*((C + noise_param^2*eye(TZ))\dC_t');
%   VarianceMatrix_delY_delt{ci} = reshape(diag(CovMat_delSM_deltime{ci}), ...
%                                             [ZZZ, TTT]);

clearvars PredCovDiag PredCov PredVarianceVect DiffWRT_Dim_t DiffWRT_Dim_z
clearvars Cxy2 GPConst GPestmd_Y GPestmd_Y2 gppar noise_param 
clearvars CovMat_delY_delz CovMat_delY_deltime head2 tempI deltaT dC_t dC_t2
clearvars delSM_deltime delY_delt2 ddC_z ddC_z2 del2Y_delz2 del2Y_delz22
clearvars dC_z dC_z2 delY_delz delY_delz2 meany Cxy

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                    Estimation of PDE Parameters          
%                    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DataTrim.z = 0.25; % fraction_of_Zdata_removed {Value in fraction in one Side}
DataTrim.t = 0.25; % fraction_of_tdata_removed {Value in fraction in one Side}
pdefun = @(PDEParameter)errorFunctOfHeatEqn(PDEParameter, ...
                      delY_delt_GPMat, del2Y_delz2_GPMat, VarianceMatrix,DataTrim);

k_Guess = 5;
PDEoptions = optimoptions('fminunc','Algorithm','quasi-newton');
estPDparam = fminunc(pdefun, k_Guess, PDEoptions);

% Contour Plotting

eqnParam01 = 0.1:0.075:2.5;
errNorm1 = NaN(length(eqnParam01), 1);
errNorm2 = NaN(length(eqnParam01), 1);
errNormInf = NaN(length(eqnParam01), 1);

for i = 1:length(eqnParam01)
    errNorm2(i) = feval(pdefun, eqnParam01(i));
    errNorm22(i) = errorFunctOfHeatEqn(eqnParam01(i), ...
                      delY_delt_GPMat, del2Y_delz2_GPMat, VarianceMatrix,DataTrim);
end

errors = cell(1, 3);
errors(:) = {NaN};
errors(1, [1 2 3]) = {errNorm2, errNorm1, errNormInf};
errtitle = {'errNorm2', 'errNorm1', 'errNormInf'};

    figure
    plot(eqnParam01, errors{1})
    [~] = myVline(k);
    [~] = myVline(estPDparam, 'color', 'red', 'LineStyle', '--');
    xlabel('\lambda (cm^2s^{-1})')
    ylabel('sqError (°C^{2}s^{-2})')

% -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x

%% ~~~~~~~~~~~~~~~~~~~~~~~~~ Plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                  ~~~~~~~~

f2 = figure();     % Derivative WRT dist
for i = 1:length(tVector)
    
  H.fTa1 = subplot(1,3,1);  % ------- theta -------------
    cla; hold on
    H.fTa1p1 = plot(zVector, Temp(:, i), 'r-');
    H.fTa1p2 = plot(zVector, temp_Observation(:,i), 'g.');
    H.fTa1p3 = plot(zVector2, GPestmd_Y2_Mat(:, i), 'b-*');
    legend([H.fTa1p1 H.fTa1p2 H.fTa1p3], {'noiseless data', 'Observations', 'Gaussian fit'})
    xlabel('{\it{x}}; in (cm)'); 
    ylabel('{\it{T}} (s)'); 
    
  H.fTa2 = subplot(1,3,2);  % ------ First Derivative wrt x -----------
    cla; hold on
    H.fTa2p2 = plot(zVector2, delY_delz_GPMat2(:, i), 'b-*');
    title('First Derivative wrt x')
    xlabel('{\it{x}} (cm)'); 
    ylabel('$\frac{\partial^2 T}{\partial z^2}$','Interpreter','latex', 'FontSize',20);

  H.fTa4 = subplot(1,3,3);  % ------ Second Derivative wrt x ------------
    cla; hold on;
    H.fTa4p2 = plot(zVector2, del2Y_delz2_GPMat2(:, i), 'b-*');
    xlabel('{\it{x}} (cm)'); 
    title('Second Derivative wrt x')
    ylabel('$\frac{\partial^2 T}{\partial x^2}$','Interpreter','latex', 'FontSize',20);

  supText = sprintf('at Time: %.2f', tVector(i));
  supHandel = suplabel(supText, 't');
  pause

end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~