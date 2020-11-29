% The code loads Soil Moisture data set for GPPDE Paper
% This script loads the following cases:
%      SyntheticSMData_1ErrorFree               Case 1
%      SyntheticSMData_1ErrorFree_12Points      Case 2
%      SyntheticSMData_1Errored                 Case 3
%      SyntheticSMData_1Errored_12Points        Case 4

clear;  %close all;
addpath('./')
addpath('./SM_Functions/')
addpath('../')
addpath('../BIC and Likelihood Functions')
addpath('../Kernel Functions')
addpath('../Supporting Codes')

% =========================================================================
% ~~~~~~~~~~~~~~~~~~~~~~ Generating Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                        ~~~~~~~~~~~~~~~
%     alpha = 0.05;   % 1/cm
%     N = 2;          % Dimention Less
%     Ks = 1.5;       % cm/h
%     
%     % Conditions1 = importConditionsfromXLS('Conditions');
%     Conditions = importConditionsfromTable('Conditions_SynData1.mat');
%     % Conditions for Synthetic Data stored in following files:
      %     SyntheticData1 = Conditions_SynData1.mat
      %     SyntheticData2{from Celia 1990} = ConditionsData2_Celia1990.mat

%     Z = 60;
%     if Z > max(Conditions.HeightVect)
%         msgbox('Depth ''Z'' can not be more than maxmum of Depth Vector in initial condition')
%     end
%     
%     z.ProfileDepth = 60;    % This may be equal or less then max(Conditions.HeightVect)
%     z.discreteZ = 0:1:60;     % Discretized values of depth(z) [cm]
%     z.deltaZ = 0.5;           
%     
%     t.discreteT = 0:0.5:30;
%     t.deltaT = 0.5;
%     
%     thta.s = 0.368;
%     thta.r = 0.102;
%     Param.alpha = alpha;
%     Param.N = N;
%     M = 1-1/N;
%     
%     head = Richards_Celia_MixForm_VanGenuchten(Conditions, z, t, thta, Param, Ks);
%     % Coverting head to Soil moisture using VanGenuchten
%     theta = SoilMoistureVanGenuchten(head, thta, Param);
% 
% 
%     % clearvars head me mydir mydir2
%     % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% 
%     SM = theta(:, 2:end);
%     tVector = t.discreteT(2:end);
%     zVector = z.discreteZ; % Depths at which GPs to be fitted
%     SM_Stdev = 0.000;
%     SM_Observation = SM + normrnd(0, SM_Stdev, size(SM));
%     Ytrain = SM_Observation(:);    % Training set

DataSetName = 'SyntheticSMData_1ErrorFree';
[zVector, tVector, SM_Observation, thta, Ks, alpha, n, hyper0 ] = ...
                                      load_PaperSyntheticData(DataSetName);

%%%%%%%%%%%%%%%%%%%%%%%%%% TRAINING SET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ytrain = SM_Observation(:)% + normrnd(0, .0040);		% Training set

TTT = length(tVector);
ZZZ = length(zVector);

[tGrid, zGrid] = meshgrid(tVector, zVector);
Xtrain = [zGrid(:), tGrid(:)];      % Training Set
% -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x

%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING SET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tVectorTest = tVector;
zVectorTest = zVector;
[tGrid_test, zGrid_test] = meshgrid(tVectorTest, zVectorTest);

X = [zGrid_test(:), tGrid_test(:)];      % Testing Set
ZZZ_Test = length(zVectorTest);
TTT_Test = length(tVectorTest);
% -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x

TZ = size(Xtrain, 1);
meany = mean(Ytrain);
% figure();
% axes;
% for i = 1:TTT
%     cla
%     plot(zVector, theta(:,i), '-o')
%     hold on
%     plot(zVector, SM_Observation(:,i), 'xr', 'MarkerSize', 8, 'linewidth', 3)
% %     plot(zVector, errored_Observation2(:,i), 'xg', 'MarkerSize', 8, 'linewidth', 2)
%     ylim([min(min(SM_Observation)), max(max(SM_Observation))]);
%     pause
% end

% figure()      % ===================
% H.p1a1f = surf(zGrid, tGrid, SM);
% title('SM(no error)');  xlabel('distance'); ylabel('time'); zlabel('Soil Moisture')

f1 = figure();
H.af1 = surf(zGrid, tGrid, SM_Observation(:, :));
xlabel('height({\it{z}}; in cm)'); ylabel('time({\it{t}}; in h)');
zlabel('soil moisture(\theta)');
title('Observed SM'); 

figure; plot3(Xtrain(:,1), Xtrain(:,2), Ytrain, '.')
figure; mesh(zGrid, tGrid, SM_Observation(:, :));

hndlf = figure();
H.a1f = subplot(1,2,1);      % ===================
H.p1a1f = surf(zGrid, tGrid, SM_Observation);
xlabel('height ({\it{z}}; in cm)'); ylabel('time ({\it{t}}; in h)');
zlabel('soil moisture (\theta; m^3m^{-3})');
title('Observed SM'); 


clearvars Conditions t 
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

kparams0 = hyper0 ;
GPoptions = optimoptions('fminunc','Algorithm','trust-region',...
                        'SpecifyObjectiveGradient',true);
fun = @(hyppar)negative_LogMarginalLikelihood(hyppar,Xtrain,Ytrain,meany);
estHyperpara = fminunc(fun, kparams0, GPoptions);
% estHyperpara = [3.42831211976612, 0.624365255593804, 0.0000128384776762633];
% -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x
% -x-x-x-x-x-x-x- OVER: Estimation of GP Hyperparameters -x-x-x-x-x-x-x

%% ~~~~~~~~~~~~~~ Implimentation of GP Hyperparameters ~~~~~~~~~~~~~~~~~~~
% 			   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gppar_temp = estHyperpara(1:2);
gpnoise_temp = estHyperpara(3);

Ctrain = RBFcov(Xtrain, Xtrain, gppar_temp);
GPConst = (Ctrain + gpnoise_temp^2*eye(TZ))\(Ytrain-meany);

Cxy = RBFcov(Xtrain, X, gppar_temp);
GPestmd_Y = meany + Cxy'*GPConst;
GPestmd_Y_Mat = reshape(GPestmd_Y, [ZZZ, TTT]);

logMarginalLiklihood = logMarginalLikelihood(estHyperpara,Xtrain,Ytrain,meany);
R_square = CoeffOfDeter(Ytrain, GPestmd_Y); % CoeffOfDeter(Observed, Predicted)
RMSE = RootMeanSqrErr(Ytrain, GPestmd_Y);   % RootMeanSqrErr(Observed, Predicted)

figure(hndlf);
H.a1f = subplot(1,2,2);
H.p1a1f = surf(zGrid, tGrid, GPestmd_Y_Mat);
xlabel('height ({\it{z}}; in cm)'); ylabel('time ({\it{t}}; in h)');
zlabel('soil moisture (\theta; m^3m^{-3})');
title('GP Estimated SM');

% ---------- Predictive Covariance --------------------
PredCovDiag = diag(diag(Ctrain)) + gpnoise_temp^2*eye(TZ);
PredCov = PredCovDiag - Cxy'*((Ctrain + gpnoise_temp^2*eye(TZ))\Cxy);
PredVarianceVect = diag(PredCov);
VarianceMatrix = reshape(PredVarianceVect, [ZZZ, TTT]);

% ======== Calculating Derivatives of Theta (Gaussian Process) ============
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DiffWRT_Dim_z = 1;          % The derivative is to be carried out 
                            % w.r.t to 1st variable z
% -------------- First Derivative wrt z -----------------------

  dC_z = dRBFcov(Xtrain, X, gppar_temp, DiffWRT_Dim_z);
  delSM_delz = dC_z'*GPConst;
  delSM_delz_GPMat = reshape(delSM_delz, [ZZZ, TTT]);

% ------------ Second Derivative wrt z -------------------------
  ddC_z = dxdxRBFcov(Xtrain, X, gppar_temp, DiffWRT_Dim_z);
  del2SM_delz2 = ddC_z'*GPConst;
  del2SM_delz2_GPMat = reshape(del2SM_delz2, [ZZZ, TTT]);

% ----------- Derivative of Theta wrt time (Gaussian Process) =============

DiffWRT_Dim_t = 2;         % The derivative is to be carried out
                           % w.r.t to 2nd variable t
  dC_t = dRBFcov(Xtrain, X, gppar_temp, DiffWRT_Dim_t);
  delSM_deltime = dC_t'*GPConst;
  delSM_deltime_GPMat = reshape(delSM_deltime, [ZZZ, TTT]);

clearvars PredCovDiag PredCov PredVarianceVect DiffWRT_Dim_t DiffWRT_Dim_z
clearvars zGrid  GPConst GPestmd_Y gppar noise_param delSM_delz ddC_z Cxy
clearvars CovMat_delSM_delz CovMat_delSM_deltime tempI deltaT dC_z meany   
% -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                    Estimation of PDE Parameters          
%                    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DataTrim.z = 0.25; % fraction_of_Zdata_removed {Value in fraction in one Side}
DataTrim.t = 0.25; % fraction_of_tdata_removed {Value in fraction in one Side}
pdefun = @(PDEParam_N)CalErrRichardGP4(PDEParam_N, ...
            delSM_delz_GPMat, del2SM_delz2_GPMat, delSM_deltime_GPMat, ...
            GPestmd_Y_Mat, Ks, thta, VarianceMatrix, DataTrim);

PDEparams0 = [1.45, 2.1];      % [alpha, N]
PDEoptions = optimoptions('fminunc','Algorithm','quasi-newton');
estPDEparams = fminunc(pdefun, PDEparams0, PDEoptions);

alpha0 =  0.03:.005:0.06;
N0 =       1.75:.015:2.25;

errMat      = NaN(length(alpha0), length(N0));
errNorm1    = NaN(length(alpha0), length(N0));
errNorm2    = NaN(length(alpha0), length(N0));
errNormInf  = NaN(length(alpha0), length(N0));

for i = 1:length(alpha0)
 PDEParam(1) = alpha0(i);
  for j = 1:length(N0)
    PDEParam(2)    = N0(j);
    errNorm2(i, j) = feval(pdefun, PDEParam);
  end
end

[X2N0, X1Alp] = meshgrid(N0, alpha0);

figure;
contour(X1Alp, X2N0, errNorm2, 'ShowText', 'On')
myVline(alpha, 'LineStyle','-'); myHline(n,'LineStyle', '-');
xlabel('\alpha (cm^{-1})',  'FontSize', 13);
ylabel('\it n',  'FontSize', 13);
hold on
plot(estPDEparams(1), estPDEparams(2), 'Marker','+', 'MarkerSize', 9, 'LineWidth', 2, 'Color', 'red')
plot(alpha, n, 'Marker','x', 'MarkerSize', 8, 'LineWidth', 2, 'Color', 'Black')

%clear X2N0 X1Alp X1Alp X2N0 errNorm1 errNorm2 errNormInf alpha0 N0

%% ~~~~~~~~~~~~~~~~~~~~~~~~~ Plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                  ~~~~~~~~
zdataselection = zVector([ceil(DataTrim.z*ZZZ)+1, floor((1-DataTrim.z)*ZZZ)-1]);
tdataselection = tVector([ceil(DataTrim.t*TTT)+1, floor((1-DataTrim.t)*TTT)-1]);

f2 = figure();     % Derivative WRT dist
for i = 1:length(tVector)
  H.fTa1 = subplot(1,3,1);  % ------- theta -------------
  hold on
%     H.fTa1p1 = plot(zVector, SM(:,i), 'g*');
    H.fTa1p2 = plot(zVector, SM_Observation(:,i), 'ko');
    H.fTa1p3 = plot(zVector, GPestmd_Y_Mat(:, i), 'b--');
    [~] = myVline(zdataselection, 'LineStyle', '--');
    box on 
    legend([H.fTa1p2 H.fTa1p3],...
          {'SM Observation', 'GP fit to SM'}, 'FontSize',10)
    xlabel('{\it{z}} (cm)',  'FontSize', 14);
    ylabel('\theta', 'FontSize',15)
    title('(a)', 'FontSize', 12)
    
  H.fTa2 = subplot(1,3,2);  % ------ First Derivative wrt x -----------
    H.fTa2p2 = plot(zVector, delSM_delz_GPMat(:, i), 'b-');
    [~] = myVline(zdataselection, 'LineStyle', '--');
    box on
    title('(b)', 'FontSize', 12)
    xlabel('{\it{z}} (cm)',  'FontSize', 14);
    lt2 = ylabel('$\frac{\partial\theta}{\partial z}$');
    set(lt2,'Interpreter','latex', 'FontSize',20);
    
  H.fTa4 = subplot(1,3,3);  % ------ Second Derivative wrt x ------------
    H.fTa4p2 = plot(zVector, del2SM_delz2_GPMat(:, i), 'b-');
    [~] = myVline(zdataselection, 'LineStyle', '--');
    box on
    xlabel('{\it{z}} (cm)',  'FontSize', 14);
    lt3 = ylabel({'$\frac{\partial^2\theta}{\partial z^2}$'});
    set(lt3,'Interpreter','latex', 'FontSize',20)
    title('(c)', 'FontSize', 12)

  supText = sprintf('at Time: %.1f h', tVector(i));
  supHandel = suplabel(supText, 't');
  
  pause
    %   flpath = 'D:\Academics\Matlab Codes\Richards Solution\SavedFigs';
    %   figstr = sprintf('RichFig%d', i);
    %   saveas(gcf, fullfile(flpath, figstr), 'png')
end


H.fT = figure();     % Derivative WRT time ~~~~~~~~~~~~~~~~~~
for j = 1:length(zVector)
  H.fT2a1 = subplot(1,2,1);  % ------- theta -------------
    cla; hold on
%     H.fTa1p1 = plot(zVector, SM(j,:), 'g*');
    H.fT2a1p2 = plot(tVector, SM_Observation(j,:), 'ko');
    H.fT2a1p3 = plot(tVector, GPestmd_Y_Mat(j,:), 'b--');
    [~] = myVline(tdataselection, 'LineStyle', '--');
    box on 
    legend([H.fT2a1p2 H.fT2a1p3],...
        {'\theta Observation', 'Gaussian fit to \theta'}, 'FontSize',10)
    xlabel('{\it{t}} (h)',  'FontSize', 13); 
    ylabel('\theta', 'FontSize', 15)
    
  H.fT2a2 = subplot(1,2,2);  % ------ First Derivative wrt x -----------
    H.fT2a2p2 = plot(tVector, delSM_deltime_GPMat(j,:), 'b-');
    box on
    xlabel('{\it{t}} (h)',  'FontSize', 13);
    ybl2 = ylabel('$\frac{\partial\theta}{\partial t}$');
    myVline(tdataselection, 'LineStyle', '--')
    set(ybl2,'Interpreter','latex', 'FontSize',20);
%     [~] = myVline(tdataselection, 'LineStyle', '--');

  supText = sprintf('at height: %.1f cm', zVector(j));
  supHandel = suplabel(supText, 't');
  pause
    %   flpath = 'D:\Academics\Matlab Codes\Richards Solution\SavedFigs';
    %   figstr = sprintf('RichFig%d', i);
    %   saveas(gcf, fullfile(flpath, figstr), 'png')
end
clearvars supText supHandel

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~