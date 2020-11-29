function errNorm2 = CalErrRichardGP4...
      (PDEparams, deltheta_delz, del2theta_delz2, deltheta_deltime, Y,Ks, thta,...
      VarianceMatrix, DataTrim)

  
% [errNorm1, errNorm2, errNormInf] = CalErrRichardGP3...
%       (PDEparams, deltheta_delz, del2theta_delz2, deltheta_deltime, Y, thta, Ks,...
%       VarianceMatrix)

VarianceMatrix = ones(size(VarianceMatrix));

alpha  = PDEparams(1);
N = PDEparams(2);
M = 1-1/N;
theta = Y;

% ------------------ RHS ------------------

delK_delTheta = delK_delTheta_fofTheta(theta, M, thta, Ks);
delk_delz = delK_delTheta.*deltheta_delz; % TERM 3
Dh_DTheta = delh_delTheta_ALT(theta, thta, [alpha, N]);

% TERM 1
TERM1 = delk_delz.*Dh_DTheta.*deltheta_delz;

% TERM 2
K = K_vanGenuchten1980(theta, thta, M, Ks);
D2h_DzDTheta = D2h_DzDThetaALT(deltheta_delz, theta, thta, [alpha, N]);
TERM2 = K.*(D2h_DzDTheta.*deltheta_delz + Dh_DTheta.*del2theta_delz2);

% TERM 3 
TERM3 = delk_delz;

% -------------- LHS and RHS ----------------
RHS = TERM1 + TERM2 + TERM3;
LHS = deltheta_deltime;

% Trimming ~~~~

fraction_of_Zdata_removed = DataTrim.z; % Value in fraction in one Side
fraction_of_tdata_removed = DataTrim.t; % Value in fraction in one Side

% Twice of fraction_of_Zdata_removed will be removed

[ZZZ, TTT] = size(theta);
zdataselection = (ceil(fraction_of_Zdata_removed*ZZZ)+1):...
                (floor((1-fraction_of_Zdata_removed)*ZZZ));

tdataselection = (ceil(fraction_of_tdata_removed*TTT)+1):...
                (floor((1-fraction_of_tdata_removed)*TTT));
    
% ############## Defining Error ####################

errMatrix = LHS - RHS;

errNorm1 = (sum(sum(abs(errMatrix))));
errNorm2 = sqrt(sum(sum((errMatrix(zdataselection,tdataselection).^2)./...
                        VarianceMatrix(zdataselection,tdataselection))));
errNormInf = max(max(abs(errMatrix)));

