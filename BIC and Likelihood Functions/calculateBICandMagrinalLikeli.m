function [BIC, logMarginalLiklihood, VID] = calculateBICandMagrinalLikeli(GP, X, Y, meany)
    
% The function calculates Bayesian Information Criteria (BIC)and 
% log marginal likelihood for GP for given data and hyperparameters
    

    VID = sortrows(combvec(GP.hypPara{:}).');
    I_NUM = size(VID,1);

    logLikli = zeros(I_NUM, 1);    
    DOF = zeros(I_NUM, 1); 
    BIC = zeros(I_NUM, 1); 
    logMarginalLiklihood = zeros(I_NUM, 1);
    N = size(X, 1);

for i = 1:I_NUM
    disp(i)
    gppar_temp = VID(i, 1:2);
    gpnoise_temp = VID(i, 3);

    C = GP.CovFun(X, X, gppar_temp);
    GramMatrix = C + gpnoise_temp^2*eye(N);
    
    % ---------------------- BIC -------------------------------
    logLikli(i) = logmvnpdf(Y', meany', GramMatrix);
    DOF(i) = trace((GramMatrix)*(GramMatrix)');
    BIC(i)  = -2*logLikli(i) + log(N)*DOF(i);
    
    % --------------------- Log Marginal Liklihood ----------------------------
    logDeterminant = logdet(GramMatrix, 'chol');
%     logMarginalLiklihood(i) = -Y'*(GramMatrix\Y)/2 - logDeterminant/2 - TZ*log(2*pi)/2;   
    
    logMarginalLiklihood(i) = -(Y-meany)'*(GramMatrix\(Y-meany))/2 - logDeterminant/2 - N*log(2*pi)/2; 
end
