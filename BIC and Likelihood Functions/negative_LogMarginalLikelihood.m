function [negative_logMrgnlLklihd, grad_nLogMargLikeli] = ...
                    negative_LogMarginalLikelihood(hyperpara,X,Y,meany)

% The function calculates negative log marginal-likelihood for GP and 
% gradient of `negative log marginal-likelihood` wrt hyperparameters; 
% given data and hyperparameters
NaNindices = isnan(Y);
Y = Y(~NaNindices);
X = X(~NaNindices, :);

kernelParam = hyperpara(1:2);
gpnoise = hyperpara(3);
C = RBFcov(X, X, kernelParam);

N = size(X, 1);
GramMatrix = C + gpnoise^2*eye(N);

% --------------------- Log Marginal Liklihood ----------------------------
logDeterminant = logdet(GramMatrix, 'chol');
logMarginalLiklihood = -(Y-meany)'*(GramMatrix\(Y-meany))/2 - logDeterminant/2 - N*log(2*pi)/2;

% ~~~~~~~~~ Taking negative of logMarginalLikelihood ~~~~~~~~~~~
negative_logMrgnlLklihd = -logMarginalLiklihood;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% ------ Derivative of log marginal-likelihood -----------------------
if nargout > 1
inv_GramMat = inv(GramMatrix);      %#ok<*MINV>
[dk_dl, dk_df] = delk_delhyper(X, kernelParam);

grad_LogMargLikeli(1) = (1/2)*(Y-meany)'*inv_GramMat*dk_dl*inv_GramMat*(Y-meany)...
                    -(1/2)*trace(inv_GramMat*dk_dl); 
grad_LogMargLikeli(2) = (1/2)*(Y-meany)'*inv_GramMat*dk_df*inv_GramMat*(Y-meany)...
                    -(1/2)*trace(inv_GramMat*dk_df); 

grad_LogMargLikeli_dsigma = (1/2)*(Y-meany)'*inv_GramMat*(2*gpnoise*eye(N))...
            *inv_GramMat*(Y-meany)-(1/2)*trace(inv_GramMat*2*gpnoise*eye(N));

% ~~~~~~~ Taking negative of gradients ~~~~~~~~~~~
grad_nLogMargLikeli = [-grad_LogMargLikeli, -grad_LogMargLikeli_dsigma];
end
% -X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X

%%  
% The function provides derivative of RBF covariance function with respect
%  to the hyperparameters
function [dk_dl, dk_df] = delk_delhyper(X, kernelParam)
    x1 = X;
    x2 = X;

    [n, D] = size(x1);
    distSquare = zeros(n, n);
    % --Calculating Euclidean distance ---------
    for i = 1:D
        x1Matrix = repmat(x1(:,i), 1, n);
        x2Matrix = repmat(x2(:,i)', n, 1);
        temp_dist = (x1Matrix - x2Matrix).^2;
        distSquare = distSquare + temp_dist;
    end

    l = kernelParam(1);    % Characteristic length scale
    f = kernelParam(2);    % Controls the vertical variation

    dk_dl = f^2*exp(-distSquare/(2*l^2)).*distSquare/l^3;
    dk_df = 2*f*exp(-distSquare/(2*l^2));

end   
end