function delK = delK_delTheta_fofTheta(theta, M, thta, Ks)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ------------------ Ks Incorporated -------------------
%                    ~~~~~~~~~~~~~~~

% Verified
% 24-10-2016, 25-10-2016
% The function has been verified atleast three times and the Function is
% written correctly, there is NO possibility of the function being wrong 

THETA = (theta-thta.r)/(thta.s-thta.r);
theta_diff = thta.s-thta.r;

Term1 = 1 - THETA.^(1/M);

Term2 = (1-Term1.^M).^2/2 ;
Term3 = 2*(THETA.^(1/M)).*(1-Term1.^M).*(Term1).^(M-1);
delK = Ks*THETA.^(-1/2).*(Term2 + Term3)/theta_diff;

