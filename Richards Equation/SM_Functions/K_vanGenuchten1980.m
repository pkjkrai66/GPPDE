function hydrauConductivity = K_vanGenuchten1980(theta, thta, m,Ks)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ------------------ Ks Incorporated -------------------
%                    ~~~~~~~~~~~~~~~

% Verified
% Hydraulic conductivity as defined by vanGenuchten 1980

THETA = (theta-thta.r)/(thta.s-thta.r);
hydrauConductivity = Ks*THETA.^(1/2).*(1-(1-THETA.^(1/m)).^m).^2;
