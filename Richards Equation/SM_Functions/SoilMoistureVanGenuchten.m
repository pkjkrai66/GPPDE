function theta = SoilMoistureVanGenuchten(h, thta, Param)
    
alpha = Param.alpha;   % 1/cm
N = Param.N;          % Dimention Less
M = 1-1/N;

theta = zeros(size(h));
tempI = h < 0;
    theta(tempI) = thta.r + (thta.s-thta.r)*(1./(1+(alpha*abs(h(tempI))).^N)).^M;
    theta(~tempI) = thta.s;
end