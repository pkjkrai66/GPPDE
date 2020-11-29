function y = delh_delTheta_ALT(theta, thta, PDEparams)

THETA = (theta-thta.r)/(thta.s-thta.r);
alpha  = PDEparams(1);
n = PDEparams(2);
m = 1-1/n;

h = -1;
C0 = -alpha*m*n*sign(h)*(thta.s-thta.r);

y = (1/C0)*(THETA.^(-1/m)).*(1-THETA.^(1/m)).^(-m);