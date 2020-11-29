function y = D2h_DzDThetaALT(deltheta_delz, theta, thta, PDEparams)

THETA = (theta-thta.r)/(thta.s-thta.r);
theta_diff = thta.s-thta.r;

alpha  = PDEparams(1);
n = PDEparams(2);
m = 1-1/n;

h = -1;
C0 = -alpha*m*n*sign(h)*(thta.s-thta.r);

C1 = 1-THETA.^(1/m);
y = (1/(C0*(thta.s-thta.r)))*deltheta_delz.*(...
      (-1/m)*THETA.^((-1-m)/m).*C1.^(-m) + THETA.^(-1).*(C1.^(-m-1))...
      );

