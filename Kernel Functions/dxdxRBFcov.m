function ddC = dxdxRBFcov(s, t, kernelParam, DiffWRT_Dimention)
%% dxdxGPcov gives the Double derivative of a Gaussian Process
%  the derivative is wrt to second argument i.e. wrt `t` in this function
%    s and t = [x1, x2 ...]
%    % GP derivative d^2{c(s,t)}/dtdt
%    Where x1, x2 ... are indipendent Variables/Dimension arranged as column vectors
%    If s and t are 1-Dimensional then they must only be Column Vectors
% ppar = Hyperparameters of Gaussian Process
% DiffWRT_Dimention = Dimension/Variable with respect to which double
%                     derivative is required
% More information available in function `RBFcov()`
%%

x1 = s;
x2 = t;

[n, D] = size(x1);
[m, d] = size(x2);
if size(x1,2) ~= size(x2,2) % D ~= d
    error('Error: Dimension mismatch of x1 and x2')
end

if (nargin < 4) && (D == 1)
    DiffWRT_Dimention = 1;
end
% ---------------------------------------------------------------
ddC = zeros(n, m);
ddC2 = zeros(n, m);
BasisVect = zeros(1, D);                % Basis Vector
BasisVect(1,DiffWRT_Dimention) = 1;

l = kernelParam(1);    % Characteristic length scale
f = kernelParam(2);    % Controls the vertical variation

C = RBFcov(x1, x2,kernelParam);
for i = 1:n
   xRwVect = x1(i,:)';
    for j = 1:m
      zRwVect = x2(j,:)';
      Y = xRwVect - zRwVect;
      ddC(i,j) = (C(i,j)/l^2)*((BasisVect*Y)^2/(l^2) - 1);
    end
end

end

