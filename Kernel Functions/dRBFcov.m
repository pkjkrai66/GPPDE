function dC = dRBFcov(s,t, kernelParam,DiffWRT_Dimention)
%% dGPcov gives the first order derivative of a Gaussian Kernel ...
%  the derivative is wrt to second argument i.e. wrt `t` in this function
%  GP derivative d{c(s,t)}/dt
%    s and t = [x1, x2 ...]
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

if (nargin < 4) && (D ==1)
    DiffWRT_Dimention = 1;
end
% ---------------------------------------------------------------
dC = zeros(n, m);
BasisVect = zeros(1, D);                % Basis Vector
BasisVect(1,DiffWRT_Dimention) = 1;
C = RBFcov(x1, x2, kernelParam);

l = kernelParam(1);    % Characteristic length scale
f = kernelParam(2);    % Controls the vertical variation

for i = 1:n
   xRwVect = x1(i,:)';
    for j = 1:m
      zRwVect = x2(j,:)';
      Y = xRwVect - zRwVect;
      dC(i,j) = (1/l^2)*C(i,j)*BasisVect*Y;
    end
end

end