function ddC = dxdzRBFcov(s,t, kernelParam, DiffWRT_Dimention)

%% Calculates auto-covariance
% GP derivative d^2{c(s,t)}/dsdt
% s and t = [x1, x2 ...]
%   Where x1, x2 ... are indipendent variables arranged as column vectors
%   If s and t are 1-Dimensional then they must only be Column Vectors
% ppar = Hyperparameters of Gaussian Process
% More information available in function `GPcov()`
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
                    % dist = zeros(n, m);
                    % % --Calculating Euclidean distance ---------
                    % for i = 1:D
                    %     x1Matrix = repmat(x1(:,i), 1, m);
                    %     x2Matrix = repmat(x2(:,i)', n, 1);
                    %     temp_dist = (x1Matrix - x2Matrix).^2;
                    %     dist = dist + temp_dist;
                    % end
                    % expterm = exp(-ppar(2)*dist);
                    % C = (2*ppar(2)-4*ppar(2)^2*dist).*GPcov(x1,x2,ppar);

%% ---------------------------------------------------------------
% ----------------------------------------------------------------
ddC = zeros(n, m);
BasisVect = zeros(1, D);                % Basis Vector
BasisVect(1,DiffWRT_Dimention) = 1;

l = kernelParam(1);    % Characteristic length scale
f = kernelParam(2);    % Controls the vertical variation

C = GPcov(x1, x2, kernelParam);
for i = 1:n
   xRwVect = x1(i,:)';
   for j = 1:m
      zRwVect = x2(j,:)';
        
      Y = xRwVect - zRwVect;
      ddC(i,j) = (C(i,j)/l^2)*(1-(BasisVect*Y)^2/(l^2));
    end
end

end

