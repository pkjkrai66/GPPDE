function C = RBFcov(s, t, kernelParam)
% s and t = [x1, x2 ...]
%   Where x1, x2 ... are indipendent variables arranged as column vectors
%   If s and t are 1-Dimensional then they must only be Column Vectors
%   The data should be arranged in this order
%              ----------------             
%              | 5i | 7j | 9k |             [ 5i, 7j, 9k ; ...
%              | 2i | 8j | 4k |     =         2i, 8j, 9k ; ...
%     x =      | 7i | 6j | 1k |               7i, 6j, 1k ; ...
%              | 9i | 3j | 8k |               9i, 3j, 8k ]
%              ----------------

%%

x1 = s;
x2 = t;

[n, D] = size(x1);
[m, d] = size(x2);
if size(x1,2) ~= size(x2,2) % D ~= d
    error('Error: Dimension mismatch of x1 and x2')
end
distSquared = zeros(n, m);
% --Calculating Euclidean distance ---------
for i = 1:D
    x1Matrix = repmat(x1(:,i), 1, m);
    x2Matrix = repmat(x2(:,i)', n, 1);
    temp_dist = (x1Matrix - x2Matrix).^2;
    distSquared = distSquared + temp_dist;
end

l = kernelParam(1);    % Characteristic length scale
f = kernelParam(2);    % Controls the vertical variation

C = f^2*exp(-distSquared/(2*l^2));

end