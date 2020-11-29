function hndl = myHline(y, varargin)
% plots a verticle line at given x, were x can be a vector
% varargin takes all values which default function line() takes
% line() properties could be looked at
% http://in.mathworks.com/help/matlab/ref/primitiveline-properties.html

% Ex:
% plot([1 2], [3,4])
% hndl = myHline([3.4 3.5], 'Color', 'r', 'Marker', '*', ...
%           'LineStyle', ':','MarkerSize',10)

HLIM = get(gca, 'Xlim');
hndl = zeros(length(y), 1);
for i = 1:length(y)
    hndl(i) = line(HLIM, [y(i) y(i)], varargin{:});
end
