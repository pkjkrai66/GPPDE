function hndl = myVline(x, varargin)
% plots a verticle line at given x, were x can be a vector
% varargin takes all values which default function line() takes
% line() properties could be looked at
% http://in.mathworks.com/help/matlab/ref/primitiveline-properties.html

% Ex:
% plot([1 2], [3,4])
% hndl = myVline([1.4 1.5], 'Color', 'b', 'Marker', '*', ...
%           'LineStyle', ':','MarkerSize',10)

YLIM = get(gca, 'Ylim');
hndl = zeros(length(x), 1);
for i = 1:length(x)
    hndl(i) = line([x(i) x(i)], YLIM, varargin{:});
end
set(gca, 'Ylim', YLIM)