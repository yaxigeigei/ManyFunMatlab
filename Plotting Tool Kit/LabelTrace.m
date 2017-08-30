function LabelTrace( yArray, color )
%PLOTWITHLABEL Summary of this function goes here
%   Detailed explanation goes here


if nargin < 2
    color = 'k';
end

yArray = yArray';

for i = 1 : size(yArray, 1)
    lastDataPointX = find(~isnan(yArray(i,:)), 1, 'last');
    lastDataPointY = yArray(i, lastDataPointX);
    text(lastDataPointX, lastDataPointY, num2str(i), 'color', color);
end




end

