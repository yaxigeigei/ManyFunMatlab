function [ imgHeight, imgWidth, imgSlices ] = GetTiffImageSize( tifObj )
%GETTIFFOBJFRAMENUM Summary of this function goes here
%   Detailed explanation goes here

% Convert input to a Tiff object
if ischar(tifObj)
    tifObj = Tiff(tifObj, 'r');
end

% Get the height and width
imgWidth = tifObj.getTag('ImageWidth');
imgHeight = tifObj.getTag('ImageLength');

% Determine the search range
probeIncrement = 1000;
upperLim = 1;
while true
    try
        tifObj.setDirectory(upperLim);
        lowerLim = upperLim;
        upperLim = upperLim + probeIncrement;
    catch
        break;
    end
end

% Dichotomizing the search range
while upperLim - lowerLim > 1
    midDir = floor((upperLim + lowerLim)/2);
    try
        tifObj.setDirectory(midDir);
        lowerLim = midDir;
    catch
        upperLim = midDir;
    end
end

try
    tifObj.setDirectory(upperLim);
    imgSlices = upperLim;
catch
    imgSlices = lowerLim;
end


end

