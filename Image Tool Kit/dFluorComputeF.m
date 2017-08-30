function [ actualF, baselineF, fRatio ] = dFluorComputeF( img, mask, baselineWindow )
%Compute actual, baseline fluorescence intensities and deltaF/F of the input image series or its masked subregion
%
% [ actualF, baselineF, fRatio ] = dFluorComputeF( img, mask, baselineWindow )
%
%
%Input discription
%
% img: 
% An three dimensional array in double
%
% mask: 
% A logic array the same size of each frame where pixels included are 1 (default is [ ] for whole image)
% 
% baselineWindow:
% Range of frames taken as baseline fluorescence intensity. (e.g. [ 1 10 ])
% (default is 'auto' under which the program takes the average of lowest 15% intensity as baseline)
%
% 
%Output discription
%
%Actual, baseline fluorescence intensities and deltaF/F


[ imgHeight, imgWidth, imgFrames ] = size(img);


if isempty(mask)
    actualF = squeeze(sum(sum(img,1),2));
else
    roi = zeros(imgHeight, imgWidth, imgFrames);
    for i = 1 : imgFrames
        roi(:,:,i) = img(:,:,i) .* mask;
    end
    pixelNum = sum(sum(mask));
    actualF = squeeze(sum(sum(roi,1),2) / pixelNum);              % Actual F as means of signal ROI
end


if ~strcmp(baselineWindow, 'auto')
    baselineF = mean(actualF(baselineWindow(1):baselineWindow(2)));
    baselineF = ones(imgFrames,1)*baselineF;
else
    sortedActF = sort(actualF);
    lowActF = sortedActF(1 : ceil(length(sortedActF)*0.15));     % Lowest 15% values selected as baseline
    baselineF = ones(imgFrames,1)*mean(lowActF);               % The length of baselineF is n
end


deltaF = actualF - baselineF;
fRatio = deltaF./baselineF;
fRatio(isinf(fRatio)) = 0;





end

