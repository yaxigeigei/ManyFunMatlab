function [ img, bgValues ] = dFluorBackgroundMinus( img, bgROI )
%Subtract background of an image or image stack frame by frame according to the specified background region
%
% [ img, bgValues ] = dFluorDebackground( img, bgROI )
%
%
%Input discription
%
% img:
% An two or three dimensional image array
%
% bgROI:
% A 1*2 cell array. The first cell contains coordinates of its contour points and the second cell is its logic mask.
% (However, this function only uses the information from the second cell)
%
%
%Output discription
%
% img:
% The input image whose background was subtracted
%
% bgValues:
% A vector containing the background values of each respective frame


    imgFrames = size(img, 3);
    
    bgValues = zeros(imgFrames, 1);

    for i = 1 : imgFrames
        oneFrame = img(:,:,i);
        bgValues(i) = mean(oneFrame(find(bgROI{2})));
        img(:,:,i) = img(:,:,i) - bgValues(i);
    end



end

