function [ img, k ] = ImageSubtract( img, imgSubtract )
%Background subtraction based on the image that only contains background
%signal. For example, you can subtract red channel from green channel in
%order to eliminate inspecific fluorescence. The intensity of background in
%each channel are not necessarily the same.
%
% img = ImageSubtract(img, imgSubtract)
%
%
%Input discription
%
% img: 
% A two or three dimensional image array containing both signal and background.
% 
% imgSubtract: 
% The corresponding image array representing the background.
% 
%
%Output discription
%
% img:
% The resulting image array.
    
    
    [ imgHeight, imgWidth, imgFrames ] = size(img);
    
    k = zeros(1,imgFrames);
    
    
    for i = 1 : imgFrames
    
        mask = zeros(imgHeight, imgWidth);
        serialImgSubtract = double(reshape(imgSubtract(:,:,i), imgHeight * imgWidth, 1));
        mask(imgSubtract(:,:,i) > median(serialImgSubtract)) = 1;
        
        maskedImg = double(img(:,:,i)) .* mask;
        maskedImgSubtract = double(imgSubtract(:,:,i)) .* mask;
        
        k(i) = sum(maskedImg(:)) / sum(maskedImgSubtract(:));
        img(:,:,i) = img(:,:,i) - k(i) * imgSubtract(:,:,i);
        
    end
    
    img(img < 0) = 0;
    
    
    
end

