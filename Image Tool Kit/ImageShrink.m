function [ imgShrinked ] = ImageShrink( img, factor )
%Change the resolution of an image or image stack
%
% imgShrinked = ImageShrink(img, factor)
%
%
%Input discription
%
% img:
% An two or three dimensional image array
%
% factor:
% One pixel in shrinked image comes from factor*factor pixels of the original image
%
%
%Output discription
%
% imgShrinked:
% Image or image stack with shrinked resolution


    factor = round(factor);
    
    [ imgHeight, imgWidth, imgFrames ] = size(img);
    
    img = img(1 : imgHeight - mod(imgHeight,factor), 1 : imgWidth - mod(imgWidth,factor), :);
    [imgHeight, imgWidth, ~] = size(img);
    
    shrinkedWidth = imgWidth / factor;
    shrinkedHeight = imgHeight / factor;
    
    img = reshape(img, [factor, shrinkedHeight, factor, shrinkedWidth, imgFrames]);
    imgShrinked = squeeze(mean(mean(img, 1), 3));



end

