function [ img, masks ] = dFluorBackgroundEven( img, radius )
%DFLUOREVENBACKGROUND Summary of this function goes here
%   Detailed explanation goes here
    
    
    if nargin < 2
        radius = 10;
    end
    
    
    se = strel('disk', radius);
    h = fspecial('gaussian', radius*2, radius/3);
    
    masks = zeros(size(img));
    for i = 1 : size(img, 3)
        masks(:,:,i) = imopen(img(:,:,i), se);
        masks(:,:,i) = imfilter(masks(:,:,i), h);
        img(:,:,i) = img(:,:,i) - masks(:,:,i);
    end
    
    img(img < 0) = 0;
    img = img / (max(max(max(img))) - min(min(min(img))));
    masks = masks / (max(max(max(img))) - min(min(min(img))));
    
    
    
end