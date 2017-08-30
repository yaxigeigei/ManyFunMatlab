function [ stackArray, stackCells ] = GeoPatchSplit( img )
%GEOPATCHSPLIT Summary of this function goes here
%   Detailed explanation goes here


identification = bwconncomp(img, 4);

stackArray = zeros([ size(img), identification.NumObjects ]);
stackCells = cell(identification.NumObjects, 1);
framePixels = size(img, 1) * size(img, 2);

for i = 1 : identification.NumObjects
    stackArray((i-1)*framePixels + identification.PixelIdxList{i}) = 1;
    stackCells{i} = stackArray(:,:,i);
end



end

