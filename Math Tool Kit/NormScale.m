function [ normPool ] = NormScale( pool, colRef )
%NORMALIZATION Summary of this function goes here
%   Detailed explanation goes here


    matRef = repmat(colRef, [1 size(pool, 2)]);
    normPool = pool ./ matRef;


end

