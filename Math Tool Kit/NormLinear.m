function [ result, correction ] = NormLinear( data, baselineMask )
%CORRECTBLEACHING Summary of this function goes here
%   Detailed explanation goes here


    if nargin < 2
        baselineMask = 1 : size(data, 2);
    end


    correction = zeros(size(data));
    fitVectX = 1 : size(data, 2);

    for i = 1 : size(data, 1)
        fitObject = fit((1:length(baselineMask))', data(i,baselineMask)', 'poly1');
        correction(i,:) = fitObject.p1 * (fitVectX - 1);
    end

    result = data - correction;



end

