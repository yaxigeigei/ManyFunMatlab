function [ stacks, accumVects, sideStacks ] = StacksAutoRegZ( stacks, zScale )
%STACKSAUTOREGZ Summary of this function goes here
%   Detailed explanation goes here



    if nargin < 2
        zScale = 10;
    end
    
    stacks = StacksArray2cell(stacks);
    
    
    % Pad Stacks with Unequal Slices
    
    maxNumSlices = 1;
    for i = 1 : length(stacks)
        maxNumSlices = max(maxNumSlices, size(stacks{i}, 3));
    end
    
    for i = 1 : length(stacks)
        [ imgHeight, imgWidth, imgSlices ] = size(stacks{i});
        stacks{i} = cat(3, stacks{i}, zeros(imgHeight, imgWidth, maxNumSlices - imgSlices));
    end
    
    
    % Registration
    
    sideStacks = cell(length(stacks), 1);
    for i = 1 : length(stacks)
        sideStacks{i} = permute(stacks{i}, [ 3, 1, 2 ]); % Create Side-viewed Stacks
    end
    maxSideStack = StacksProjMax(sideStacks);
    maxSideStack = im2double(maxSideStack);
    maxSideStack = StackInterp(maxSideStack, 1, zScale);
    
    accumVects = zeros(length(stacks), 1);
    imgRef = maxSideStack(:,:,1);
    msSize = [ size(maxSideStack,1), ceil(size(maxSideStack,2)/zScale/10) ];
    for k = 2 : size(maxSideStack, 3)
        ms = XMeanSquare2(maxSideStack(:,:,k), imgRef, msSize);
        [ ~, index ] = min(ms(:));
        [ i, ~ ] = ind2sub(msSize, index);
        accumVects(k) = i - size(ms,1)/2;
    end
    
    accumVects = round(accumVects / zScale);
    stacks = StacksCrop(stacks, 'zaxis', accumVects);
    
    for i = 1 : length(stacks)
        sideStacks{i} = permute(stacks{i}, [ 3, 1, 2 ]); % Create Side-viewed Stacks
    end


end

