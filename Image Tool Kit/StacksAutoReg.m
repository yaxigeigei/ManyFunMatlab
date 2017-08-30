function [ stacks, accumVects ] = StacksAutoReg( stacks, varargin )
%Automatic registration of a three dimensional image array
%
% [ img ] = StackRegister( img, template, multiCores, optimizer, metric )
%
%
%Input discription
%
% img (required): 
% A three dimensional array to be registered
% 
% template:
% The index of frame taken as reference. (default is the frame at the middle
% of the entire image stack)
%
% multiCore:
% true      use multi-core processing (suitable for large stack, e.g. frame number > 50)
% false     use single core (default, faster startup)
%
% optimizer, metric:
% Configurations generated by MATLAB function imregconfig() and customized
% for best registration.
%
%
%Output discription
%
% img: the registered array
    

    p = inputParser;
    
    addParamValue(p, 'ref', [ ], @isnumeric);
    addParamValue(p, 'quality', 1, @isnumeric);
    addParamValue(p, 'size', 'same', @(x) any(validatestring(x, { 'same', 'cut' })));
    addParamValue(p, 'multiCore', false, @islogical);
    
    parse(p, varargin{:});
    
    ref = p.Results.ref;
    fineFactor = p.Results.quality;
    resultSize = p.Results.size;
    multiCore = p.Results.multiCore;
    
    
    stacks = StacksArray2cell(stacks);
    if length(stacks) > 1
        maxStack = StacksProjMax(stacks);
    else
        maxStack = stacks{1};
    end
    maxStack = im2double(maxStack);
    
    
    if isempty(ref)
        ref = ceil(size(maxStack,3)/2);
    end
    disp(ref);
    
    if multiCore
        matlabpool
    end
    
    tic
    
    accumVects = zeros(size(maxStack,3), 2);
    imgRef = maxStack(:,:,ref);
    parfor k = 1 : size(maxStack,3)
        accumVects(k,:) = ImageShift(maxStack(:,:,k), imgRef, 'fineFactor', fineFactor);
        disp(k);
    end
    
    toc

    if multiCore
        matlabpool close
    end
    
    accumVects = round(accumVects);
    
    if length(stacks) > 1
        stacks = StacksCrop(stacks, resultSize, accumVects);
    else
        stacks = StackCrop(stacks{1}, resultSize, accumVects);
    end
    
    



end

