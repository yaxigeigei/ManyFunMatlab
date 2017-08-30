function [ snROI, roiFMasks, roiCMasks ] = dFluorAutoROI( gridFRatio, edgePixels, estiROI, maxIteration )
%This function is still under construction
    
    

    % Find Best Threshold by Iteration
    maxGridFRatio = max(gridFRatio, [ ], 3);
    initThreshold = max(max(maxGridFRatio));
    thresholdGradient = initThreshold * 0.95.^(1:maxIteration);
    numPatchVect = zeros(size(thresholdGradient));
    evoMaps = zeros([ size(maxGridFRatio), maxIteration ]);
    
    for i = 1 : length(thresholdGradient)
        roiMap = maxGridFRatio;
        roiMap(roiMap < thresholdGradient(i)) = 0;
        roiMap(roiMap ~= 0) = 1;
        evoMaps(:,:,i) = roiMap;
        
        identification = bwconncomp(roiMap, 4);
        for j = 1 : identification.NumObjects
            if length(identification.PixelIdxList{j}) >= edgePixels^2
                numPatchVect(i) = numPatchVect(i) + 1;
            end
        end
    end
    
    distanceVect = abs(numPatchVect - estiROI);
    [ ~, index ] = min(distanceVect);
    bestThreshold = thresholdGradient(index);
    
    
    
    % Use Best Threshold to Identify Patches
    roiMap = maxGridFRatio;
    roiMap(roiMap < bestThreshold) = 0;
    roiMap(roiMap ~= 0) = 1;
    
    identification = bwconncomp(roiMap, 4);
    for i = 1 : identification.NumObjects
        if length(identification.PixelIdxList{i}) < edgePixels^2
            roiMap(identification.PixelIdxList{i}) = 0;
        end
    end
    
    [ roiFMasks, roiCMasks ] = GeoPatchSplit(roiMap);
    
    
    % Split Joint Patches and Remove Non-correlated Patches
    roiMapsCell = { };
    
    for i = 1 : length(roiCMasks)
        temp = cat(3, zeros(size(roiCMasks{i})), roiCMasks{i});
        while size(temp, 3) == 2
            temp = dFluorRoiCorrSplit(temp(:,:,2), gridFRatio, edgePixels);
            if ~isempty(temp)
                roiMapsCell{end+1} = temp(:,:,1);
            end
        end
    end
    
    
    
    % Refine and Output Fully Split ROIs
    snROI = cell(length(roiMapsCell), 2);
    roiCMasks = zeros([ size(roiMapsCell{1}), length(roiMapsCell) ]);
    
    for i = 1 : length(roiMapsCell)
        roiMapsCell{i} = imfill(roiMapsCell{i}, 'holes');
        roiMapsCell{i} = bwmorph(roiMapsCell{i}, 'erode', ceil(edgePixels/10));
        roiMapsCell{i} = bwmorph(roiMapsCell{i}, 'dilate', ceil(edgePixels/10));
        identification = bwconncomp(roiMapsCell{i}, 4);
        for j = 1 : identification.NumObjects
            if length(identification.PixelIdxList{j}) < edgePixels
                roiMapsCell{i}(identification.PixelIdxList{j}) = 0;
            end
        end
        
        roiCMasks(:,:,i) = roiMapsCell{i};
        roiContours = find(bwmorph(roiMapsCell{i}, 'remove'));
        [ y, x ] = ind2sub(size(roiMapsCell{i}), roiContours);
        snROI(i, [1 2]) = { [ x, y ], roiCMasks(:,:,i) };
    end
    
    
    

end










