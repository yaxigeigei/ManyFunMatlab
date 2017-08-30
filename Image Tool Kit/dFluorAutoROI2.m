function [ snROI, roiFMasks, roiCMasks ] = dFluorAutoROI2( gridFRatio, edgePixels, estiROI, maxIteration )
%This function is still under construction
    
    

    maxGridFRatio = max(gridFRatio, [ ], 3);
    initThreshold = max(max(maxGridFRatio));
    thresholdGradient = initThreshold * 0.95.^(1:maxIteration);
    evoMaps = zeros([ size(maxGridFRatio), maxIteration ]);
    
    for i = 1 : length(thresholdGradient)
        
        roiMap = maxGridFRatio;
        roiMap(roiMap < thresholdGradient(i)) = 0;
        roiMap(roiMap ~= 0) = 1;
        
        identification = bwconncomp(roiMap, 4);
        for j = 1 : identification.NumObjects
            if length(identification.PixelIdxList{i}) < edgePixels^2
                roiMap(identification.PixelIdxList{i}) = 0;
            end
        end
        
        identification = bwconncomp(roiMap, 4);
        if identification.NumObjects > 0
            [ ~, roiCMasks ] = GeoPatchSplit(roiMap);
            
            roiMapsCell = { };
            for j = 1 : length(roiCMasks)
                temp = cat(3, zeros(size(roiCMasks{j})), roiCMasks{j});
                while size(temp, 3) == 2
                    temp = dFluorRoiCorrSplit(temp(:,:,2), gridFRatio, edgePixels);
                    if ~isempty(temp)
                        roiMapsCell{end+1} = temp(:,:,1);
                    end
                end
            end
            
            if ~isempty(roiMapsCell)
                roiCMasks = zeros([ size(roiMapsCell{1}), length(roiMapsCell) ]);
                for j = 1 : length(roiMapsCell)
                    roiMapsCell{j} = imfill(roiMapsCell{j}, 'holes');
                    roiMapsCell{j} = bwmorph(roiMapsCell{j}, 'erode', ceil(edgePixels/10));
                    roiMapsCell{j} = bwmorph(roiMapsCell{j}, 'dilate', ceil(edgePixels/10));
                    identification = bwconncomp(roiMapsCell{j}, 4);
                    for k = 1 : identification.NumObjects
                        if length(identification.PixelIdxList{k}) < edgePixels
                            roiMapsCell{j}(identification.PixelIdxList{k}) = 0;
                        end
                    end
                end
                roiMap = max(roiCMasks, [ ], 3);
            end
        end

        evoMaps(:,:,i) = roiMap;
        
    end
    
    
    

end










