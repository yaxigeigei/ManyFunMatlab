function [ roiMaps ] = dFluorRoiCorrSplit( roiMap, gridFRatio, edgePixels )
%DFLUORROICORRSPLIT Summary of this function goes here
%   Detailed explanation goes here
    
    
    linearIndices = find(roiMap);
    [ rowIndices, columnIndices ] = find(roiMap);
    
    linearFRatio = zeros(length(linearIndices), size(gridFRatio, 3));
    for i = 1 : length(linearIndices)
        linearFRatio(i,:) = gridFRatio(rowIndices(i), columnIndices(i), :);
    end
    
    
    corrMap = corrcoef(linearFRatio');
    [ ~, maxCorrIndex ] = max(sum(corrMap(:,1:round(edgePixels^2 * 0.75))));
    [ ~, order ] = sort(corrMap(maxCorrIndex,:), 'descend');
    corrMapSort = corrcoef(linearFRatio(order,:)');
    
    
    warning off
    templateVect = mean(corrMapSort(1:round(edgePixels^2 * 0.75), :), 1);
    fitObject1 = fit((1:length(linearIndices))', templateVect', 'poly1');
    fitObject5 = fit((1:length(linearIndices))', templateVect', 'poly5');
    fitVectX = 1 : length(linearIndices);
    fitVect5 = fitObject5.p1*fitVectX.^5 + fitObject5.p2*fitVectX.^4 + fitObject5.p3*fitVectX.^3 ...
        + fitObject5.p4*fitVectX.^2 + fitObject5.p5*fitVectX + fitObject5.p6;
    
    
    maxCorrValue = min(fitObject1(1), fitObject5(1));
    if maxCorrValue < 0.4
        roiMaps = [ ];
    else
        cutoffValue = maxCorrValue / 2;
        cutoffIndex = find(fitVect5 < cutoffValue, 1);
        extractIndices = linearIndices(order(1:cutoffIndex));
        if cutoffIndex ~= length(linearIndices)
            leftoverIndices = linearIndices(order((cutoffIndex+1):end));
        else
            leftoverIndices = 0;
        end
        
        if length(leftoverIndices) >= edgePixels^2 * 0.75
            subMap1 = zeros(size(roiMap));
            subMap1(extractIndices) = 1;
            subMap2 = zeros(size(roiMap));
            subMap2(leftoverIndices) = 1;
%             dFluorPlotIdentityMap(subMap1);
%             dFluorPlotIdentityMap(subMap2);
            roiMaps = cat(3, subMap1, subMap2);
        else
            roiMaps = roiMap;
        end
    end
    
    
%     if maxCorrValue < 0.4 || length(leftoverIndices) >= edgePixels^2
%         
%         figure
%         subplot(1,2,1)
%         pcolor(corrMapSort);
%         axis equal tight ij
%         shading flat
%         
%         subplot(1,2,2)
%         hold on
%         plot(templateVect);
%         plot(fitObject1);
%         plot(fitObject5);
%         xlim([ 1 length(linearIndices) ]);
%         axis square
%         hold off
%         
%     end
    
    
end
