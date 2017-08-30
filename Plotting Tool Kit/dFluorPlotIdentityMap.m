function dFluorPlotIdentityMap( roiMaps )
% This function is going to be merged with others in the future
    

    roiMap = max(roiMaps, [ ], 3);
    
    figure('Color', 'w');
    pcolor(roiMap(end:-1:1, 1:end));
    axis equal tight
    colormap gray
    
    
    
end

